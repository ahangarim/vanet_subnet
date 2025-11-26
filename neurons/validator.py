"""
VANET Validator - Genomics Variant Calling Subnet

Sends BAM files with synthetic mutations to miners,
collects VCF results, and scores using hap.py validation.
"""

import sys
import os
from pathlib import Path

# Add parent directory to path so we can import base and utils
sys.path.insert(0, str(Path(__file__).parent.parent))

import time
import typing
import bittensor as bt
import argparse
import torch
import random
import asyncio
import json
import tempfile
from datetime import datetime
from dotenv import load_dotenv

from base import TaskSynapse, ComputeSynapse, GenomicsTaskSynapse
from base import GENOMICS_CONFIG, VALIDATOR_CONFIG, get_manifest_path, is_docker_available, BASE_DIR
from utils import (
    GenomicsTaskGenerator,
    MutationInjector,
    HappyScorer,
    ScoreTracker,
    AdvancedScorer
)

sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)
load_dotenv()


class Validator:
    """VANET validator for genomics variant calling tasks."""

    def __init__(self, config=None):
        self.config = config or self.get_config()

        bt.logging.info("Setting up validator...")
        bt.logging.set_trace(self.config.logging.trace)
        bt.logging.set_debug(self.config.logging.debug)

        bt.logging.info(f"Loading wallet: {self.config.wallet.name}/{self.config.wallet.hotkey}")
        self.wallet = bt.wallet(config=self.config)
        bt.logging.info(f"Wallet loaded: {self.wallet.hotkey.ss58_address}")

        bt.logging.info(f"Connecting to network: {self.config.subtensor.network}")
        self.subtensor = bt.subtensor(config=self.config)
        bt.logging.info(f"Connected to {self.config.subtensor.network}")

        bt.logging.info(f"Loading metagraph for netuid: {self.config.netuid}")
        self.metagraph = self.subtensor.metagraph(self.config.netuid)
        bt.logging.info(f"Metagraph loaded: {len(self.metagraph.hotkeys)} neurons")

        if self.wallet.hotkey.ss58_address not in self.metagraph.hotkeys:
            bt.logging.error(
                f"Validator not registered: {self.wallet}. "
                f"Run btcli subnets register and try again."
            )
            exit()

        self.my_subnet_uid = self.metagraph.hotkeys.index(
            self.wallet.hotkey.ss58_address
        )
        bt.logging.info(f"Validator registered with UID: {self.my_subnet_uid}")

        bt.logging.info("Setting up dendrite...")
        self.dendrite = bt.dendrite(wallet=self.wallet)

        self.scores = torch.zeros(len(self.metagraph.hotkeys), dtype=torch.float32)
        self.score_tracker = ScoreTracker(
            num_miners=len(self.metagraph.hotkeys),
            alpha=GENOMICS_CONFIG["ema_alpha"],
            boost_factor=GENOMICS_CONFIG["boost_factor"]
        )
        bt.logging.info(f"Scores initialized for {len(self.metagraph.hotkeys)} miners")

        self.setup_genomics_components()

        self.task_mode = self.config.task_mode

        bt.logging.info(f"Validator initialization complete")
        bt.logging.info(f"Network: {self.config.subtensor.network}, Netuid: {self.config.netuid}")
        bt.logging.info(f"Task mode: {self.task_mode}, Docker: {is_docker_available()}")

    def setup_genomics_components(self):
        """Setup genomics-specific components."""
        manifest_path = get_manifest_path()
        if manifest_path and os.path.exists(manifest_path):
            bt.logging.info(f"Loading manifest from: {manifest_path}")
            self.task_generator = GenomicsTaskGenerator(manifest_path=manifest_path)
        else:
            bt.logging.warning("No manifest found, using synthetic tasks")
            self.task_generator = GenomicsTaskGenerator()

        self.mutation_injector = MutationInjector(seed=int(time.time()))
        self.happy_scorer = HappyScorer(use_docker=is_docker_available())

        self.current_tasks = {}
        self.task_truth_mapping = {}

        self._cleanup_old_files()

    def _cleanup_old_files(self, max_age_hours: int = 24):
        """Clean up old mutated BAMs and truth VCFs (24h+)."""
        import time

        cutoff_time = time.time() - (max_age_hours * 3600)
        dirs_to_clean = [
            BASE_DIR / "output" / "mutated_bams",
            BASE_DIR / "output" / "merged_truth",
        ]

        total_cleaned = 0
        total_bytes = 0

        for cleanup_dir in dirs_to_clean:
            if not cleanup_dir.exists():
                continue

            for file_path in cleanup_dir.iterdir():
                try:
                    if file_path.is_file() and file_path.stat().st_mtime < cutoff_time:
                        file_size = file_path.stat().st_size
                        file_path.unlink()
                        total_cleaned += 1
                        total_bytes += file_size
                except Exception as e:
                    bt.logging.debug(f"Cleanup failed for {file_path}: {e}")

        if total_cleaned > 0:
            bt.logging.info(f"Cleaned {total_cleaned} files ({total_bytes / (1024**3):.2f} GB)")

    @staticmethod
    def get_config():
        """Get configuration from argparse and environment."""
        parser = argparse.ArgumentParser(
            description="VANET Validator",
            allow_abbrev=False,
        )

        parser.add_argument(
            "--netuid",
            type=int,
            default=int(os.getenv("NETUID", 2)),
            help="Subnet UID to validate on"
        )

        parser.add_argument(
            "--task_mode",
            type=str,
            choices=["compute", "task", "genomics", "mixed"],
            default="genomics",
            help="Type of tasks to send to miners"
        )
        parser.add_argument(
            "--use_docker",
            action="store_true",
            default=True,
            help="Use Docker for hap.py scoring"
        )
        bt.subtensor.add_args(parser)
        bt.logging.add_args(parser)
        bt.wallet.add_args(parser)
        bt.axon.add_args(parser)

        config = bt.config(parser)
        env_network = os.getenv("NETWORK")
        if env_network:
            config.subtensor.network = env_network

        env_netuid = os.getenv("NETUID")
        if env_netuid:
            try:
                config.netuid = int(env_netuid)
            except ValueError:
                print(f"⚠️ Invalid NETUID env var: {env_netuid}", flush=True)

        env_wallet_name = os.getenv("WALLET_NAME")
        if env_wallet_name:
            config.wallet.name = env_wallet_name

        env_wallet_hotkey = os.getenv("WALLET_HOTKEY")
        if env_wallet_hotkey:
            config.wallet.hotkey = env_wallet_hotkey

        return config

    async def query_miners_genomics(self, axons):
        """Query miners with genomics tasks following GenomeNet workflow."""
        task_data = self.task_generator.generate_task(use_manifest=True)
        task_id = task_data["task_id"]
        region = task_data["regions"][0] if task_data.get("regions") else "chr20:10000000-15000000"

        bt.logging.info(f"Generated task {task_id} for region {region}")

        self.current_tasks[task_id] = task_data
        manifest_paths = self.task_generator.get_truth_for_task(task_id) or {}

        bam_local = BASE_DIR / "datasets" / "bams" / "sample.GRCh38.300x_chr20.bam"
        ref_local = BASE_DIR / "datasets" / "reference" / "chr20.fa"
        truth_vcf_local = BASE_DIR / "datasets" / "truth" / "sample_GRCh38_1_22_v4.2.1_benchmark.chr20.vcf.gz"
        truth_bed_local = BASE_DIR / "datasets" / "truth" / "sample_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.chr20.bed"

        bt.logging.debug(f"Dataset files - BAM: {bam_local.exists()}, Ref: {ref_local.exists()}, "
                        f"Truth VCF: {truth_vcf_local.exists()}, Truth BED: {truth_bed_local.exists()}")
        truth_vcf_path = str(truth_vcf_local) if truth_vcf_local.exists() else manifest_paths.get("truth_vcf_s3")
        truth_bed_path = str(truth_bed_local) if truth_bed_local.exists() else manifest_paths.get("truth_bed_s3")
        ref_path = str(ref_local) if ref_local.exists() else manifest_paths.get("reference_fasta_s3")

        self.task_truth_mapping[task_id] = {
            "truth_vcf": truth_vcf_path,
            "truth_bed": truth_bed_path,
            "reference_fasta": ref_path,
            "region": region
        }

        if not truth_vcf_path:
            bt.logging.error(f"No truth VCF for task {task_id}")

        mutated_bam_path = None
        if GENOMICS_CONFIG.get("num_synthetic_mutations", 0) > 0 and is_docker_available():
            if bam_local.exists() and ref_local.exists():
                output_dir = BASE_DIR / "output" / "mutated_bams"
                output_dir.mkdir(parents=True, exist_ok=True)
                round_timestamp = int(time.time())
                mutated_bam = output_dir / f"{task_id}_{round_timestamp}_mutated.bam"

                new_seed = self.mutation_injector.reseed()
                mutations = self.mutation_injector.inject_mutations(
                    task_id, region, GENOMICS_CONFIG["num_synthetic_mutations"],
                    reference_path=ref_local
                )
                bt.logging.info(f"Generated {len(mutations)} mutations (seed={new_seed})")

                bt.logging.info("Injecting mutations with BamSurgeon (may take 5-15 minutes)")

                if GENOMICS_CONFIG.get("use_bamsurgeon", True):
                    success = self.mutation_injector.inject_mutations_into_bam(
                        bam_local, ref_local, mutations, mutated_bam, region=region
                    )

                    if success and mutated_bam.exists():
                        bt.logging.info("BamSurgeon mutation injection complete")
                        mutated_bam_path = mutated_bam

                        if truth_vcf_path:
                            merged_truth_dir = BASE_DIR / "output" / "merged_truth"
                            merged_truth_dir.mkdir(parents=True, exist_ok=True)
                            merged_truth_vcf = merged_truth_dir / f"{task_id}_{round_timestamp}_truth_merged.vcf.gz"

                            merge_success = self.mutation_injector.create_merged_truth_vcf(
                                Path(truth_vcf_path), mutations, merged_truth_vcf, ref_local
                            )

                            if merge_success and merged_truth_vcf.exists():
                                self.task_truth_mapping[task_id]["truth_vcf"] = str(merged_truth_vcf)
                                bt.logging.info("Merged truth VCF created")
                            else:
                                bt.logging.warning("Merge failed, using original truth")
                    else:
                        bt.logging.warning("BamSurgeon failed, using original BAM")
            else:
                bt.logging.warning("BAM/reference not found, skipping mutations")
        else:
            bt.logging.info("Mutation injection disabled")

        if mutated_bam_path and mutated_bam_path.exists():
            data_locator = {"path": str(mutated_bam_path)}
            data_source = "mutated BAM"
        elif bam_local.exists():
            data_locator = {"path": str(bam_local)}
            data_source = "original BAM"
        else:
            ref_paths = self.task_generator.get_reference_paths()
            bam_s3_url = ref_paths.get("bam", {}).get("s3", "")
            data_locator = {"uri": bam_s3_url} if bam_s3_url else task_data.get("data_locator", {})
            data_source = "S3"

        synapse = GenomicsTaskSynapse(
            task_id=task_id,
            task_type="variant_calling",
            data_locator=data_locator,
            regions=task_data.get("regions", []),
            ref_build=task_data.get("ref_build", "GRCh38"),
            input_type=task_data.get("input_type", "BAM"),
            output_format="VCF"
        )

        bt.logging.info(f"Sending task to {len(axons)} miners (source: {data_source}, timeout: {GENOMICS_CONFIG['variant_calling_timeout']//60} min)")

        try:
            query_start = time.time()
            responses = await self.dendrite.forward(
                axons=axons,
                synapse=synapse,
                deserialize=True,
                timeout=GENOMICS_CONFIG["variant_calling_timeout"]
            )
            query_elapsed = time.time() - query_start

            valid_responses = sum(1 for r in responses if r is not None and isinstance(r, dict) and r.get("vcf_output"))
            bt.logging.info(f"Received {valid_responses}/{len(responses)} valid responses in {query_elapsed:.1f}s")

            return synapse, responses
        except Exception as e:
            bt.logging.error(f"Query error: {e}")
            return synapse, []

    def score_genomics_responses(self, synapse, responses):
        """Score miner VCF outputs using hap.py validation."""
        rewards = torch.zeros(len(responses), dtype=torch.float32)
        task_id = synapse.task_id
        region = synapse.regions[0] if synapse.regions else None
        truth_data = self.task_truth_mapping.get(task_id)

        bt.logging.info(f"Scoring {len(responses)} responses for task {task_id}, region {region}")

        if not truth_data or not truth_data.get("truth_vcf"):
            bt.logging.error(f"No truth data for task {task_id}")
            return rewards

        for i, response in enumerate(responses):
            if response is None or not isinstance(response, dict):
                rewards[i] = 0.0
                bt.logging.debug(f"Miner {i}: no response")
                continue

            vcf_output = response.get("vcf_output")
            processing_error = response.get("processing_error")
            metadata = response.get("metadata", {})

            if processing_error or not vcf_output:
                rewards[i] = 0.0
                bt.logging.debug(f"Miner {i}: {processing_error or 'No VCF'}")
                continue

            try:
                variant_count = sum(1 for line in vcf_output.split('\n') if line and not line.startswith('#'))

                output_dir = BASE_DIR / "output" / "scoring" / task_id
                output_dir.mkdir(parents=True, exist_ok=True)

                query_vcf_path = output_dir / f"miner_{i}.vcf"
                with open(query_vcf_path, 'w') as f:
                    f.write(vcf_output)

                scoring_start = time.time()
                metrics = self.happy_scorer._score_with_happy(
                    truth_vcf=truth_data.get("truth_vcf"),
                    query_vcf=str(query_vcf_path),
                    reference_fasta=truth_data.get("reference_fasta"),
                    confident_bed=truth_data.get("truth_bed"),
                    region=region
                )
                scoring_elapsed = time.time() - scoring_start

                score = AdvancedScorer.compute_advanced_score(metrics)
                rewards[i] = score / 100.0

                bt.logging.info(f"Miner {i}: {variant_count} variants, F1_snp={metrics['f1_snp']:.3f}, "
                               f"F1_indel={metrics['f1_indel']:.3f}, score={score:.2f} ({scoring_elapsed:.1f}s)")

            except Exception as e:
                bt.logging.error(f"Scoring error miner {i}: {e}")
                rewards[i] = 0.0

        bt.logging.info(f"Scoring complete: {(rewards > 0).sum().item()}/{len(responses)} valid, "
                       f"avg={rewards.mean().item():.4f}, best={rewards.max().item():.4f}")

        return rewards

    async def query_and_score(self):
        """Main query and scoring loop for one round."""
        # Get active miners
        active_miners = self.get_active_miners()

        if len(active_miners) == 0:
            bt.logging.warning("No active miners found")
            return

        # Select task mode
        if self.task_mode == "genomics":
            synapse, responses = await self.query_miners_genomics(active_miners)
            rewards = self.score_genomics_responses(synapse, responses)
        elif self.task_mode == "mixed":
            # Randomly choose between genomics and compute tasks
            if random.random() < 0.7:  # 70% genomics, 30% compute
                synapse, responses = await self.query_miners_genomics(active_miners)
                rewards = self.score_genomics_responses(synapse, responses)
            else:
                synapse, responses = await self.query_miners_compute(active_miners)
                rewards = self.score_responses_compute(synapse, responses)
        else:
            # Fall back to compute tasks
            synapse, responses = await self.query_miners_compute(active_miners)
            rewards = self.score_responses_compute(synapse, responses)

        # Update scores with EMA and innovation boost
        timestamp = time.time()
        for i, (axon, reward) in enumerate(zip(active_miners, rewards)):
            uid = self.metagraph.axons.index(axon)

            # Update with tracker (includes EMA and boost)
            effective_score = self.score_tracker.update(uid, reward.item(), timestamp)
            self.scores[uid] = effective_score

            bt.logging.debug(f"Miner {uid}: raw={reward.item():.3f}, effective={effective_score:.3f}")

    async def query_miners_compute(self, axons):
        """Query miners with compute tasks (fallback)."""
        operations = ["add", "subtract", "multiply", "divide", "power", "sqrt"]
        operation = random.choice(operations)

        num1 = round(random.uniform(1, 100), 2)
        num2 = round(random.uniform(1, 100), 2) if operation != "sqrt" else None

        if operation == "add":
            problem = f"{num1} + {num2}"
        elif operation == "subtract":
            problem = f"{num1} - {num2}"
        elif operation == "multiply":
            problem = f"{num1} × {num2}"
        elif operation == "divide":
            problem = f"{num1} ÷ {num2}"
        elif operation == "power":
            num2 = round(random.uniform(1, 4), 1)
            problem = f"{num1} ^ {num2}"
        else:  # sqrt
            problem = f"√{num1}"

        synapse = ComputeSynapse(
            operation=operation,
            problem=problem,
            num1=num1,
            num2=num2
        )

        bt.logging.info(f"Sending compute task: {problem}")

        responses = await self.dendrite.forward(
            axons=axons,
            synapse=synapse,
            deserialize=True,
            timeout=12
        )

        return synapse, responses

    def score_responses_compute(self, synapse, responses):
        """Score compute responses."""
        rewards = torch.zeros(len(responses), dtype=torch.float32)

        # Calculate correct answer
        if synapse.operation == "add":
            correct = synapse.num1 + synapse.num2
        elif synapse.operation == "subtract":
            correct = synapse.num1 - synapse.num2
        elif synapse.operation == "multiply":
            correct = synapse.num1 * synapse.num2
        elif synapse.operation == "divide":
            correct = synapse.num1 / synapse.num2 if synapse.num2 != 0 else float('inf')
        elif synapse.operation == "power":
            correct = synapse.num1 ** synapse.num2
        else:  # sqrt
            correct = synapse.num1 ** 0.5

        for i, response in enumerate(responses):
            if response is None or response == 0.0:
                rewards[i] = 0.0
            else:
                error = abs(response - correct)
                relative_error = error / (abs(correct) + 1e-10)

                if error < 0.01:
                    rewards[i] = 1.0
                else:
                    rewards[i] = max(0.0, 1.0 - relative_error * 2)

        return rewards

    def get_active_miners(self):
        """Get list of active miners to query."""
        # Get all axons except our own
        axons = []
        for i, axon in enumerate(self.metagraph.axons):
            if i != self.my_subnet_uid:  # Skip ourselves
                # Check minimum stake requirement
                if self.metagraph.S[i] >= VALIDATOR_CONFIG["min_stake_requirement"]:
                    axons.append(axon)

        # Limit to batch size
        if len(axons) > VALIDATOR_CONFIG["query_batch_size"]:
            axons = random.sample(axons, VALIDATOR_CONFIG["query_batch_size"])

        return axons

    def set_weights(self):
        """Set weights on the blockchain."""
        try:
            weights = self.score_tracker.get_normalized_weights()
            weights_tensor = torch.FloatTensor(weights)
            uids = torch.arange(len(self.metagraph.hotkeys))

            bt.logging.info(f"Setting weights: min={weights_tensor.min():.4f}, "
                           f"max={weights_tensor.max():.4f}, "
                           f"mean={weights_tensor.mean():.4f}")

            self.subtensor.set_weights(
                netuid=self.config.netuid,
                wallet=self.wallet,
                uids=uids,
                weights=weights_tensor,
                wait_for_finalization=False,
                wait_for_inclusion=False,
            )

            bt.logging.success("Weights set successfully")

        except Exception as e:
            bt.logging.error(f"Failed to set weights: {e}")

    async def run(self):
        """Main validator loop."""
        task_interval_mins = GENOMICS_CONFIG["task_interval"] // 60

        bt.logging.info(f"Validator running - Network: {self.config.subtensor.network}, "
                       f"Netuid: {self.config.netuid}, Task mode: {self.task_mode}")
        bt.logging.info(f"Task interval: {task_interval_mins} min, Miner timeout: {GENOMICS_CONFIG['variant_calling_timeout']//60} min, "
                       f"Docker: {is_docker_available()}")

        step = 0

        while True:
            try:
                round_start = time.time()
                bt.logging.info(f"Round {step} started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

                active_miners = self.get_active_miners()
                bt.logging.info(f"Active miners: {len(active_miners)}")

                if len(active_miners) == 0:
                    bt.logging.warning("No active miners found, skipping round")
                else:
                    for i, axon in enumerate(active_miners):
                        uid = self.metagraph.axons.index(axon)
                        bt.logging.debug(f"Miner {uid}: {axon.ip}:{axon.port}")

                    await self.query_and_score()

                if step % VALIDATOR_CONFIG["scoring_interval"] == 0 and step > 0:
                    bt.logging.info("Updating weights on chain...")
                    self.set_weights()

                    stats = self.score_tracker.get_stats()
                    bt.logging.info(f"Performance stats - Global best: {stats['global_best']:.4f}, "
                                   f"Active boosts: {stats['active_boosts']}")

                if step % 100 == 0:
                    self.metagraph.sync(subtensor=self.subtensor)
                    bt.logging.info(f"Metagraph synced - {len(self.metagraph.hotkeys)} neurons")

                round_elapsed = time.time() - round_start
                step += 1

                next_round = datetime.now().timestamp() + GENOMICS_CONFIG["task_interval"]
                next_round_str = datetime.fromtimestamp(next_round).strftime('%H:%M:%S')

                bt.logging.info(f"Round {step-1} complete in {round_elapsed:.1f}s. "
                               f"Next round at {next_round_str} (in {task_interval_mins} min)")

                await asyncio.sleep(GENOMICS_CONFIG["task_interval"])

            except KeyboardInterrupt:
                bt.logging.info(f"Validator shutting down - Completed {step} rounds, "
                               f"processed {len(self.current_tasks)} tasks")
                break
            except Exception as e:
                bt.logging.error(f"Error in main loop: {e}. Retrying in 10 seconds...")
                await asyncio.sleep(10)


def main():
    """Main entry point."""
    validator = Validator()
    asyncio.run(validator.run())


if __name__ == "__main__":
    main()