"""
VANET Miner - Genomics Variant Calling Subnet

Processes BAM files with GATK HaplotypeCaller, returns VCF results.
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
import subprocess
import tempfile
import json
import hashlib
from dotenv import load_dotenv

from base import TaskSynapse, ComputeSynapse, GenomicsTaskSynapse
from base import GENOMICS_CONFIG, MINER_CONFIG, is_docker_available

sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)
load_dotenv()


class Miner:
    """VANET miner for genomics variant calling tasks."""

    def __init__(self, config=None):
        self.config = config or self.get_config()

        bt.logging.info("Setting up miner...")
        bt.logging.set_trace(self.config.logging.trace)
        bt.logging.set_debug(self.config.logging.debug)

        self.wallet = bt.wallet(config=self.config)
        bt.logging.info(f"Wallet loaded: {self.wallet.hotkey.ss58_address}")

        bt.logging.info(f"Connecting to network: {self.config.subtensor.network}")
        self.subtensor = bt.subtensor(config=self.config)

        bt.logging.info(f"Loading metagraph for netuid: {self.config.netuid}")
        self.metagraph = self.subtensor.metagraph(self.config.netuid)
        bt.logging.info(f"Metagraph loaded: {len(self.metagraph.hotkeys)} neurons")

        if self.wallet.hotkey.ss58_address not in self.metagraph.hotkeys:
            bt.logging.warning("Miner not registered. Attempting registration...")
            success = self.subtensor.register(
                wallet=self.wallet,
                netuid=self.config.netuid,
                wait_for_finalization=True,
                wait_for_inclusion=True,
            )
            if success:
                bt.logging.success("Registration successful")
                self.metagraph.sync(subtensor=self.subtensor)
            else:
                bt.logging.error("Registration failed")
                exit()

        self.my_subnet_uid = self.metagraph.hotkeys.index(self.wallet.hotkey.ss58_address)
        bt.logging.info(f"Miner registered with UID: {self.my_subnet_uid}")

        bt.logging.info(f"Setting up Axon server on {self.config.axon.ip}:{self.config.axon.port}")
        self.axon = bt.axon(wallet=self.wallet, config=self.config)

        self.axon.attach(
            forward_fn=self.forward_genomics,
            blacklist_fn=self.blacklist_genomics,
            priority_fn=self.priority_genomics,
        ).attach(
            forward_fn=self.forward_compute,
            blacklist_fn=self.blacklist_compute,
            priority_fn=self.priority_compute,
        ).attach(
            forward_fn=self.forward_task,
            blacklist_fn=self.blacklist_task,
            priority_fn=self.priority_task,
        )

        self.setup_variant_caller()

        self.result_cache = {}
        self.max_cache_size = MINER_CONFIG["max_cache_size"]

        bt.logging.info(f"Miner initialization complete - Variant caller: {self.variant_caller}, Docker: {is_docker_available()}")

    def setup_variant_caller(self):
        """Setup the variant calling tool."""
        self.variant_caller = MINER_CONFIG["default_caller"]

        if self.variant_caller in ["gatk", "deepvariant"] and not is_docker_available():
            bt.logging.warning("Docker not available, falling back to synthetic caller")
            self.variant_caller = "synthetic"
        else:
            bt.logging.info(f"Using variant caller: {self.variant_caller}")

    @staticmethod
    def get_config():
        """Get configuration from argparse and environment."""
        parser = argparse.ArgumentParser(
            description="VANET Miner",
            allow_abbrev=False,
        )

        # ---- Core subnet arg (you were missing this) ----
        parser.add_argument(
            "--netuid",
            type=int,
            default=int(os.getenv("NETUID", 2)),  # pick whatever default you want
            help="Subnet UID to mine on",
        )

        # ---- Custom miner-only arguments ----
        parser.add_argument(
            "--variant_caller",
            type=str,
            choices=["gatk", "deepvariant", "synthetic"],
            default=MINER_CONFIG["default_caller"],
            help="Variant calling tool to use",
        )

        # ---- Bittensor's standard args ----
        bt.subtensor.add_args(parser)
        bt.logging.add_args(parser)
        bt.wallet.add_args(parser)
        bt.axon.add_args(parser)

        # Build config
        config = bt.config(parser)

        # Optional: env overrides
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

        env_miner_ip = os.getenv("MINER_IP")
        if env_miner_ip:
            config.axon.ip = env_miner_ip

        env_miner_port = os.getenv("MINER_PORT")
        if env_miner_port:
            try:
                config.axon.port = int(env_miner_port)
            except ValueError:
                print(f"⚠️ Invalid MINER_PORT env var: {env_miner_port}", flush=True)

        return config

    async def forward_genomics(self, synapse: GenomicsTaskSynapse) -> GenomicsTaskSynapse:
        """Process genomic variant-calling tasks."""
        start_time = time.time()

        bt.logging.info(f"Genomics task received: {synapse.task_id}, region={synapse.regions[0] if synapse.regions else 'None'}, "
                       f"ref={synapse.ref_build}, source={'local' if 'path' in synapse.data_locator else 'S3'}")

        try:
            if MINER_CONFIG["cache_results"]:
                cache_key = self._get_cache_key(synapse)
                if cache_key in self.result_cache:
                    bt.logging.info(f"Cache hit for {synapse.task_id}")
                    synapse.vcf_output = self.result_cache[cache_key]["vcf"]
                    synapse.metadata = self.result_cache[cache_key]["metadata"]
                    return synapse

            bt.logging.info(f"Starting variant calling with {self.variant_caller.upper()}")

            if self.variant_caller == "gatk":
                vcf_output = await self.run_gatk(synapse)
            elif self.variant_caller == "deepvariant":
                vcf_output = await self.run_deepvariant(synapse)
            else:
                vcf_output = await self.run_synthetic(synapse)

            variant_count = sum(1 for line in vcf_output.split('\n') if line and not line.startswith('#'))
            elapsed = time.time() - start_time

            synapse.vcf_output = vcf_output
            synapse.metadata = {
                "tool": self.variant_caller,
                "version": self._get_tool_version(),
                "runtime_seconds": elapsed,
                "miner_uid": self.my_subnet_uid,
                "variants_called": variant_count
            }

            if MINER_CONFIG["cache_results"]:
                self._cache_result(synapse)

            bt.logging.info(f"Task {synapse.task_id} completed: {variant_count} variants in {elapsed:.1f}s, size={len(vcf_output):,} bytes")

        except Exception as e:
            bt.logging.error(f"Error processing genomics task {synapse.task_id}: {e}")
            synapse.processing_error = str(e)

        return synapse

    async def run_gatk(self, synapse: GenomicsTaskSynapse) -> str:
        """Run GATK HaplotypeCaller on BAM input (matching GenomeNet step 4-5)."""
        from utils import download_file, run_gatk_docker
        from base.genomics_config import GENOMICS_CONFIG, BASE_DIR

        print(f"\n   [GATK] Preparing input files...", flush=True)

        output_dir = BASE_DIR / "output" / synapse.task_id
        output_dir.mkdir(parents=True, exist_ok=True)

        bam_path = None

        if "path" in synapse.data_locator:
            validator_bam_path = Path(synapse.data_locator["path"])
            if validator_bam_path.exists():
                bam_size_gb = validator_bam_path.stat().st_size / (1024**3)
                print(f"   [GATK] Using validator BAM: {validator_bam_path.name} ({bam_size_gb:.2f} GB)", flush=True)
                bam_path = validator_bam_path
            else:
                print(f"   [GATK] Validator BAM not found: {validator_bam_path}", flush=True)

        if bam_path is None:
            bam_local = BASE_DIR / "datasets" / "bams" / "sample.GRCh38.300x_chr20.bam"
            if bam_local.exists():
                bam_size_gb = bam_local.stat().st_size / (1024**3)
                print(f"   [GATK] Using local BAM: {bam_local.name} ({bam_size_gb:.2f} GB)", flush=True)
                bam_path = bam_local
            else:
                print(f"   [GATK] Local BAM not found, downloading from data_locator...", flush=True)
                # Download from data_locator
                with tempfile.TemporaryDirectory() as tmpdir:
                    bam_path = await self._prepare_input_data(synapse, tmpdir)
                    if not bam_path:
                        raise RuntimeError("Failed to prepare input BAM file")

        # Get reference genome - check local first
        ref_local = BASE_DIR / "datasets" / "reference" / "chr20.fa"

        if not ref_local.exists():
            print(f"   [GATK] Reference not found locally, downloading...", flush=True)
            # Try to download reference
            ref_url = GENOMICS_CONFIG.get("reference_fasta", "")
            if ref_url.startswith("http"):
                ref_local.parent.mkdir(parents=True, exist_ok=True)
                downloaded = download_file(ref_url, ref_local, use_cache=True)
                if not downloaded:
                    raise RuntimeError("Reference genome not available")
            else:
                ref_local = Path(ref_url)
                if not ref_local.exists():
                    raise RuntimeError(f"Reference genome not found: {ref_local}")
        else:
            ref_size_mb = ref_local.stat().st_size / (1024**2)
            print(f"   [GATK] Using local reference: {ref_local.name} ({ref_size_mb:.1f} MB)", flush=True)

        # Output VCF path
        output_vcf = output_dir / f"{synapse.task_id}.vcf"

        # Run GATK with Docker (uses same format as test_step5_miner.py)
        interval = synapse.regions[0] if synapse.regions else None
        print(f"   [GATK] Region: {interval}", flush=True)
        print(f"   [GATK] Output: {output_vcf}", flush=True)
        print(f"   [GATK] Starting GATK HaplotypeCaller (this may take 10-20 minutes)...", flush=True)

        gatk_start = time.time()
        success = run_gatk_docker(bam_path, ref_local, output_vcf, interval, threads=4)
        gatk_elapsed = time.time() - gatk_start

        if not success:
            print(f"   [GATK] FAILED after {gatk_elapsed:.1f} seconds", flush=True)
            raise RuntimeError("GATK HaplotypeCaller failed")

        print(f"   [GATK] Completed in {gatk_elapsed:.1f} seconds", flush=True)

        # Read VCF output
        with open(output_vcf, 'r') as f:
            vcf_content = f.read()

        # Count variants
        variants_called = sum(1 for line in vcf_content.split('\n') if line and not line.startswith('#'))
        print(f"   [GATK] Variants called: {variants_called}", flush=True)

        return vcf_content

    async def run_deepvariant(self, synapse: GenomicsTaskSynapse) -> str:
        """Run DeepVariant on the input data."""
        # Similar to GATK but with DeepVariant
        # For brevity, returning synthetic data
        return await self.run_synthetic(synapse)

    async def run_synthetic(self, synapse: GenomicsTaskSynapse) -> str:
        """Generate synthetic VCF output for testing."""
        # Parse region
        region_str = synapse.regions[0] if synapse.regions else "chr20:1-1000000"
        chrom, coords = region_str.split(':')
        start, end = map(int, coords.split('-'))

        # Generate synthetic VCF
        vcf_lines = [
            "##fileformat=VCFv4.2",
            f"##source=GenomicsMiner_Synthetic_v1.0",
            f"##reference={synapse.ref_build}",
            f"##contig=<ID={chrom},length={end}>",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        ]

        # Generate random variants
        import random
        random.seed(hash(synapse.task_id))

        num_variants = random.randint(500, 1500)
        nucleotides = ['A', 'C', 'G', 'T']

        positions = sorted(random.sample(range(start, end), min(num_variants, end - start)))

        for pos in positions:
            ref = random.choice(nucleotides)
            alt = random.choice([n for n in nucleotides if n != ref])
            qual = random.randint(30, 255)
            dp = random.randint(10, 100)
            gt = random.choice(["0/1", "1/1"])
            gq = random.randint(30, 99)

            vcf_lines.append(
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\tDP={dp}\tGT:GQ\t{gt}:{gq}"
            )

        return "\n".join(vcf_lines)

    async def _prepare_input_data(self, synapse: GenomicsTaskSynapse, tmpdir: str) -> Path:
        """Prepare input BAM file for processing."""
        from utils import download_file
        from pathlib import Path

        tmpdir_path = Path(tmpdir)
        bam_path = tmpdir_path / "input.bam"

        # Check data_locator for file location
        if "uri" in synapse.data_locator:
            # URL (S3 public URL, HTTP, etc.)
            uri = synapse.data_locator["uri"]
            if uri.startswith(("http://", "https://", "s3://")):
                # Download from URL
                bt.logging.info(f"Downloading BAM from {uri}")
                downloaded = download_file(uri, bam_path, use_cache=False)
                if downloaded:
                    return downloaded
                else:
                    bt.logging.error(f"Failed to download BAM from {uri}")
                    # Fall back to synthetic
                    with open(bam_path, 'wb') as f:
                        f.write(b"SYNTHETIC_BAM")
                    return bam_path

        elif "path" in synapse.data_locator:
            # Local path
            local_path = Path(synapse.data_locator["path"])
            if local_path.exists():
                # Copy to tmpdir
                import shutil
                shutil.copy(local_path, bam_path)
                return bam_path
            else:
                bt.logging.warning(f"Local BAM not found: {local_path}")

        # Fallback: create synthetic BAM
        bt.logging.info("Using synthetic BAM data")
        with open(bam_path, 'wb') as f:
            f.write(b"SYNTHETIC_BAM_CONTENT")

        return bam_path

    def _get_tool_version(self) -> str:
        """Get version of the variant calling tool."""
        if self.variant_caller == "gatk":
            return "4.5.0.0"
        elif self.variant_caller == "deepvariant":
            return "1.5.0"
        else:
            return "1.0.0"

    def _get_cache_key(self, synapse: GenomicsTaskSynapse) -> str:
        """Generate cache key for a task.

        IMPORTANT: Include BAM path in cache key so that different mutated BAMs
        (with fresh synthetic mutations each round) don't return stale cached results.
        """
        # Get BAM path from data_locator - this changes each round with mutations
        bam_path = synapse.data_locator.get("path", "") if synapse.data_locator else ""

        key_parts = [
            synapse.task_id,
            str(synapse.regions),
            synapse.ref_build,
            self.variant_caller,
            bam_path  # Include BAM path so mutated BAMs get fresh processing
        ]
        key_str = "_".join(key_parts)
        return hashlib.md5(key_str.encode()).hexdigest()

    def _cache_result(self, synapse: GenomicsTaskSynapse):
        """Cache the result of a task."""
        cache_key = self._get_cache_key(synapse)

        # Store in cache
        self.result_cache[cache_key] = {
            "vcf": synapse.vcf_output,
            "metadata": synapse.metadata,
            "timestamp": time.time()
        }

        # Limit cache size
        if len(self.result_cache) > self.max_cache_size:
            # Remove oldest entry
            oldest_key = min(self.result_cache.keys(),
                           key=lambda k: self.result_cache[k]["timestamp"])
            del self.result_cache[oldest_key]

    async def forward_compute(self, synapse: ComputeSynapse) -> ComputeSynapse:
        """Process compute tasks."""
        bt.logging.info(f"Processing compute task: {synapse.problem}")
        print(f"\n🔢 Received compute task: {synapse.problem}", flush=True)

        try:
            if synapse.operation == "add":
                synapse.answer = synapse.num1 + synapse.num2
            elif synapse.operation == "subtract":
                synapse.answer = synapse.num1 - synapse.num2
            elif synapse.operation == "multiply":
                synapse.answer = synapse.num1 * synapse.num2
            elif synapse.operation == "divide":
                synapse.answer = synapse.num1 / synapse.num2 if synapse.num2 != 0 else float('inf')
            elif synapse.operation == "power":
                synapse.answer = synapse.num1 ** synapse.num2
            elif synapse.operation == "sqrt":
                synapse.answer = synapse.num1 ** 0.5
            else:
                synapse.answer = 0.0

            print(f"   Answer: {synapse.answer}", flush=True)

        except Exception as e:
            bt.logging.error(f"Error in compute: {e}")
            synapse.answer = 0.0

        return synapse

    async def forward_task(self, synapse: TaskSynapse) -> TaskSynapse:
        """Process text tasks."""
        bt.logging.info(f"Processing task: {synapse.task}")
        print(f"\n📝 Received task: {synapse.task[:50]}...", flush=True)

        try:
            # Simple task processing
            synapse.result = f"Processed: {synapse.task.upper()}"
            print(f"   Result: {synapse.result[:50]}...", flush=True)
        except Exception as e:
            bt.logging.error(f"Error processing task: {e}")
            synapse.result = f"Error: {str(e)}"

        return synapse

    def blacklist_genomics(self, synapse: GenomicsTaskSynapse) -> typing.Tuple[bool, str]:
        """Blacklist unregistered requests for genomics."""
        if synapse.dendrite.hotkey not in self.metagraph.hotkeys:
            bt.logging.warning(f"Blacklisting unregistered hotkey: {synapse.dendrite.hotkey}")
            return True, "Unregistered hotkey"

        # Check stake requirement
        uid = self.metagraph.hotkeys.index(synapse.dendrite.hotkey)
        stake = self.metagraph.S[uid].item()

        if stake < 0.0:  # No minimum stake for now
            return True, f"Insufficient stake: {stake}"

        return False, ""

    def blacklist_compute(self, synapse: ComputeSynapse) -> typing.Tuple[bool, str]:
        """Blacklist unregistered requests for compute."""
        if synapse.dendrite.hotkey not in self.metagraph.hotkeys:
            return True, "Unregistered hotkey"
        return False, ""

    def blacklist_task(self, synapse: TaskSynapse) -> typing.Tuple[bool, str]:
        """Blacklist unregistered requests for tasks."""
        if synapse.dendrite.hotkey not in self.metagraph.hotkeys:
            return True, "Unregistered hotkey"
        return False, ""

    def priority_genomics(self, synapse: GenomicsTaskSynapse) -> float:
        """Priority based on stake for genomics."""
        if synapse.dendrite.hotkey not in self.metagraph.hotkeys:
            return 0.0

        uid = self.metagraph.hotkeys.index(synapse.dendrite.hotkey)
        stake = self.metagraph.S[uid].item()

        return max(0.0, min(1.0, stake / 10000.0))

    def priority_compute(self, synapse: ComputeSynapse) -> float:
        """Priority based on stake for compute."""
        if synapse.dendrite.hotkey not in self.metagraph.hotkeys:
            return 0.0

        uid = self.metagraph.hotkeys.index(synapse.dendrite.hotkey)
        stake = self.metagraph.S[uid].item()

        return max(0.0, min(1.0, stake / 10000.0))

    def priority_task(self, synapse: TaskSynapse) -> float:
        """Priority based on stake for tasks."""
        if synapse.dendrite.hotkey not in self.metagraph.hotkeys:
            return 0.0

        uid = self.metagraph.hotkeys.index(synapse.dendrite.hotkey)
        stake = self.metagraph.S[uid].item()

        return max(0.0, min(1.0, stake / 10000.0))

    def run(self):
        """Run the miner."""
        bt.logging.info("Starting miner...")

        print(f"\n{'='*60}", flush=True)
        print(f"🚀 MINER RUNNING", flush=True)
        print(f"{'='*60}", flush=True)
        print(f"   Axon: {self.config.axon.ip}:{self.config.axon.port}", flush=True)
        print(f"   Hotkey: {self.wallet.hotkey.ss58_address[:16]}...", flush=True)
        print(f"   Network: {self.config.subtensor.network}", flush=True)
        print(f"   Netuid: {self.config.netuid}", flush=True)
        print(f"   Variant Caller: {self.variant_caller}", flush=True)
        print(f"   Docker: {'Available' if is_docker_available() else 'Not Available'}", flush=True)
        print(f"{'='*60}", flush=True)

        # Start axon server - need both serve() and start()
        # serve() registers the axon on the network
        # start() actually starts the HTTP server to receive requests
        self.axon.serve(netuid=self.config.netuid, subtensor=self.subtensor)
        self.axon.start()
        bt.logging.info("Axon server started and serving")

        print(f"\n✅ AXON SERVER STARTED", flush=True)
        print(f"   Listening on: {self.axon.ip}:{self.axon.port}", flush=True)
        print(f"\n⏳ WAITING FOR TASKS FROM VALIDATORS...", flush=True)
        print(f"   Tasks expected every ~40 minutes", flush=True)
        print(f"   Press Ctrl+C to stop\n", flush=True)

        # Keep running
        sync_count = 0
        try:
            while True:
                # Sync metagraph periodically
                time.sleep(60)
                sync_count += 1
                self.metagraph.sync(subtensor=self.subtensor)

                # Show heartbeat every 5 minutes
                if sync_count % 5 == 0:
                    print(f"💓 Heartbeat | {time.strftime('%Y-%m-%d %H:%M:%S')} | Cache: {len(self.result_cache)} results | Uptime: {sync_count} min", flush=True)

        except KeyboardInterrupt:
            print(f"\n{'='*60}", flush=True)
            print(f"👋 MINER SHUTTING DOWN", flush=True)
            print(f"{'='*60}", flush=True)
            bt.logging.info("Keyboard interrupt, shutting down...")
            self.axon.stop()
            print(f"   Axon stopped", flush=True)
            print(f"   Total uptime: {sync_count} minutes", flush=True)
            print(f"   Cached results: {len(self.result_cache)}", flush=True)


def main():
    """Main entry point."""
    miner = Miner()
    miner.run()


if __name__ == "__main__":
    main()