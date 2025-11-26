"""
Configuration for VANET genomics subnet.

Centralizes all configuration with sensible defaults and environment
variable overrides.
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional


def _get_base_dir() -> Path:
    """Get base directory for data files."""
    env_dir = os.environ.get("VANET_DATA_DIR")
    if env_dir:
        return Path(env_dir)
    return Path(__file__).parent


@dataclass
class NetworkConfig:
    """Bittensor network configuration."""

    netuid: int = 2
    network: str = "test"
    min_stake_threshold: float = 0.0

    @classmethod
    def from_env(cls) -> "NetworkConfig":
        return cls(
            netuid=int(os.environ.get("NETUID", 2)),
            network=os.environ.get("SUBTENSOR_NETWORK", "test"),
            min_stake_threshold=float(os.environ.get("MIN_STAKE", 0.0)),
        )


@dataclass
class ValidatorConfig:
    """Validator-specific configuration."""

    interval_minutes: int = 240  # 4 hours between tasks
    variant_calling_timeout: int = 3600  # 1 hour for miner to complete
    weight_update_interval: int = 100  # rounds
    use_mutations: bool = True
    num_mutations: int = 10

    @classmethod
    def from_env(cls) -> "ValidatorConfig":
        return cls(
            interval_minutes=int(os.environ.get("VALIDATOR_INTERVAL", 240)),
            variant_calling_timeout=int(os.environ.get("VARIANT_TIMEOUT", 3600)),
            weight_update_interval=int(os.environ.get("WEIGHT_UPDATE_INTERVAL", 100)),
            use_mutations=os.environ.get("USE_MUTATIONS", "true").lower() == "true",
            num_mutations=int(os.environ.get("NUM_MUTATIONS", 10)),
        )


@dataclass
class MinerConfig:
    """Miner-specific configuration."""

    tool: str = "gatk"
    threads: int = 4
    timeout: int = 3600  # 1 hour for GATK to complete

    @classmethod
    def from_env(cls) -> "MinerConfig":
        return cls(
            tool=os.environ.get("MINER_TOOL", "gatk"),
            threads=int(os.environ.get("MINER_THREADS", 4)),
            timeout=int(os.environ.get("MINER_TIMEOUT", 3600)),
        )


@dataclass
class DataConfig:
    """Data file paths configuration."""

    base_dir: Path = field(default_factory=_get_base_dir)
    s3_base_url: str = "https://vanetdata.s3.us-east-1.amazonaws.com"

    @property
    def datasets_dir(self) -> Path:
        return self.base_dir / "datasets"

    @property
    def reference_dir(self) -> Path:
        return self.datasets_dir / "reference"

    @property
    def bams_dir(self) -> Path:
        return self.datasets_dir / "bams"

    @property
    def truth_dir(self) -> Path:
        return self.datasets_dir / "truth"

    @property
    def output_dir(self) -> Path:
        return self.base_dir / "output"

    @property
    def manifest_path(self) -> Path:
        return self.base_dir / "base" / "s3_manifest.json"

    # Common file paths
    @property
    def reference_fasta(self) -> Path:
        return self.reference_dir / "chr20.fa"

    @property
    def sample_bam(self) -> Path:
        return self.bams_dir / "sample.GRCh38.300x_chr20.bam"

    @property
    def truth_vcf(self) -> Path:
        return self.truth_dir / "sample_GRCh38_1_22_v4.2.1_benchmark.chr20.vcf.gz"

    @property
    def truth_bed(self) -> Path:
        return self.truth_dir / "sample_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.chr20.bed"

    @classmethod
    def from_env(cls) -> "DataConfig":
        base_dir = os.environ.get("VANET_DATA_DIR")
        return cls(
            base_dir=Path(base_dir) if base_dir else _get_base_dir(),
            s3_base_url=os.environ.get(
                "S3_BASE_URL",
                "https://vanetdata.s3.us-east-1.amazonaws.com"
            ),
        )


@dataclass
class Config:
    """Complete subnet configuration."""

    network: NetworkConfig = field(default_factory=NetworkConfig)
    validator: ValidatorConfig = field(default_factory=ValidatorConfig)
    miner: MinerConfig = field(default_factory=MinerConfig)
    data: DataConfig = field(default_factory=DataConfig)

    @classmethod
    def from_env(cls) -> "Config":
        """Load configuration from environment variables."""
        return cls(
            network=NetworkConfig.from_env(),
            validator=ValidatorConfig.from_env(),
            miner=MinerConfig.from_env(),
            data=DataConfig.from_env(),
        )

    def validate(self) -> bool:
        """Validate configuration."""
        # Check required directories exist
        if not self.data.base_dir.exists():
            return False
        return True


# Global config instance
_config: Optional[Config] = None


def get_config() -> Config:
    """Get the global configuration instance."""
    global _config
    if _config is None:
        _config = Config.from_env()
    return _config


def set_config(config: Config):
    """Set the global configuration instance."""
    global _config
    _config = config
