"""
File utilities for downloading and managing genomics data files.

Handles S3 downloads, local file resolution, and data caching.
"""

from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger(__name__)

try:
    import boto3
    from botocore.exceptions import ClientError, NoCredentialsError
    HAS_BOTO3 = True
except ImportError:
    HAS_BOTO3 = False
    logger.warning("boto3 not installed - S3 features disabled")


def download_file(url: str, local_path: Path, use_cache: bool = True) -> Optional[Path]:
    """
    Download file from URL (supports S3 public URLs, HTTP/HTTPS).

    Args:
        url: URL to download from
        local_path: Path to save the file
        use_cache: Whether to use cached version if exists

    Returns:
        Path to downloaded file or None if failed
    """
    import urllib.request

    local_path = Path(local_path)

    # Check cache
    if use_cache and local_path.exists():
        logger.info(f"Using cached file: {local_path}")
        return local_path

    # Create directory if needed
    local_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        logger.info(f"Downloading {url} to {local_path}")
        urllib.request.urlretrieve(url, local_path)
        logger.info(f"Downloaded successfully")
        return local_path
    except Exception as e:
        logger.error(f"Failed to download {url}: {e}")
        return None


def ensure_s3_file(file_path: str, file_type: str = "data") -> Path:
    """
    Ensure a file exists locally, downloading from S3 if necessary.

    Args:
        file_path: Local file path or filename
        file_type: Type of file ("reference", "bam", "truth", "data")

    Returns:
        Path to local file (downloaded if necessary)
    """
    from base.genomics_config import GENOMICS_CONFIG
    from base.config import DataConfig

    data_config = DataConfig()
    BASE_DIR = data_config.base_dir

    local_path = Path(file_path)

    # If absolute path and exists, return it
    if local_path.is_absolute() and local_path.exists():
        return local_path

    # Try relative to BASE_DIR
    if not local_path.is_absolute():
        local_path = BASE_DIR / file_path

    # If exists locally, return it
    if local_path.exists():
        logger.info(f"Found local file: {local_path}")
        return local_path

    # Determine S3 URL based on file type and name
    s3_base = GENOMICS_CONFIG["s3_base_url"]
    filename = local_path.name

    if file_type == "reference" or "chr20.fa" in filename:
        s3_url = f"{s3_base}/reference/{filename}"
    elif file_type == "bam" or filename.endswith(".bam") or filename.endswith(".bai"):
        s3_url = f"{s3_base}/bams/{filename}"
    elif file_type == "truth" or filename.endswith(".vcf") or filename.endswith(".bed"):
        s3_url = f"{s3_base}/truth/{filename}"
    else:
        # Default to data type
        s3_url = f"{s3_base}/{file_type}/{filename}"

    # Download from S3
    logger.info(f"File not found locally, downloading from S3: {s3_url}")
    downloaded = download_file(s3_url, local_path, use_cache=True)

    if downloaded:
        # Also download index files if this is a BAM or reference
        if filename.endswith(".bam"):
            # Download BAM index
            index_url = f"{s3_url}.bai"
            index_path = local_path.with_suffix(".bam.bai")
            download_file(index_url, index_path, use_cache=True)
        elif filename.endswith(".fa"):
            # Download reference indices
            for ext in [".fai", ".dict", ".amb", ".ann", ".bwt", ".pac", ".sa"]:
                index_url = f"{s3_url}{ext}"
                index_path = local_path.with_suffix(f".fa{ext}")
                download_file(index_url, index_path, use_cache=True)
            # Also try to download SDF directory
            sdf_url = f"{s3_base}/reference/chr20.sdf"
            sdf_path = local_path.with_suffix(".sdf")
            # Note: SDF is a directory, needs special handling

    return downloaded
