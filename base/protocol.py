"""
VANET Protocol Definition

This module defines the communication protocol (synapse) between validators and miners.
The synapse is used to send tasks from validators to miners and receive responses.
"""

import bittensor as bt
from typing import Optional, Dict, List, Any


class TaskSynapse(bt.Synapse):
    """
    A simple synapse protocol for sending tasks to miners.

    This synapse carries a task string from the validator to the miner,
    and the miner responds with a result string.

    Attributes:
        task (str): The task description sent by the validator to the miner
        result (Optional[str]): The result computed by the miner (initially None)
    """

    # Required field: task sent by the validator
    task: str

    # Optional field: result populated by the miner
    result: Optional[str] = None

    def deserialize(self) -> str:
        """
        Deserialize the response from the miner.

        Returns:
            str: The result string from the miner, or empty string if None
        """
        return self.result or ""


class ComputeSynapse(bt.Synapse):
    """
    A synapse for mathematical computation tasks.

    The validator sends a math problem, and the miner solves it and returns the answer.

    Attributes:
        operation (str): The operation to perform (e.g., "add", "subtract", "multiply", "divide", "power", "sqrt")
        problem (str): Human-readable math problem (e.g., "45.23 + 12.87")
        num1 (float): First number in the operation
        num2 (Optional[float]): Second number in the operation (None for unary operations like sqrt)
        answer (Optional[float]): The computed result from the miner
    """

    operation: str
    problem: str
    num1: float
    num2: Optional[float] = None
    answer: Optional[float] = None

    def deserialize(self) -> float:
        """
        Deserialize the response from the miner.

        Returns:
            float: The computed answer, or 0.0 if None
        """
        return self.answer if self.answer is not None else 0.0


class GenomicsTaskSynapse(bt.Synapse):
    """
    A synapse for genomic variant calling tasks.

    The validator sends a genomic analysis task to the miner, including:
    - Task ID for tracking
    - BAM file location (S3 or local path)
    - Genomic region to analyze (e.g., chr20:10000000-15000000)
    - Reference genome build (e.g., GRCh38)

    The miner responds with:
    - VCF file content or path containing predicted variants
    - Processing metadata (optional)

    Attributes:
        task_id (str): Unique identifier for this task
        task_type (str): Type of genomics task ("variant_calling")
        data_locator (Dict): Location of input data (BAM file)
        regions (List[str]): Genomic regions to analyze
        ref_build (str): Reference genome build (default: GRCh38)
        input_type (str): Type of input data (BAM or FASTQ)
        output_format (str): Expected output format (VCF)
        vcf_output (Optional[str]): VCF content or file path returned by miner
        metadata (Optional[Dict]): Processing metadata from miner
    """

    # Required fields from validator
    task_id: str
    task_type: str = "variant_calling"
    data_locator: Dict[str, Any]  # {"uri": "s3://bucket/path.bam"} or {"path": "/local/path.bam"}
    regions: List[str]  # ["chr20:10000000-15000000"]
    ref_build: str = "GRCh38"
    input_type: str = "BAM"
    output_format: str = "VCF"

    # Optional fields populated by miner
    vcf_output: Optional[str] = None  # VCF content or path
    metadata: Optional[Dict[str, Any]] = None  # {"tool": "GATK", "version": "4.5.0", "runtime_seconds": 120}
    processing_error: Optional[str] = None  # Error message if processing failed

    def deserialize(self) -> Dict[str, Any]:
        """
        Deserialize the response from the miner.

        Returns:
            Dict: Contains VCF output and metadata, or error information
        """
        return {
            "vcf_output": self.vcf_output,
            "metadata": self.metadata or {},
            "processing_error": self.processing_error,
            "task_id": self.task_id
        }
