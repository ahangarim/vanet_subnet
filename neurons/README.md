# VANET Neurons

This folder contains the two main participants in the VANET subnet: **Miners** and **Validators**.

---

## Overview

VANET is a Bittensor subnet that rewards miners for accurately calling genetic variants (finding mutations in DNA). Think of it like a competition where:

- **Validators** are the judges - they create challenges and score the results
- **Miners** are the competitors - they solve the challenges to earn rewards

---

## Miner ([miner.py](miner.py))

### What Does a Miner Do?

A miner receives genomic data (DNA sequencing files) and finds genetic variants (mutations) in them. It's like being given a puzzle and finding all the differences from the original.

### How It Works (Simple Steps)

1. **Start Up**
   - Connect to the Bittensor network
   - Register your hotkey to get a unique ID
   - Start a server to receive tasks from validators

2. **Receive a Task**
   - Get a BAM file (contains DNA sequencing data)
   - Get a genomic region to analyze (e.g., "chromosome 20, position 10M to 15M")

3. **Process the Data**
   - Download the BAM file if needed
   - Run GATK HaplotypeCaller (a specialized tool for finding mutations)
   - This analyzes the DNA data and identifies variants

4. **Return Results**
   - Send back a VCF file (contains the list of mutations found)
   - Include metadata (tool used, runtime, number of variants found)

5. **Get Rewarded**
   - Validators score your accuracy
   - More accurate results = higher rewards
   - Earn TAO tokens based on performance

### Key Features

- **Caching**: Remembers previous results to avoid duplicate work
- **Docker Integration**: Runs GATK in containers for consistency
- **Flexible**: Can use different variant calling tools (GATK, DeepVariant, etc.)

### Running the Miner

```bash
python neurons/miner.py \
  --netuid 1 \
  --subtensor.network test \
  --wallet.name my_wallet \
  --wallet.hotkey my_hotkey
```

---

## Validator ([validator.py](validator.py))

### What Does a Validator Do?

A validator creates challenges for miners, scores their answers, and decides how rewards are distributed. It's like a teacher who:
- Makes up test questions
- Grades the answers
- Assigns grades fairly

### How It Works (Simple Steps)

1. **Start Up**
   - Connect to the Bittensor network
   - Register as a validator
   - Load reference genomic data (the "answer key")

2. **Create a Challenge** (Every 4 Hours)
   - Pick a random region of chromosome 20
   - Inject synthetic mutations
   - Create a "mutated" BAM file with these hidden changes
   - Keep track of what mutations were added (this is the answer key)

3. **Send Tasks to Miners**
   - Send the mutated BAM file to all active miners
   - Give them 1 hour to process it
   - Collect their VCF results (their answers)

4. **Score the Results**
   - Use hap.py (a specialized validation tool) to compare:
     - What the miner found (their VCF)
     - What should have been found (truth VCF = original + synthetic mutations)
   - Calculate accuracy scores:
     - **Precision**: How many found variants are correct?
     - **Recall**: Did they find all the variants?
     - **F1 Score**: Overall accuracy (combines precision and recall)

5. **Update Weights**
   - Track each miner's performance over time (using EMA - exponential moving average)
   - Reward miners who improve the global best with a temporary boost
   - Submit updated weights to the blockchain every 100 rounds
   - Higher weights = more TAO rewards for that miner

### Why Inject Mutations?

This is the **anti-cheating mechanism**:

- **Problem**: Miners could memorize answers or copy each other
- **Solution**: Each round has fresh, random mutations that have never been seen before
- **Result**: Miners must actually run real variant calling - they can't cheat!

### Scoring System

The validator uses a sophisticated scoring system:

1. **SNP Score (70%)**: How well did they find single-letter mutations?
2. **INDEL Score (30%)**: How well did they find insertions/deletions?
3. **Completeness**: Did they analyze the whole region?
4. **Quality**: Are the genotypes biologically realistic?

### Innovation Boost

Miners who achieve a breakthrough (exceed the global best by 5%) get:
- A temporary 15% boost to their score
- The boost decays over 1 hour
- This rewards miners who improve the state of the art!

### Running the Validator

```bash
python neurons/validator.py \
  --netuid 1 \
  --subtensor.network test \
  --wallet.name my_wallet \
  --wallet.hotkey my_hotkey
```

---

## The Complete Cycle

```
┌─────────────┐
│  VALIDATOR  │
└──────┬──────┘
       │
       │ 1. Generate random genomic region
       │ 2. Inject 10-200 synthetic mutations
       │ 3. Create mutated BAM file
       │
       ├───────────────────────────────────┐
       │                                   │
       ▼                                   ▼
  ┌─────────┐                        ┌─────────┐
  │ MINER 1 │                        │ MINER 2 │
  └────┬────┘                        └────┬────┘
       │                                   │
       │ 4. Receive mutated BAM            │
       │ 5. Run GATK HaplotypeCaller       │
       │ 6. Return VCF results             │
       │                                   │
       └───────────────┬───────────────────┘
                       │
                       ▼
              ┌─────────────┐
              │  VALIDATOR  │
              └──────┬──────┘
                     │
                     │ 7. Score with hap.py
                     │ 8. Update EMA scores
                     │ 9. Apply innovation boost
                     │ 10. Set weights on-chain
                     │
                     ▼
              ┌──────────────┐
              │  BLOCKCHAIN  │
              └──────────────┘
                     │
                     │ Miners receive TAO
                     │ based on their weights
                     ▼
```

---

## Configuration

Both miner and validator use configuration from the `base/` folder:

- **GENOMICS_CONFIG**: Window size, mutation counts, timeouts, Docker images
- **MINER_CONFIG**: Which variant calling tool to use, number of threads
- **VALIDATOR_CONFIG**: How often to run rounds, when to update weights

---

## Key Differences

| Aspect | Miner | Validator |
|--------|-------|-----------|
| **Role** | Solve challenges | Create challenges & judge |
| **Input** | Mutated BAM files | Original reference data |
| **Processing** | Run GATK/DeepVariant | Run hap.py scoring |
| **Output** | VCF files | Weights on blockchain |
| **Reward** | Earn TAO based on accuracy | Earn TAO for validating |
| **Runs** | On demand (when tasks arrive) | Every 4 hours |

---

## Requirements

### Both Need:
- Bittensor wallet with registered hotkey
- Docker installed (for running genomics tools)
- Network connection to Bittensor testnet or mainnet

### Miner Needs:
- GATK Docker image: `broadinstitute/gatk:4.5.0.0`
- ~4GB RAM for processing
- Storage for temporary BAM files

### Validator Needs:
- Reference genomic data (~15GB)
- BamSurgeon Docker image: `quay.io/biocontainers/bamsurgeon:1.4.1`
- hap.py Docker image: `mgibio/hap.py:v0.3.12`
- samtools, bcftools Docker images
- ~50GB storage for datasets and temporary files

---

## Tips

### For Miners:
- Use faster variant callers for quicker results
- Make sure you have enough disk space for BAM files
- Monitor your accuracy scores to improve performance

### For Validators:
- Download reference data before starting (datasets/reference, datasets/truth)
- Ensure all Docker images are pulled beforehand
- Monitor disk usage - old mutation files can accumulate

---

## Troubleshooting

**Miner Issues:**
- "BAM index not found" → Run `samtools index` on your BAM file
- "GATK timeout" → Increase timeout in config or use faster hardware
- "No response from validator" → Check network connection and firewall

**Validator Issues:**
- "BamSurgeon failed" → Ensure pysam is installed: `pip install pysam`
- "hap.py returned zero scores" → Check that truth VCF matches reference build
- "Mutations file kept for inspection" → BamSurgeon couldn't create mutated BAM

---

## Learn More

- See `utils/` folder for the genomics processing tools
- See `base/` folder for configuration options
- Check logs for detailed execution information
