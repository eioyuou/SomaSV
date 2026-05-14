<p align="center">
  <img src="logo/logo.svg" alt="SomaSV Logo" style="max-width: 100%; height: auto;" width="360">
</p>
<p align="center">
  <strong> Systematic integration of long- and short-read sequencing improves somatic structural variant detection </strong>
</p>

<p align="center">
  <a href="https://github.com/eioyuou/SomaSV/releases"><img src="https://img.shields.io/badge/version-0.0.1-blue.svg" alt="Version"></a>
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-≥3.8-green.svg" alt="Python"></a>
  <a href="https://github.com/eioyuou/SomaSV/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-MIT-purple.svg" alt="License"></a>
  <a href="https://github.com/eioyuou/SomaSV"><img src="https://img.shields.io/badge/platform-ONT%20%2F%20HiFi-orange.svg" alt="Platform"></a>
  <a href="https://github.com/eioyuou/SomaSV"><img src="https://img.shields.io/badge/mode-long--read%20%2F%20hybrid-red.svg" alt="Mode"></a>
</p>
---

## Overview

**SomaSV** is a hybrid method for high-accuracy somatic structural variant (SSV) detection from long-read sequencing data. By integrating tumor long-read data with matched normal data composed of both long-read and short-read sequencing, SomaSV achieves superior detection performance while significantly reducing sequencing costs.

SomaSV supports two operating modes:

- **Hybrid mode** — combines 30× tumor long-read, 10× normal long-read, and 30× normal short-read data for cost-effective, high-accuracy somatic SV calling.
- **Long-read-only mode** — uses tumor and matched normal long-read data without short-read input.


## Key Features

- **Hybrid sequencing design** — leverages inexpensive short-read data as matched normal to replace a substantial portion of long-read sequencing, making somatic SV detection more accessible and cost-effective.
- **Multi-platform support** — compatible with Oxford Nanopore Technologies (ONT) and PacBio HiFi long-read platforms.
- **High accuracy** — achieves state-of-the-art F1 scores across multiple benchmark datasets and SV types (DEL, DUP, INV, INS, BND).
- **Clinical relevance** — capable of identifying somatic SVs in cancer-associated genes with potential diagnostic and screening value.

## Installation


### Install SomaSV

```bash
git clone https://github.com/eioyuou/SomaSV.git
cd SomaSV
pip install -e .
```

## Recommended Sequencing Design

### Hybrid mode (recommended)

| Data Type | Coverage | Cost Ratio |
| :--- | :---: | :---: |
| Tumor long-read | 30× | High |
| Normal long-read | 10× | Medium |
| Normal short-read | 30× | Low |

This hybrid design reduces overall sequencing cost by approximately **19%** compared to a full long-read-only approach while maintaining or improving detection accuracy.

### Long-read-only mode

| Data Type | Coverage | Cost Ratio |
| :--- | :---: | :---: |
| Tumor long-read | — | High |
| Normal long-read | — | High |

This mode uses standard matched tumor–normal long-read sequencing without short-read data. Coverage can be adjusted based on project requirements and budget.

## Quick Start
The commands below show the standard usage of SomaSV with user-provided input files.
For a minimal runnable example using the bundled demo dataset, see [Demo](#demo).
### Required Data

| File | Description |
| :--- | :--- |
| Tumor long-read BAM | Long-read sequencing alignment of the tumor sample |
| Normal long-read BAM | Long-read sequencing alignment of the matched normal sample |
| Reference genome | Reference FASTA file (e.g., GRCh38) |
| PoN VCF | Panel of Normals VCF for germline SV filtering (e.g., [gnomAD SV v4.1](https://gnomad.broadinstitute.org/downloads#v4-structural-variants)) |
| Normal short-read BAM | *(Hybrid mode only)* Short-read sequencing alignment of the matched normal sample |

### Hybrid mode

```bash
python main.py \
    --tumor-bam tumor_long_read.bam \
    --normal-bam normal_long_read.bam \
    --reference ref.fasta \
    --output-dir results/ \
    --sample-platform HIFI \
    --mode hybrid \
    --short-read-normal-bam normal_short_read.bam \
    --short-read-coverage 30 \
    --pon-vcf gnomad.v4.1.sv.sites.vcf.gz
```

### Long-read-only mode

```bash
python main.py \
    --tumor-bam tumor_long_read.bam \
    --normal-bam normal_long_read.bam \
    --reference ref.fasta \
    --output-dir results/ \
    --sample-platform ONT \
    --mode long-read-only 
```

### Parameters

| Parameter | Required | Description |
| :--- | :---: | :--- |
| `--tumor-bam` | ✓ | Path to tumor long-read BAM file |
| `--normal-bam` | ✓ | Path to normal long-read BAM file |
| `--reference` | ✓ | Path to reference genome FASTA file |
| `--output-dir` | ✓ | Output directory for results |
| `--sample-platform` | ✓ | Sequencing platform: `HIFI` or `ONT` |
| `--mode` | ✓ | Running mode: `hybrid` or `long-read-only` |
| `--pon-vcf` | Hybrid | Path to Panel of Normals (PoN) VCF for germline filtering |
| `--short-read-normal-bam` | Hybrid | Path to normal short-read BAM file |
| `--short-read-coverage` | Hybrid | Coverage of normal short-read data (e.g., `30`) |

### Output

The main output is a standard VCF file located at:

```
results/final_somatic_variants.vcf
```

## Demo

A lightweight demo dataset is provided in `data/demo/` for quick testing and validation of SomaSV.

These demo files can be used directly with the example commands below to verify installation and basic pipeline execution.

### Demo data

The demo files are small subset BAMs derived from chromosome 22 of sample HG008 and are intended for demonstration purposes only. They do not represent complete sequencing datasets.

Example demo files:

- `data/demo/HG008_30X_HiFi_chr22_tumor_demo_subset.bam`
- `data/demo/HG008_10X_HiFi_chr22_normal_demo_subset.bam`
- `data/demo/HG008_30X_illumina_chr22_normal_demo_subset.bam`

### Run the demo

#### Hybrid mode

    python main.py \
        --tumor-bam data/demo/HG008_30X_HiFi_chr22_tumor_demo_subset.bam \
        --normal-bam data/demo/HG008_10X_HiFi_chr22_normal_demo_subset.bam \
        --reference ref.fasta \
        --output-dir demo_results/ \
        --sample-platform HIFI \
        --mode hybrid \
        --short-read-normal-bam data/demo/HG008_30X_illumina_chr22_normal_demo_subset.bam \
        --short-read-coverage 30 \
        --pon-vcf gnomad.v4.1.sv.sites.vcf.gz

#### Long-read-only mode

    python main.py \
        --tumor-bam data/demo/HG008_30X_HiFi_chr22_tumor_demo_subset.bam \
        --normal-bam data/demo/HG008_10X_HiFi_chr22_normal_demo_subset.bam \
        --reference ref.fasta \
        --output-dir demo_results_long_read_only/ \
        --sample-platform HIFI \
        --mode long-read-only

### Demo output

The main output VCF will be written to:

    demo_results/final_somatic_variants.vcf

### Notes

- The demo BAM files contain only a small subset region from chromosome 22.
- These files are intended for demonstration purposes only and do not represent complete datasets.
- For real analyses, please use complete datasets and an appropriate reference genome and PoN resource.


## Citation

If you use SomaSV in your research, please cite:

> **Systematic integration of long- and short-read sequencing improves somatic structural variant detection**
>
> *Manuscript in preparation*

## License

This project is licensed under the MIT License.

## Contact

For questions, bug reports, or feature requests, please open an [issue](https://github.com/eioyuou/SomaSV/issues) on GitHub.
