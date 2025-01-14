# ugenome-maf-pdf

*All code is written by Adam Kurth (ASU) as a part of a UGenome AI internship.*

This project is a simple demonstration of how to process and analyze genomic data for clinical report generation using MAF files. The project is divided into two scripts: `process-maf-data.py` and `generate-clinical-report.py`. The first script processes MAF files and generates a JSON file of the genetic analysis. The second script processes the clinical report and generates a summary report and visualizations.

## Directory Structure

Directory structure for processing and analyzing genomic data for a clinical report generation.


```bash 
project_root/
├── data/
│   ├── maf/           # Put your .maf files here
│   └── results/       # Output will be created here
└── scripts/
    ├── process_maf_data.py
    └── generate_clinical_report.py
```

## `scripts/`

1. `process-maf-data.py` - Process MAF files and generate .json files of the genetic analysis.
   
2. `process-clinical-report.py` - Process clinical report and generate a summary report and visualizations.

```bash
# Step 1: Process MAF files
python scripts/process_maf_data.py

# Step 2: Generate clinical reports
python scripts/generate_clinical_report.py
```

### `process-maf-data.py`

Please run this script from either the project root, use `scripts/process-maf-data.py` in the command-line or from within `scripts/` directory.

The output of this script will be in `data/results/processed_data`, and will contain: (a) `**_analysis.json` files and (b) all `processed_data/*.csv` files.

To use the script: 

```bash
> python scripts/process-maf-data.py
2025-01-14 11:55:36,756 - INFO - Processing MMCID-26B
2025-01-14 11:55:37,978 - INFO - Processing MMCID-30B
2025-01-14 11:55:38,683 - INFO - Processing TCMK1-14B
2025-01-14 11:55:39,364 - INFO - Processing TMCK1-23B
2025-01-14 11:55:43,944 - INFO - 
Analysis Summary:
2025-01-14 11:55:43,944 - INFO - 
MMCID-26B:
2025-01-14 11:55:43,944 - INFO -   Total variants: 535
2025-01-14 11:55:43,944 - INFO -   Variant types: 5
2025-01-14 11:55:43,944 - INFO -   Mean VAF: 0.164
2025-01-14 11:55:43,944 - INFO -   Total genes affected: 438
# ...
```

```bash
# Ensure your data structure is:
data/
  results/
    processed_data/
      samplename_analysis.json
      samplename_complete.csv
      samplename_del.csv
      samplename_dnp.csv
      samplename_ins.csv
      samplename_snp.csv
      samplename_tnp.csv
      ...
```
The file structure here is important for the `process-clincal-report.py` script to run properly. All of the `.csv` files are from the main `.json` file. If you wish to change anything about the data processing, please change the first script `process-maf-data.py` first before changing the clincal report script. It's important to maintain the structure within the main `.json` file.

### `process-clinical-report.py`

Within the `data/results/reports` directory, the script will generate a summary report and visualizations for each sample. The output will be in the `data/results/reports` directory.

The required `.json` and `.csv` files will be duplicated due to sensitivity and robustness for code changes and comparisons. This will be locatedd within the `data/results/reports/data` directory.

# Run the script

```bash
> python scripts/process-clinical-report.py
 2025-01-14 12:07:23,846 - INFO - Loading complete data from /Users/adamkurth/Documents/vscode/research/ugenome/ugenome-maf-pdf/data/results/processed_data/MMCID-26B_complete.csv
2025-01-14 12:07:24,078 - WARNING - Visualization directory not found at: data/results/visualizations
2025-01-14 12:07:24,084 - INFO - Generated clinical report for MMCID-26B
2025-01-14 12:07:24,084 - INFO - Loading complete data from /Users/adamkurth/Documents/vscode/research/ugenome/ugenome-maf-pdf/data/results/processed_data/MMCID-30B_complete.csv
2025-01-14 12:07:24,313 - WARNING - Visualization directory not found at: data/results/visualizations
2025-01-14 12:07:24,318 - INFO - Generated clinical report for MMCID-30B
2025-01-14 12:07:24,318 - INFO - Report generation complete. Data files copied to reports/data/
```

Tweaking the formatting of the report can be done within the `ClinicalReportGenerator` class in the `scripts/process-clinical-report.py` script. 

## Classes and Methods

### `ClinicalReportGenerator`

- `create_variant_table(self, df: pd.DataFrame, variant_type: str) -> Table`
  - Create an enhanced table for variants with modern styling.

- `create_summary_section(self, analysis_data: dict) -> list`
  - Create summary section with modern styling.

### `GeneVisualizer`

- Class for generating gene-specific visualizations.

## Dependencies

- pandas
- json
- shutil
- pathlib
- logging
- datetime
- matplotlib
- seaborn
- reportlab
- PyPDF2
