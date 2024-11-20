# ugenome

## `scripts/`

1. `process-maf-data.py` - Process MAF files and generate a summary report.
   
2. `process-clinical-report.py` - Process clinical report and generate a summary report.


### `process-maf-data.py`

Place this script in the root of the repo. 

Ensure `data/results` directory contains: (a) `complete_analysis.json` and (b) `processed_data/*.csv` files.

Use the script: 

```bash
> python scripts/process-maf-data.py

2024-11-20 10:28:44,493 - INFO - Analysis complete. Processed 4 MAF files.
2024-11-20 10:28:44,493 - INFO - Analysis Summary:
2024-11-20 10:28:44,493 - INFO - MMCID-26B: 535 variants processed
2024-11-20 10:28:44,493 - INFO - MMCID-30B: 275 variants processed
2024-11-20 10:28:44,493 - INFO - TCMK1-14B: 327 variants processed
2024-11-20 10:28:44,493 - INFO - TMCK1-23B: 322 variants processed
```

```bash
# Ensure your data structure is:
data/
  results/
    complete_analysis.json
    processed_data/
      sample_name_variant_type.csv
      ...

# Run the script
python generate_clinical_report.py

# Check the output in:
data/results/reports/
```