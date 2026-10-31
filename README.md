# Cancer staging of miRNA expression

Group per-sample miRNA expression measurements by TCGA clinical tumor stage. The included `reader.py` links expression sample barcodes to patient annotations and exports expression values labeled as normal, stage I, II, III, or IV.

The bundled example uses **hsa-mir-490** in **LUAD** (lung adenocarcinoma) and **LUSC** (lung squamous cell carcinoma). This is a data preparation workflow; it does not infer a patient's stage from miRNA expression or perform statistical testing.

## Files

| File | Purpose |
| --- | --- |
| `reader.py` | Reads expression CSVs, matches clinical annotations, and writes stage-grouped CSVs. |
| `TCGA_clinical_data.tsv` | Clinical annotations: 11,315 rows; the script uses `submitter_id` and `tumor_stage`. |
| `TCGA_biospecimen_sample_sheet.tsv` | Bundled biospecimen metadata: 49,360 rows; not read by the script. |
| `TCGA_clinical_sample_sheet.tsv` | Bundled clinical file metadata: 12,444 rows; not read by the script. |
| `hsa-mir-490-output_LUAD.csv` | LUAD expression input: 567 rows. |
| `hsa-mir-490-output_LUSC.csv` | LUSC expression input: 523 rows. |
| `staging-hsa-mir-490_LUAD.csv` | Included LUAD stage-grouped output: 560 rows. |
| `staging-hsa-mir-490_LUSC.csv` | Included LUSC stage-grouped output: 519 rows. |

## Requirements

Use Python 3 with pandas, NumPy, and SciPy. All three packages are imported by the current script, although NumPy and SciPy are not used in the processing logic.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install pandas numpy scipy
```

Package versions are not pinned in the original project.

## Usage

1. Open `reader.py` and change `pathname` to the absolute path of your local project folder. The checked-in script retains its original path:

   ```python
   pathname = '/Users/darmirador/Desktop/Python Scripts/cancer_staging_miRNA'
   ```

2. Keep `TCGA_clinical_data.tsv` and the expression input files in that folder.
3. Run the script, supplying the miRNA identifier exactly as it appears in the input filenames:

   ```bash
   python3 reader.py -i hsa-mir-490
   ```

The script prints progress for each cancer type and writes output into `pathname`. Existing outputs for the selected miRNA and discovered cancer types are overwritten. The usage message inside the original script refers to `cancer_staging_reader_miRNA.py`; the actual bundled filename is `reader.py`.

## Input format

Expression filenames must follow `<miRNA>-output_<CANCER>.csv`, for example `hsa-mir-490-output_LUAD.csv`. The parser uses underscores to extract the cancer type, so keep the naming pattern simple.

Required expression columns:

| Column | Meaning |
| --- | --- |
| `tissue` | `normal` or `cancer` (case-insensitive). Other values are omitted. |
| `barcode` | TCGA sample barcode; its first 12 characters are used as the patient identifier. |
| `rpm` | Expression value copied unchanged to the output. |

The clinical TSV requires `submitter_id` and `tumor_stage`. The script does not download data or calculate/normalize expression values.

## Stage assignment

- Samples labeled `normal` are retained as normal without requiring a clinical match.
- For `cancer` samples, the first 12 barcode characters are matched against clinical `submitter_id` values.
- If multiple clinical records match, the first occurrence is used.
- Stage annotations are converted to lowercase and checked for the substrings `stage iv`, `stage iii`, `stage ii`, then `stage i`. Substages such as stage IIA are therefore grouped with stage II.
- Cancer samples with no matching patient or no recognized stage are omitted. Missing stages, stage 0, and unrecognized labels have no separate output group.
- Output rows are ordered normal, stage I, stage II, stage III, then stage IV. Input order is preserved within each group.

## Output format and included results

Each output is named `staging-<miRNA>_<CANCER>.csv` and contains `barcode,stage,rpm`. Stage values are lowercase: `normal`, `stage i`, `stage ii`, `stage iii`, and `stage iv`.

| Group | LUAD | LUSC |
| --- | ---: | ---: |
| Normal | 46 | 45 |
| Stage I | 283 | 230 |
| Stage II | 123 | 158 |
| Stage III | 84 | 80 |
| Stage IV | 24 | 6 |
| **Total output rows** | **560** | **519** |
| Input rows absent from output | 7 | 4 |

These counts describe the bundled CSV files, not a new execution of the script. The script does not produce a report explaining individual omissions.

## Scope and reproducibility

This repository preserves the supplied script and data files. It contains no model training, survival analysis, stage comparisons, plots, automated test suite, or dependency lockfile. The clinical table includes annotations beyond the two example cancer types; the script matches by patient ID rather than filtering by project.

The filenames identify the annotation tables as TCGA data, but download commands, original source URLs, release identifiers, and upstream expression-generation code are not included. Confirm provenance and applicable data-use terms before redistributing or extending the dataset. No software license was supplied with the original folder.
