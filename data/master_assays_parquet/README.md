# Master Primer Assays (Atlas-Screened)

Converted from `master_primer_assays_atlas_screened.xlsx` (120 MB) to Parquet
format to fit within GitHub's file size limit. All data preserved.

## Files

| File | Rows | Size | Description |
|------|------|------|-------------|
| `All_Primers.parquet` | 186,848 | 14.8 MB | All target-detecting screened primer pairs across 82 cell types |
| `Best_Per_DMR.parquet` | 40,190 | 4.0 MB | Best target-detecting primer per DMR (ranked by penalty) |
| `DMR_Blocks.parquet` | 46,977 | 1.2 MB | DMR block metadata (coordinates, DMS scores, genes) |
| `Selected_Panel.parquet` | 519 | 127 KB | Final selected panel (target-detecting, ePCR PASS, <=20 per cell type) |
| `Summary.parquet` | 82 | 9 KB | Per-cell-type summary statistics |
| `Data_Dictionary.parquet` | 81 | 13 KB | Column definitions |

## Reading the Data

```python
import pandas as pd

# Read a single table
df = pd.read_parquet('data/master_assays_parquet/All_Primers.parquet')

# Convert back to Excel
all_primers = pd.read_parquet('data/master_assays_parquet/All_Primers.parquet')
all_primers.to_excel('master_primer_assays.xlsx', index=False)
```

## Requirements

```bash
pip install pandas pyarrow
```
