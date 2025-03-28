# Ploidy Assay Analysis Project

This repository contains code and data for analyzing and comparing different methods of ploidy determination in *S. cerevisiae*. The project focuses on comparing benomyl sensitivity assays with flow cytometry data to evaluate their agreement in ploidy classification.

## Project Structure

```
/your/path/here/GitHub/ploidy_assay/
├── raw_data/
│   ├── benamil_assay/
│   │   └── SM_benomyl_assay_df.csv
│   └── flow_cytometry/
│       ├── WGS1_ploidy_assay.xlsx
│       ├── WGS2_ploid_assay.xlsx
│       └── WGShybrid_ploidy_assay.xlsx
├── analysis_output/
│   ├── figures/
│   ├── flow_cytometry_with_ambiguous_gate.csv
│   ├── ploidy_comparison_summary.txt
│   ├── ploidy_comparison_summary.csv
│   ├── ploidy_comparison_with_ambiguous_gate.csv
└── code/
    └── ploidy_comparison_with_ambiguous_gate.py
```

## Overview

This project analyzes two different methods of determining ploidy in yeast:

1. **Benomyl Sensitivity Assay**: Growth inhibition in the presence of benomyl indicates diploid cells, while resistance suggests haploid cells.

2. **Flow Cytometry**: DNA content measurement using SYBR Green staining to classify cells as haploid or diploid based on fluorescence intensity.

The analysis includes calculating the number of cells in the "ambiguous gate" of flow cytometry data (cells that don't clearly fall into haploid or diploid categories) and a detailed analysis of discrepancies between the two methods.

## Key Findings

- Agreement between methods: ~50%
-Flow cytometry method could be revised to fix cells at a certain growth stage (i.e, prevent G2 stage)
- Primary sources of disagreement: 
  - Borderline benomyl results
  - Clean disagreements (when both methods show high confidence)

## Getting Started

To run the analysis:

1. Ensure that Python 3.6+ and required libraries are installed
2. Clone this repository
3. Run the analysis script:

```bash
cd /your/path/here/GitHub/ploidy_assay/code
python3 ploidy_discrepancy_analysis.py
```

## Required Libraries

- pandas
- numpy
- openpyxl
- matplotlib
- scipy

## See Also

- The `code/` directory contains detailed information about the analysis scripts
- The `analysis_output/` directory contains documentation on the output files and figures
