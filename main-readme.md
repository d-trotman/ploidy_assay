# Ploidy Assay Analysis Project

This repository contains code and data for analyzing and comparing different methods of ploidy determination in *S. cerevisiae*. The project focuses on comparing benomyl sensitivity assays with flow cytometry data to evaluate their agreement in ploidy classification.

## Project Structure

```
/Users/dawsontrotman/Documents/GitHub/ploidy_assay/
├── raw_data/
│   ├── benamil_assay/
│   │   └── SM_benomyl_assay_df.csv
│   └── flow_cytometry/
│       ├── WGS1_ploidy_assay.xlsx
│       ├── WGS2_ploid_assay.xlsx
│       └── WGShybrid_ploidy_assay.xlsx
├── analysis_output/
│   ├── figures/
│   ├── ploidy_method_comparison.csv
│   ├── ploidy_method_discrepancy_summary.txt
│   ├── discrepancy_borderline_benomyl.csv
│   ├── discrepancy_mixed_flow.csv
│   └── discrepancy_clean.csv
└── code/
    └── ploidy_discrepancy_analysis.py
```

## Overview

This project analyzes two different methods of determining ploidy in yeast:

1. **Benomyl Sensitivity Assay**: Growth inhibition in the presence of benomyl indicates diploid cells, while resistance suggests haploid cells.

2. **Flow Cytometry**: DNA content measurement using SYBR Green staining to classify cells as haploid or diploid based on fluorescence intensity.

The analysis includes calculating the number of cells in the "ambiguous gate" of flow cytometry data (cells that don't clearly fall into haploid or diploid categories) and a detailed analysis of discrepancies between the two methods.

## Key Findings

- Agreement between methods: ~50%
- Primary sources of disagreement: 
  - Borderline benomyl results
  - Mixed populations in flow cytometry
  - Clean disagreements (when both methods show high confidence)

## Getting Started

To run the analysis:

1. Ensure that Python 3.6+ and required libraries are installed
2. Clone this repository
3. Run the analysis script:

```bash
cd /Users/dawsontrotman/Documents/GitHub/ploidy_assay/code
python ploidy_discrepancy_analysis.py
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
