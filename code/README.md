# Ploidy Assay Analysis Code

This directory contains the Python scripts used for analyzing ploidy determination methods.

## Scripts

### ploidy-comparison-with-ambiguous-gate.py

A  script that compares ploidy results between flow cytometry and benomyl assay methods, with special attention to the ambiguous (G1 diploid and G2 haploid) gate in flow cytometry data.

#### Key Features

- Processes flow cytometry data from multiple plate formats (WGS1, WGS2, WGS hybrid)
- Calculates and tracks cells that fall into the "ambiguous gate" (G1 diploid or G2 diploid)
- Compares results between benomyl and flow cytometry methods
- Generates detailed statistics on agreement between methods
- Creates specialized output files for disagreements and flow cytometry results
- Produces summary statistics on ambiguous gate percentages

#### Usage

```bash
python3 ploidy-comparison-with-ambiguous-gate.py
```

#### Functions

- `format_well_id(plate, row, column)`: Formats well ID in a standardized way (e.g., "WGS1A1")
- `extract_plate_row_col(well_id)`: Extracts plate, row, and column from a well ID
- `load_benomyl_data()`: Loads the benomyl assay data from CSV
- `extract_flow_data(file_path)`: Processes flow cytometry data from Excel files with ambiguous gate calculation
- `load_flow_cytometry_data()`: Loads and combines all flow cytometry data
- `compare_methods(benomyl_data, flow_data)`: Compares the results of both methods
- `generate_summary_statistics(comparison_results)`: Generates comprehensive statistics including ambiguous gate metrics
- `main()`: Orchestrates the entire comparison workflow

## General Functions and Features

- `format_well_id(plate, row, column)`: Formats well ID in a standardized way (e.g., "WGS1A1")
- `extract_plate_row_col(well_id)`: Extracts plate, row, and column from a well ID
- `load_benomyl_data()`: Loads the benomyl assay data from CSV
- `extract_flow_data(file_path)`: Processes flow cytometry data from Excel files
- `load_flow_cytometry_data()`: Loads and combines all flow cytometry data
- `compare_methods(benomyl_data, flow_data)`: Compares the results of both methods
- `analyze_discrepancies(comparison_results)`: Analyzes patterns in method disagreements
- `create_discrepancy_visualizations(...)`: Creates visualizations for discrepancy analysis
- `generate_summary_statistics(...)`: Generates comprehensive statistics

## Dependencies

This code requires the following Python packages:

- pandas: For data manipulation
- numpy: For numerical operations
- openpyxl: For reading Excel files
- matplotlib: For creating visualizations
- scipy: For statistical analysis
- re: For regular expression operations

## Path Configuration

The scripts use the following path configuration:

```python
BENOMYL_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/raw_data/benamil_assay'
FLOW_CYTOMETRY_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/raw_data/flow_cytometry'
OUTPUT_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/analysis_output'
FIGURES_PATH = os.path.join(OUTPUT_PATH, 'figures')  # Used in ploidy_discrepancy_analysis.py
CODE_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/code'  # Used in ploidy-comparison-with-ambiguous-gate.py
```

To use the scripts in a different environment, update these paths to match your file system.

## Quality Filtering

The scripts apply filtering to flow cytometry data:
- Excludes samples with fewer than 5,000 events
- Flags potentially mixed populations based on cell distribution

## Confidence Metrics

The scripts calculate confidence scores for each method:
- **Benomyl confidence**: 0.5 for borderline results, 1.0 for clear results
- **Flow cytometry confidence**: Based on the percentage of cells in the dominant gate (haploid or diploid)

## Ambiguous Gate Analysis

The ploidy-comparison-with-ambiguous-gate.py script specifically focuses on calculating and analyzing the "ambiguous gate" in flow cytometry data:
- Tracks cells that don't fall into either haploid or diploid gates
- Calculates percentages of ambiguous cells relative to total singlets
- Provides summary statistics on ambiguous gate percentages (mean, median, min, max, standard deviation)
- Helps identify potential mixed populations or technical issues in flow cytometry analysis
