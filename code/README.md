# Ploidy Assay Analysis Code

This directory contains the Python scripts used for analyzing ploidy determination methods.

## Scripts

### ploidy_discrepancy_analysis.py

The main analysis script that compares benomyl assay results with flow cytometry data, calculates the ambiguous gate in flow cytometry, and performs detailed analysis of discrepancies between methods.

#### Key Features

- Calculates the number of cells in the "ambiguous gate" (cells that don't fall into either haploid or diploid gates)
- Classifies disagreement types and provides potential explanations
- Identifies patterns in disagreements by plate, row, or column
- Generates visualizations to help understand the patterns
- Creates comprehensive summary reports and specialized output files

#### Usage

```bash
python ploidy_discrepancy_analysis.py
```

#### Functions

- `format_well_id(plate, row, column)`: Formats well ID in a standardized way (e.g., "WGS1A1")
- `extract_plate_row_col(well_id)`: Extracts plate, row, and column from a well ID
- `load_benomyl_data()`: Loads the benomyl assay data from CSV
- `extract_flow_data(file_path)`: Processes flow cytometry data from Excel files
- `load_flow_cytometry_data()`: Loads and combines all flow cytometry data
- `compare_methods(benomyl_data, flow_data)`: Compares the results of both methods
- `analyze_discrepancies(comparison_results)`: Analyzes patterns in method disagreements
- `create_discrepancy_visualizations(...)`: Creates visualizations for discrepancy analysis
- `generate_summary_statistics(...)`: Generates comprehensive statistics
- `main()`: Orchestrates the entire analysis workflow

## Dependencies

This code requires the following Python packages:

- pandas: For data manipulation
- numpy: For numerical operations
- openpyxl: For reading Excel files
- matplotlib: For creating visualizations
- scipy: For statistical analysis

## Path Configuration

The script uses the following path configuration:

```python
BENOMYL_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/raw_data/benamil_assay'
FLOW_CYTOMETRY_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/raw_data/flow_cytometry'
OUTPUT_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/analysis_output'
FIGURES_PATH = os.path.join(OUTPUT_PATH, 'figures')
```

To use the script in a different environment, update these paths to match your file system.

## Quality Filtering

The script applies filtering to flow cytometry data:
- Excludes samples with fewer than 5,000 events
- Flags potentially mixed populations based on cell distribution

## Confidence Metrics

The script calculates confidence scores for each method:
- **Benomyl confidence**: 0.5 for borderline results, 1.0 for clear results
- **Flow cytometry confidence**: Based on the percentage of cells in the dominant gate (haploid or diploid)
