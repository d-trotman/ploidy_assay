#!/usr/bin/env python3
"""
Ploidy Data Mapping Validation Script

This script validates the mapping between benomyl assay and flow cytometry data,
performing multiple cross-checks to ensure data is correctly aligned by well ID.

Validation methods:
1. Sample ID format validation
2. Cross-reference check between datasets
3. Visual mapping representation

Usage:
    python ploidy_mapping_validation.py

Output:
    CSV files and visualization of data mapping between assays
"""

import os
import pandas as pd
import numpy as np
import openpyxl
import re
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec

# Define paths (same as original script)
BENOMYL_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/raw_data/benamil_assay'
FLOW_CYTOMETRY_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/raw_data/flow_cytometry'
OUTPUT_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/analysis_output'
VALIDATION_PATH = os.path.join(OUTPUT_PATH, 'validation')

# Create output directories if they don't exist
os.makedirs(OUTPUT_PATH, exist_ok=True)
os.makedirs(VALIDATION_PATH, exist_ok=True)

def format_well_id(plate, row, column):
    """Format well ID from plate, row, and column"""
    return f"WGS{plate}{row}{column}"

def extract_plate_row_col(well_id):
    """Extract plate, row, and column from well ID
    
    Args:
        well_id: Sample ID in format WGSxYZ (x=plate, Y=row, Z=column)
        
    Returns:
        tuple: (plate, row, column) or (None, None, None) if invalid format
    """
    if not isinstance(well_id, str):
        return None, None, None
        
    match = re.match(r'WGS(\d+)([A-Z]+)(\d+)', well_id)
    if match:
        try:
            plate = int(match.group(1))
            row = match.group(2)
            column = int(match.group(3))
            return plate, row, column
        except (ValueError, TypeError):
            return None, None, None
    return None, None, None

def load_benomyl_data():
    """Load benomyl assay data
    
    Returns:
        benomyl_data: Dictionary of benomyl data keyed by well ID
        benomyl_df: Raw DataFrame of benomyl data
    """
    print("Loading benomyl assay data...")
    
    benomyl_file = os.path.join(BENOMYL_PATH, 'SM_benomyl_assay_df.csv')
    benomyl_df = pd.read_csv(benomyl_file)
    
    # Save a copy of the raw data
    benomyl_raw_file = os.path.join(VALIDATION_PATH, 'raw_benomyl_data.csv')
    benomyl_df.to_csv(benomyl_raw_file, index=False)
    
    # Create a dictionary for easy lookup
    benomyl_data = {}
    invalid_format_count = 0
    
    for _, row in benomyl_df.iterrows():
        # Validate data format first
        try:
            plate = row['plate']
            row_letter = row['row']
            column = row['column']
            
            # Check for non-standard formats
            if not isinstance(plate, (int, float)) or pd.isna(plate):
                invalid_format_count += 1
                continue
                
            if not isinstance(row_letter, str) or not row_letter.isalpha():
                invalid_format_count += 1
                continue
                
            if not isinstance(column, (int, float)) or pd.isna(column):
                invalid_format_count += 1
                continue
            
            # Format the well ID
            well_id = format_well_id(int(plate), row_letter, int(column))
            
            benomyl_data[well_id] = {
                'is_diploid': row['growth inhibited (diploid)'],
                'is_borderline': row['borderline'],
                'notes': row['notes'] if not pd.isna(row['notes']) else "",
                'double_checked': row['doublechecked'] if not pd.isna(row['doublechecked']) else False,
                'raw_plate': plate,
                'raw_row': row_letter,
                'raw_column': column
            }
        except Exception as e:
            print(f"  Error processing benomyl row: {e}")
            invalid_format_count += 1
    
    print(f"  Loaded {len(benomyl_data)} benomyl data points")
    if invalid_format_count > 0:
        print(f"  WARNING: {invalid_format_count} benomyl entries had invalid format and were skipped")
    
    return benomyl_data, benomyl_df

def extract_flow_data(file_path):
    """Extract flow cytometry data with format validation
    
    Returns:
        flow_data: Dictionary of flow cytometry data keyed by well ID
        raw_df: DataFrame with raw data for inspection
    """
    file_name = os.path.basename(file_path)
    print(f"Processing flow cytometry data from {file_name}...")
    
    flow_data = {}
    excluded_count = 0
    invalid_format_count = 0
    
    try:
        # Load the Excel file
        wb = openpyxl.load_workbook(file_path, data_only=True)
        sheet = wb.active
        
        # Find column indices (assume header is in row 1)
        header_row = [cell.value for cell in sheet[1]]
        
        # Common columns in standardized format
        sample_name_col = header_row.index('sample name') + 1 if 'sample name' in header_row else -1
        
        # Check if required columns exist
        required_columns = ['All Events Count', 'singlets Count', 'haploid Count', 
                          'diploid 2x Count', 'haploid % Parent', 'diploid 2x % Parent']
        
        for col in required_columns:
            if col not in header_row:
                print(f"  ERROR: Required column '{col}' not found in {file_name}")
                return {}, pd.DataFrame()
        
        all_events_count_col = header_row.index('All Events Count') + 1
        singlets_count_col = header_row.index('singlets Count') + 1
        
        # Find hap/dip column if it exists (WGS2 has it, others don't)
        ploidy_col = -1
        if 'hap/dip' in header_row:
            ploidy_col = header_row.index('hap/dip') + 1
        elif 'dip_hap_na' in header_row:
            ploidy_col = header_row.index('dip_hap_na') + 1
        
        # Common columns in all files
        haploid_count_col = header_row.index('haploid Count') + 1
        diploid_count_col = header_row.index('diploid 2x Count') + 1
        haploid_percent_col = header_row.index('haploid % Parent') + 1
        diploid_percent_col = header_row.index('diploid 2x % Parent') + 1
        
        # Create a list to store raw data for inspection
        raw_data_rows = []
        
        # Extract data starting from row 2
        for row_idx, row in enumerate(sheet.iter_rows(min_row=2), start=2):
            # Get sample ID based on file format
            if sample_name_col > 0:
                sample_id = row[sample_name_col - 1].value
            else:
                # Fall back to alternative sample ID location if needed
                sample_id = None
            
            # Store raw data for this row
            raw_data_rows.append({
                'row_number': row_idx,
                'sample_id': sample_id,
                'all_events_count': row[all_events_count_col - 1].value if all_events_count_col > 0 else None,
                'file_name': file_name
            })
            
            # Validate sample ID format
            if not sample_id or not isinstance(sample_id, str) or not sample_id.startswith('WGS'):
                invalid_format_count += 1
                continue
            
            plate, row_letter, column = extract_plate_row_col(sample_id)
            if plate is None or row_letter is None or column is None:
                invalid_format_count += 1
                continue
            
            try:
                # Get event counts
                all_events_count = row[all_events_count_col - 1].value or 0
                if isinstance(all_events_count, str):
                    all_events_count = 0
                
                if all_events_count < 5000:
                    excluded_count += 1
                    continue  # Skip samples with All Events Count < 5000
                
                singlets_count = row[singlets_count_col - 1].value or 0
                haploid_count = row[haploid_count_col - 1].value or 0
                diploid_count = row[diploid_count_col - 1].value or 0
                
                # Handle potential string values
                if isinstance(singlets_count, str): singlets_count = 0
                if isinstance(haploid_count, str): haploid_count = 0
                if isinstance(diploid_count, str): diploid_count = 0
                
                # Calculate ambiguous gate
                ambiguous_count = singlets_count - (haploid_count + diploid_count)
                
                # Calculate percentages
                haploid_percent = (haploid_count / singlets_count * 100) if singlets_count > 0 else 0
                diploid_percent = (diploid_count / singlets_count * 100) if singlets_count > 0 else 0
                ambiguous_percent = (ambiguous_count / singlets_count * 100) if singlets_count > 0 else 0
                
                # Get declared ploidy if available
                declared_ploidy = None
                if ploidy_col > 0:
                    declared_ploidy = row[ploidy_col - 1].value
                
                # Determine ploidy from data
                if declared_ploidy:
                    is_haploid = (declared_ploidy == "haploid" or 
                                 (haploid_count > diploid_count and haploid_percent > diploid_percent))
                    is_diploid = (declared_ploidy == "diploid" or 
                                 (diploid_count > haploid_count and diploid_percent > haploid_percent))
                else:
                    is_haploid = (haploid_count > diploid_count and haploid_percent > diploid_percent)
                    is_diploid = (diploid_count > haploid_count and diploid_percent > haploid_percent)
                
                # Store data
                flow_data[sample_id] = {
                    'declared_ploidy': declared_ploidy,
                    'is_haploid': is_haploid,
                    'is_diploid': is_diploid,
                    'all_events_count': all_events_count,
                    'singlets_count': singlets_count,
                    'haploid_count': haploid_count,
                    'diploid_count': diploid_count,
                    'ambiguous_count': ambiguous_count,
                    'haploid_percent': haploid_percent,
                    'diploid_percent': diploid_percent,
                    'ambiguous_percent': ambiguous_percent,
                    'source_file': file_name,
                    'row_in_file': row_idx,
                    'raw_sample_id': sample_id,
                    'plate': plate,
                    'row': row_letter,
                    'column': column
                }
            except Exception as e:
                print(f"  Error processing row {row_idx} in {file_name}: {e}")
                invalid_format_count += 1
        
        print(f"  Extracted {len(flow_data)} samples from {file_name}")
        print(f"  Excluded {excluded_count} samples with All Events Count < 5000")
        print(f"  Skipped {invalid_format_count} samples with invalid format")
        
    except Exception as e:
        print(f"  ERROR processing file {file_name}: {e}")
        return {}, pd.DataFrame()
    
    # Create a DataFrame with raw data for inspection
    raw_df = pd.DataFrame(raw_data_rows)
    
    # Save raw data for inspection
    raw_file = os.path.join(VALIDATION_PATH, f'raw_flow_data_{os.path.basename(file_path).split(".")[0]}.csv')
    raw_df.to_csv(raw_file, index=False)
    
    return flow_data, raw_df

def load_flow_cytometry_data():
    """Load all flow cytometry data
    
    Returns:
        combined_data: Dictionary of flow cytometry data from all files keyed by well ID
        combined_raw: Combined DataFrame of raw flow cytometry data for inspection
    """
    print("Loading flow cytometry data...")
    
    wgs1_file = os.path.join(FLOW_CYTOMETRY_PATH, 'WGS1_ploidy_assay.xlsx')
    wgs2_file = os.path.join(FLOW_CYTOMETRY_PATH, 'WGS2_ploid_assay.xlsx')
    wgs_hybrid_file = os.path.join(FLOW_CYTOMETRY_PATH, 'WGShybrid_ploidy_assay.xlsx')
    
    # Extract data from each file
    wgs1_data, wgs1_raw = extract_flow_data(wgs1_file)
    wgs2_data, wgs2_raw = extract_flow_data(wgs2_file)
    wgs_hybrid_data, wgs_hybrid_raw = extract_flow_data(wgs_hybrid_file)
    
    # Combine all data
    combined_data = {**wgs1_data, **wgs2_data, **wgs_hybrid_data}
    
    # Combine raw data
    dfs = []
    if not wgs1_raw.empty: dfs.append(wgs1_raw)
    if not wgs2_raw.empty: dfs.append(wgs2_raw)
    if not wgs_hybrid_raw.empty: dfs.append(wgs_hybrid_raw)
    
    combined_raw = pd.concat(dfs) if dfs else pd.DataFrame()
    
    # Save combined raw data
    combined_raw_file = os.path.join(VALIDATION_PATH, 'combined_raw_flow_data.csv')
    combined_raw.to_csv(combined_raw_file, index=False)
    
    print(f"  Combined {len(combined_data)} flow cytometry data points")
    
    return combined_data, combined_raw

def validate_sample_id_formats(benomyl_data, flow_data):
    """Validate sample ID formats between datasets
    
    Returns:
        benomyl_format_df: DataFrame with benomyl well ID format analysis
        flow_format_df: DataFrame with flow cytometry well ID format analysis
        benomyl_plate_counts: Series with count of samples per plate in benomyl data
        flow_plate_counts: Series with count of samples per plate in flow cytometry data
    """
    print("Validating sample ID formats...")
    
    # Extract all well IDs
    benomyl_well_ids = list(benomyl_data.keys())
    flow_well_ids = list(flow_data.keys())
    
    # Parse sample IDs to check format consistency
    benomyl_format_check = []
    for well_id in benomyl_well_ids:
        plate, row, column = extract_plate_row_col(well_id)
        benomyl_format_check.append({
            'well_id': well_id,
            'plate': plate,
            'row': row,
            'column': column,
            'is_valid': (plate is not None and row is not None and column is not None)
        })
    
    flow_format_check = []
    for well_id in flow_well_ids:
        plate, row, column = extract_plate_row_col(well_id)
        flow_format_check.append({
            'well_id': well_id,
            'plate': plate,
            'row': row,
            'column': column,
            'is_valid': (plate is not None and row is not None and column is not None)
        })
    
    # Create DataFrames for analysis
    benomyl_format_df = pd.DataFrame(benomyl_format_check)
    flow_format_df = pd.DataFrame(flow_format_check)
    
    # Calculate statistics
    benomyl_valid_count = benomyl_format_df['is_valid'].sum()
    flow_valid_count = flow_format_df['is_valid'].sum()
    
    benomyl_invalid = benomyl_format_df[~benomyl_format_df['is_valid']]
    flow_invalid = flow_format_df[~flow_format_df['is_valid']]
    
    print(f"  Benomyl: {benomyl_valid_count}/{len(benomyl_well_ids)} valid sample IDs")
    print(f"  Flow Cytometry: {flow_valid_count}/{len(flow_well_ids)} valid sample IDs")
    
    # Save invalid IDs for inspection
    if not benomyl_invalid.empty:
        invalid_benomyl_file = os.path.join(VALIDATION_PATH, 'invalid_benomyl_ids.csv')
        benomyl_invalid.to_csv(invalid_benomyl_file, index=False)
        print(f"  Invalid benomyl IDs saved to {invalid_benomyl_file}")
    
    if not flow_invalid.empty:
        invalid_flow_file = os.path.join(VALIDATION_PATH, 'invalid_flow_ids.csv')
        flow_invalid.to_csv(invalid_flow_file, index=False)
        print(f"  Invalid flow IDs saved to {invalid_flow_file}")
    
    # Save all ID format data
    benomyl_format_file = os.path.join(VALIDATION_PATH, 'benomyl_id_formats.csv')
    flow_format_file = os.path.join(VALIDATION_PATH, 'flow_id_formats.csv')
    
    benomyl_format_df.to_csv(benomyl_format_file, index=False)
    flow_format_df.to_csv(flow_format_file, index=False)
    
    # Summary of plate distributions
    benomyl_plate_counts = benomyl_format_df[benomyl_format_df['is_valid']]['plate'].value_counts().sort_index()
    flow_plate_counts = flow_format_df[flow_format_df['is_valid']]['plate'].value_counts().sort_index()
    
    print("\n  Plate distribution in benomyl data:")
    for plate, count in benomyl_plate_counts.items():
        print(f"    Plate {plate}: {count} samples")
    
    print("\n  Plate distribution in flow cytometry data:")
    for plate, count in flow_plate_counts.items():
        print(f"    Plate {plate}: {count} samples")
    
    # Create a mapping distribution check and save results
    format_summary_file = os.path.join(VALIDATION_PATH, 'id_format_validation_summary.txt')
    with open(format_summary_file, 'w') as f:
        f.write("Sample ID Format Validation Summary\n")
        f.write("==================================\n\n")
        f.write(f"Benomyl: {benomyl_valid_count}/{len(benomyl_well_ids)} valid sample IDs\n")
        f.write(f"Flow Cytometry: {flow_valid_count}/{len(flow_well_ids)} valid sample IDs\n\n")
        
        f.write("Plate distribution in benomyl data:\n")
        for plate, count in benomyl_plate_counts.items():
            f.write(f"  Plate {plate}: {count} samples\n")
        
        f.write("\nPlate distribution in flow cytometry data:\n")
        for plate, count in flow_plate_counts.items():
            f.write(f"  Plate {plate}: {count} samples\n")
    
    print(f"  Format validation summary saved to {format_summary_file}")
    
    return benomyl_format_df, flow_format_df, benomyl_plate_counts, flow_plate_counts

def cross_reference_datasets(benomyl_data, flow_data):
    """Cross-reference benomyl and flow cytometry datasets
    
    Identifies common and unique well IDs between datasets
    
    Returns:
        reference_results: Dictionary with cross-reference analysis results
    """
    print("Cross-referencing benomyl and flow cytometry datasets...")
    
    benomyl_wells = set(benomyl_data.keys())
    flow_wells = set(flow_data.keys())
    
    # Find common and unique well IDs
    common_wells = benomyl_wells.intersection(flow_wells)
    only_benomyl = benomyl_wells - flow_wells
    only_flow = flow_wells - benomyl_wells
    
    print(f"  Common well IDs: {len(common_wells)}")
    print(f"  Wells only in benomyl: {len(only_benomyl)}")
    print(f"  Wells only in flow cytometry: {len(only_flow)}")
    
    # Group by plate for detailed analysis
    common_by_plate = {}
    benomyl_only_by_plate = {}
    flow_only_by_plate = {}
    
    # Process common wells
    for well in common_wells:
        plate, _, _ = extract_plate_row_col(well)
        if plate is None:
            continue
        if plate not in common_by_plate:
            common_by_plate[plate] = []
        common_by_plate[plate].append(well)
    
    # Process benomyl-only wells
    for well in only_benomyl:
        plate, _, _ = extract_plate_row_col(well)
        if plate is None:
            continue
        if plate not in benomyl_only_by_plate:
            benomyl_only_by_plate[plate] = []
        benomyl_only_by_plate[plate].append(well)
    
    # Process flow-only wells
    for well in only_flow:
        plate, _, _ = extract_plate_row_col(well)
        if plate is None:
            continue
        if plate not in flow_only_by_plate:
            flow_only_by_plate[plate] = []
        flow_only_by_plate[plate].append(well)
    
    # Create detailed breakdown by plate
    plates = sorted(set(list(common_by_plate.keys()) + 
                      list(benomyl_only_by_plate.keys()) + 
                      list(flow_only_by_plate.keys())))
    
    plate_breakdown = []
    for plate in plates:
        common_count = len(common_by_plate.get(plate, []))
        benomyl_only_count = len(benomyl_only_by_plate.get(plate, []))
        flow_only_count = len(flow_only_by_plate.get(plate, []))
        total_plate = common_count + benomyl_only_count + flow_only_count
        
        plate_breakdown.append({
            'plate': plate,
            'common_count': common_count,
            'benomyl_only_count': benomyl_only_count,
            'flow_only_count': flow_only_count,
            'total_count': total_plate,
            'common_percent': (common_count / total_plate * 100) if total_plate > 0 else 0,
            'benomyl_only_percent': (benomyl_only_count / total_plate * 100) if total_plate > 0 else 0,
            'flow_only_percent': (flow_only_count / total_plate * 100) if total_plate > 0 else 0
        })
    
    # Save well ID lists
    common_wells_file = os.path.join(VALIDATION_PATH, 'common_wells.csv')
    benomyl_only_file = os.path.join(VALIDATION_PATH, 'benomyl_only_wells.csv')
    flow_only_file = os.path.join(VALIDATION_PATH, 'flow_only_wells.csv')
    
    pd.DataFrame({'well_id': list(common_wells)}).to_csv(common_wells_file, index=False)
    pd.DataFrame({'well_id': list(only_benomyl)}).to_csv(benomyl_only_file, index=False)
    pd.DataFrame({'well_id': list(only_flow)}).to_csv(flow_only_file, index=False)
    
    # Save plate breakdown
    plate_breakdown_df = pd.DataFrame(plate_breakdown)
    plate_breakdown_file = os.path.join(VALIDATION_PATH, 'plate_breakdown.csv')
    plate_breakdown_df.to_csv(plate_breakdown_file, index=False)
    
    # Create a summary report
    summary_file = os.path.join(VALIDATION_PATH, 'cross_reference_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("Dataset Cross-Reference Summary\n")
        f.write("==============================\n\n")
        f.write(f"Total unique well IDs across both datasets: {len(benomyl_wells.union(flow_wells))}\n")
        f.write(f"Common well IDs found in both datasets: {len(common_wells)} ({len(common_wells)/len(benomyl_wells.union(flow_wells))*100:.1f}%)\n")
        f.write(f"Well IDs only in benomyl dataset: {len(only_benomyl)} ({len(only_benomyl)/len(benomyl_wells)*100:.1f}% of benomyl data)\n")
        f.write(f"Well IDs only in flow cytometry dataset: {len(only_flow)} ({len(only_flow)/len(flow_wells)*100:.1f}% of flow data)\n\n")
        
        f.write("Breakdown by plate:\n")
        for plate_info in plate_breakdown:
            f.write(f"\nPlate {plate_info['plate']}:\n")
            f.write(f"  Total wells: {plate_info['total_count']}\n")
            f.write(f"  Common wells: {plate_info['common_count']} ({plate_info['common_percent']:.1f}%)\n")
            f.write(f"  Benomyl only: {plate_info['benomyl_only_count']} ({plate_info['benomyl_only_percent']:.1f}%)\n")
            f.write(f"  Flow only: {plate_info['flow_only_count']} ({plate_info['flow_only_percent']:.1f}%)\n")
    
    print(f"  Cross-reference summary saved to {summary_file}")
    
    # Prepare results for return
    reference_results = {
        'common_wells': common_wells,
        'only_benomyl': only_benomyl,
        'only_flow': only_flow,
        'common_by_plate': common_by_plate,
        'benomyl_only_by_plate': benomyl_only_by_plate,
        'flow_only_by_plate': flow_only_by_plate,
        'plate_breakdown': plate_breakdown
    }
    
    return reference_results

def create_well_mapping_visualization(benomyl_data, flow_data, cross_ref_results):
    """Create visual mapping of well coverage between datasets
    
    Creates plate layout visualizations showing which wells have data in each dataset
    
    Args:
        benomyl_data: Dictionary of benomyl data
        flow_data: Dictionary of flow cytometry data
        cross_ref_results: Results from cross_reference_datasets
    """
    print("Creating well mapping visualization...")
    
    try:
        # Get all plates across both datasets
        benomyl_plates = set()
        flow_plates = set()
        
        for well_id in benomyl_data:
            plate, _, _ = extract_plate_row_col(well_id)
            if plate is not None:
                benomyl_plates.add(plate)
        
        for well_id in flow_data:
            plate, _, _ = extract_plate_row_col(well_id)
            if plate is not None:
                flow_plates.add(plate)
        
        all_plates = sorted(benomyl_plates.union(flow_plates))
        
        # Define standard plate dimensions
        # Assuming 8 rows (A-H) and 12 columns (1-12)
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        columns = list(range(1, 13))
        
        # Create a figure for each plate
        for plate in all_plates:
            fig = plt.figure(figsize=(15, 8))
            gs = GridSpec(1, 2, width_ratios=[1, 1])
            
            # Create data matrix for this plate
            # 0 = No data, 1 = Benomyl only, 2 = Flow only, 3 = Both
            data_matrix = np.zeros((len(rows), len(columns)))
            
            # Fill in the matrix
            for r_idx, row in enumerate(rows):
                for c_idx, col in enumerate(columns):
                    well_id = format_well_id(plate, row, col)
                    
                    in_benomyl = well_id in benomyl_data
                    in_flow = well_id in flow_data
                    
                    if in_benomyl and in_flow:
                        data_matrix[r_idx, c_idx] = 3  # Both
                    elif in_benomyl:
                        data_matrix[r_idx, c_idx] = 1  # Benomyl only
                    elif in_flow:
                        data_matrix[r_idx, c_idx] = 2  # Flow only
                    # 0 remains for no data
            
            # Create the heatmap
            ax1 = plt.subplot(gs[0])
            cmap = ListedColormap(['#f0f0f0', '#ffcccc', '#ccccff', '#99ff99'])
            sns.heatmap(data_matrix, ax=ax1, cmap=cmap, cbar=False,
                       xticklabels=columns, yticklabels=rows)
            
            # Add legend
            legend_elements = [
                mpatches.Patch(facecolor='#f0f0f0', label='No Data'),
                mpatches.Patch(facecolor='#ffcccc', label='Benomyl Only'),
                mpatches.Patch(facecolor='#ccccff', label='Flow Cytometry Only'),
                mpatches.Patch(facecolor='#99ff99', label='Both Datasets')
            ]
            ax1.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1))
            
            # Set title and labels
            ax1.set_title(f'Plate {plate} - Data Coverage Matrix')
            ax1.set_xlabel('Column')
            ax1.set_ylabel('Row')
            
            # Create pie chart showing distribution
            ax2 = plt.subplot(gs[1])
            
            # Count occurrences
            no_data = np.sum(data_matrix == 0)
            benomyl_only = np.sum(data_matrix == 1)
            flow_only = np.sum(data_matrix == 2)
            both = np.sum(data_matrix == 3)
            
            # Create pie chart
            sizes = [no_data, benomyl_only, flow_only, both]
            labels = ['No Data', 'Benomyl Only', 'Flow Only', 'Both Datasets']
            colors = ['#f0f0f0', '#ffcccc', '#ccccff', '#99ff99']
            
            # Only include non-zero segments
            non_zero_indices = [i for i, size in enumerate(sizes) if size > 0]
            if non_zero_indices:  # Check if we have any non-zero values
                ax2.pie([sizes[i] for i in non_zero_indices],
                       labels=[labels[i] for i in non_zero_indices],
                       colors=[colors[i] for i in non_zero_indices],
                       autopct='%1.1f%%', startangle=90)
                ax2.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
            
            ax2.set_title(f'Plate {plate} - Data Distribution')
            
            # Add plate statistics
            total_wells = len(rows) * len(columns)
            coverage = (benomyl_only + flow_only + both) / total_wells * 100
            matched = both / (benomyl_only + flow_only + both) * 100 if (benomyl_only + flow_only + both) > 0 else 0
            
            stat_text = f"""Plate {plate} Statistics:
            Total Wells: {total_wells}
            Overall Coverage: {coverage:.1f}%
            Wells in Both Datasets: {both}
            Wells Only in Benomyl: {benomyl_only}
            Wells Only in Flow: {flow_only}
            Matching Rate: {matched:.1f}%"""
            
            plt.figtext(0.5, 0.01, stat_text, ha='center', fontsize=10, bbox=dict(facecolor='white', alpha=0.5))
            
            # Adjust layout and save figure
            plt.tight_layout(rect=[0, 0.05, 1, 0.95])
            fig_file = os.path.join(VALIDATION_PATH, f'plate_{plate}_mapping.png')
            plt.savefig(fig_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"  Created visualization for Plate {plate}")
        
        # Create a summary visualization across all plates
        create_plate_summary_visualization(cross_ref_results, all_plates)
        
    except Exception as e:
        print(f"Error creating visualization: {e}")
        import traceback
        traceback.print_exc()

def create_plate_summary_visualization(cross_ref_results, all_plates):
    """Create a summary visualization across all plates
    
    Args:
        cross_ref_results: Results from cross_reference_datasets
        all_plates: List of all plate numbers
    """
    try:
        plate_breakdown = cross_ref_results['plate_breakdown']
        
        if not plate_breakdown:
            print("  No plate breakdown data available for summary visualization")
            return
        
        # Prepare data for plotting
        plates = [item['plate'] for item in plate_breakdown]
        common = [item['common_count'] for item in plate_breakdown]
        benomyl_only = [item['benomyl_only_count'] for item in plate_breakdown]
        flow_only = [item['flow_only_count'] for item in plate_breakdown]
        
        # Create stacked bar chart
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        
        # Convert plates to strings for x-axis labels
        x_pos = np.arange(len(plates))
        
        # Plot stacked bars
        ax.bar(x_pos, common, label='Both Datasets', color='#99ff99')
        ax.bar(x_pos, benomyl_only, bottom=common, label='Benomyl Only', color='#ffcccc')
        ax.bar(x_pos, flow_only, bottom=np.array(common) + np.array(benomyl_only), 
               label='Flow Only', color='#ccccff')
        
        # Set x-ticks and labels
        ax.set_xticks(x_pos)
        ax.set_xticklabels([str(p) for p in plates])
        
        ax.set_title('Well Distribution by Plate')
        ax.set_xlabel('Plate Number')
        ax.set_ylabel('Number of Wells')
        ax.legend()
        
        # Add percentages on top of each bar
        for i, plate in enumerate(plates):
            total = common[i] + benomyl_only[i] + flow_only[i]
            if total > 0:
                common_pct = common[i] / total * 100
                
                # Place text above each bar
                ax.text(i, total + 2, f"{common_pct:.1f}%\nmatched", 
                        ha='center', va='bottom', fontsize=9)
        
        # Save the figure
        plt.tight_layout()
        summary_file = os.path.join(VALIDATION_PATH, 'all_plates_summary.png')
        plt.savefig(summary_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  Created summary visualization for all plates")
        
    except Exception as e:
        print(f"Error creating plate summary visualization: {e}")
        import traceback
        traceback.print_exc()

def main():
    """Main function to execute the validation workflow"""
    print("Starting ploidy data mapping validation workflow...")
    
    # Load data
    benomyl_data, benomyl_df = load_benomyl_data()
    flow_data, flow_raw = load_flow_cytometry_data()
    
    # 1. Validate sample ID formats
    benomyl_format_df, flow_format_df, benomyl_plate_counts, flow_plate_counts = validate_sample_id_formats(benomyl_data, flow_data)
    
    # 2. Cross-reference datasets
    cross_ref_results = cross_reference_datasets(benomyl_data, flow_data)
    
    # 3. Create visual mapping
    create_well_mapping_visualization(benomyl_data, flow_data, cross_ref_results)
    
    print("Ploidy data mapping validation workflow completed successfully!")
    print(f"All validation outputs saved to: {VALIDATION_PATH}")

if __name__ == "__main__":
    main()