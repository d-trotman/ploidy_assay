#!/usr/bin/env python3
"""
Revised Ploidy Comparison Script with Ambiguous Gate

This script compares ploidy results from flow cytometry data and benomyl assay data,
and calculates the number of cells in the ambiguous gate (cells not falling into 
either haploid or diploid gates) for informational purposes.

Usage:
    python revised_ploidy_comparison.py

Output:
    CSV file with comparison results including ambiguous gate counts
"""

import os
import pandas as pd
import numpy as np
import openpyxl
import re

# Define paths
BENOMYL_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/raw_data/benamil_assay'
FLOW_CYTOMETRY_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/raw_data/flow_cytometry'
OUTPUT_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/analysis_output'
CODE_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/code'

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_PATH, exist_ok=True)

def format_well_id(plate, row, column):
    """Format well ID from plate, row, and column"""
    return f"WGS{plate}{row}{column}"

def extract_plate_row_col(well_id):
    """Extract plate, row, and column from well ID"""
    match = re.match(r'WGS(\d+)([A-Z]+)(\d+)', well_id)
    if match:
        plate = int(match.group(1))
        row = match.group(2)
        column = int(match.group(3))
        return plate, row, column
    return None, None, None

def load_benomyl_data():
    """Load benomyl assay data"""
    print("Loading benomyl assay data...")
    benomyl_file = os.path.join(BENOMYL_PATH, 'SM_benomyl_assay_df.csv')
    benomyl_df = pd.read_csv(benomyl_file)
    
    # Create a dictionary for easy lookup
    benomyl_data = {}
    for _, row in benomyl_df.iterrows():
        well_id = format_well_id(row['plate'], row['row'], row['column'])
        benomyl_data[well_id] = {
            'is_diploid': row['growth inhibited (diploid)'],
            'is_borderline': row['borderline'],
            'notes': row['notes'] if not pd.isna(row['notes']) else "",
            'double_checked': row['doublechecked'] if not pd.isna(row['doublechecked']) else False
        }
    
    print(f"  Loaded {len(benomyl_data)} benomyl data points")
    return benomyl_data

def extract_flow_data(file_path):
    """Extract flow cytometry data with ambiguous gate calculation"""
    file_name = os.path.basename(file_path)
    print(f"Processing flow cytometry data from {file_name}...")
    
    flow_data = {}
    excluded_count = 0
    
    # Load the Excel file
    wb = openpyxl.load_workbook(file_path, data_only=True)
    sheet = wb.active
    
    # Find column indices (assume header is in row 1)
    header_row = [cell.value for cell in sheet[1]]
    
    # Common columns in standardized format
    sample_name_col = header_row.index('sample name') + 1 if 'sample name' in header_row else -1
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
    
    # Extract data starting from row 2
    for row in sheet.iter_rows(min_row=2):
        # Get sample ID based on file format
        if sample_name_col > 0:
            sample_id = row[sample_name_col - 1].value
        else:
            # Fall back to alternative sample ID location if needed
            sample_id = None
            # Add any specific handling for non-standardized files here
        
        if not sample_id or not isinstance(sample_id, str) or not sample_id.startswith('WGS'):
            continue
        
        # Get event counts
        all_events_count = row[all_events_count_col - 1].value or 0
        if all_events_count < 5000:
            excluded_count += 1
            continue  # Skip samples with All Events Count < 5000
        
        singlets_count = row[singlets_count_col - 1].value or 0
        haploid_count = row[haploid_count_col - 1].value or 0
        diploid_count = row[diploid_count_col - 1].value or 0
        
        # Calculate ambiguous gate (cells that don't fall into either haploid or diploid gates)
        ambiguous_count = singlets_count - (haploid_count + diploid_count)
        
        # Calculate percentages relative to singlets
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
            # If no declared ploidy, use counts and percentages
            is_haploid = (haploid_count > diploid_count and haploid_percent > diploid_percent)
            is_diploid = (diploid_count > haploid_count and diploid_percent > haploid_percent)
        
        # Store the data with the well ID as key
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
            'source_file': file_name
        }
    
    print(f"  Extracted {len(flow_data)} samples from {file_name} (excluded {excluded_count} samples with All Events Count < 5000)")
    return flow_data

def load_flow_cytometry_data():
    """Load all flow cytometry data"""
    print("Loading flow cytometry data...")
    
    wgs1_file = os.path.join(FLOW_CYTOMETRY_PATH, 'WGS1_ploidy_assay.xlsx')
    wgs2_file = os.path.join(FLOW_CYTOMETRY_PATH, 'WGS2_ploid_assay.xlsx')
    wgs_hybrid_file = os.path.join(FLOW_CYTOMETRY_PATH, 'WGShybrid_ploidy_assay.xlsx')
    
    # Extract data from each file
    wgs1_data = extract_flow_data(wgs1_file)
    wgs2_data = extract_flow_data(wgs2_file)
    wgs_hybrid_data = extract_flow_data(wgs_hybrid_file)
    
    # Combine all data
    combined_data = {**wgs1_data, **wgs2_data, **wgs_hybrid_data}
    print(f"  Combined {len(combined_data)} flow cytometry data points")
    
    return combined_data

def compare_methods(benomyl_data, flow_data):
    """Compare benomyl and flow cytometry results"""
    print("Comparing benomyl and flow cytometry results...")
    
    comparison_results = []
    
    # Get all unique well IDs across both datasets
    all_well_ids = set(benomyl_data.keys()).union(set(flow_data.keys()))
    print(f"  Found {len(all_well_ids)} unique well IDs")
    
    # Compare data for each well
    for well_id in all_well_ids:
        benomyl_info = benomyl_data.get(well_id)
        flow_info = flow_data.get(well_id)
        
        plate, row, column = extract_plate_row_col(well_id)
        
        if benomyl_info and flow_info:
            # Both methods have data for this well
            benomyl_says_diploid = benomyl_info['is_diploid']
            flow_says_diploid = flow_info['is_diploid']
            
            # Check agreement between methods
            agreement = benomyl_says_diploid == flow_says_diploid
            
            comparison_results.append({
                'well_id': well_id,
                'plate': plate,
                'row': row,
                'column': column,
                'benomyl_says_diploid': benomyl_says_diploid,
                'flow_says_diploid': flow_says_diploid,
                'benomyl_is_borderline': benomyl_info['is_borderline'],
                'agreement': agreement,
                'benomyl_notes': benomyl_info['notes'],
                'benomyl_double_checked': benomyl_info['double_checked'],
                'flow_all_events_count': flow_info['all_events_count'],
                'flow_singlets_count': flow_info['singlets_count'],
                'flow_haploid_count': flow_info['haploid_count'],
                'flow_diploid_count': flow_info['diploid_count'],
                'flow_ambiguous_count': flow_info['ambiguous_count'],
                'flow_haploid_percent': flow_info['haploid_percent'],
                'flow_diploid_percent': flow_info['diploid_percent'],
                'flow_ambiguous_percent': flow_info['ambiguous_percent'],
                'flow_source_file': flow_info['source_file'],
                'in_benomyl': True,
                'in_flow': True
            })
        elif benomyl_info:
            # Only benomyl data available
            comparison_results.append({
                'well_id': well_id,
                'plate': plate,
                'row': row,
                'column': column,
                'benomyl_says_diploid': benomyl_info['is_diploid'],
                'flow_says_diploid': None,
                'benomyl_is_borderline': benomyl_info['is_borderline'],
                'agreement': None,
                'benomyl_notes': benomyl_info['notes'],
                'benomyl_double_checked': benomyl_info['double_checked'],
                'flow_all_events_count': None,
                'flow_singlets_count': None,
                'flow_haploid_count': None,
                'flow_diploid_count': None,
                'flow_ambiguous_count': None,
                'flow_haploid_percent': None,
                'flow_diploid_percent': None,
                'flow_ambiguous_percent': None,
                'flow_source_file': None,
                'in_benomyl': True,
                'in_flow': False
            })
        elif flow_info:
            # Only flow cytometry data available
            comparison_results.append({
                'well_id': well_id,
                'plate': plate,
                'row': row,
                'column': column,
                'benomyl_says_diploid': None,
                'flow_says_diploid': flow_info['is_diploid'],
                'benomyl_is_borderline': None,
                'agreement': None,
                'benomyl_notes': None,
                'benomyl_double_checked': None,
                'flow_all_events_count': flow_info['all_events_count'],
                'flow_singlets_count': flow_info['singlets_count'],
                'flow_haploid_count': flow_info['haploid_count'],
                'flow_diploid_count': flow_info['diploid_count'],
                'flow_ambiguous_count': flow_info['ambiguous_count'],
                'flow_haploid_percent': flow_info['haploid_percent'],
                'flow_diploid_percent': flow_info['diploid_percent'],
                'flow_ambiguous_percent': flow_info['ambiguous_percent'],
                'flow_source_file': flow_info['source_file'],
                'in_benomyl': False,
                'in_flow': True
            })
    
    print(f"  Completed {len(comparison_results)} comparisons")
    return comparison_results

def generate_summary_statistics(comparison_results):
    """Generate summary statistics"""
    print("Generating summary statistics...")
    
    # Filter for wells with data from both methods
    both_methods = [r for r in comparison_results if r['in_benomyl'] and r['in_flow']]
    
    # Count agreements and disagreements
    agreements = [r for r in both_methods if r['agreement']]
    disagreements = [r for r in both_methods if r['agreement'] is False]
    
    # Count borderline cases
    borderline_cases = [r for r in both_methods if r['benomyl_is_borderline']]
    
    # Count by plate
    by_plate = {}
    for r in both_methods:
        plate = r['plate']
        if plate not in by_plate:
            by_plate[plate] = {'total': 0, 'agreements': 0, 'disagreements': 0, 'borderline': 0}
        
        by_plate[plate]['total'] += 1
        if r['agreement']:
            by_plate[plate]['agreements'] += 1
        elif r['agreement'] is False:
            by_plate[plate]['disagreements'] += 1
        
        if r['benomyl_is_borderline']:
            by_plate[plate]['borderline'] += 1
    
    # Calculate percentages
    total_both = len(both_methods)
    percent_agreement = len(agreements) / total_both * 100 if total_both > 0 else 0
    percent_disagreement = len(disagreements) / total_both * 100 if total_both > 0 else 0
    percent_borderline = len(borderline_cases) / total_both * 100 if total_both > 0 else 0
    
    # Basic statistics on ambiguous gate percentages
    flow_results = [r for r in comparison_results if r['in_flow']]
    ambiguous_percentages = [r['flow_ambiguous_percent'] for r in flow_results]
    
    summary = {
        'total_wells': len(comparison_results),
        'wells_with_both_methods': total_both,
        'wells_with_only_benomyl': len([r for r in comparison_results if r['in_benomyl'] and not r['in_flow']]),
        'wells_with_only_flow': len([r for r in comparison_results if r['in_flow'] and not r['in_benomyl']]),
        'agreements': len(agreements),
        'disagreements': len(disagreements),
        'borderline_cases': len(borderline_cases),
        'percent_agreement': percent_agreement,
        'percent_disagreement': percent_disagreement,
        'percent_borderline': percent_borderline,
        'by_plate': by_plate,
        'ambiguous_gate_stats': {
            'count': len(ambiguous_percentages),
            'mean': np.mean(ambiguous_percentages) if ambiguous_percentages else 0,
            'median': np.median(ambiguous_percentages) if ambiguous_percentages else 0,
            'min': np.min(ambiguous_percentages) if ambiguous_percentages else 0,
            'max': np.max(ambiguous_percentages) if ambiguous_percentages else 0,
            'std': np.std(ambiguous_percentages) if ambiguous_percentages else 0
        }
    }
    
    print("  Summary statistics generated")
    return summary

def main():
    """Main function to execute the comparison workflow"""
    print("Starting ploidy comparison workflow with ambiguous gate calculation...")
    
    # Load data
    benomyl_data = load_benomyl_data()
    flow_data = load_flow_cytometry_data()
    
    # Compare methods
    comparison_results = compare_methods(benomyl_data, flow_data)
    
    # Generate summary statistics
    summary = generate_summary_statistics(comparison_results)
    
    # Convert comparison results to DataFrame
    df = pd.DataFrame(comparison_results)
    
    # Sort by plate, row, column
    df.sort_values(['plate', 'row', 'column'], inplace=True)
    
    # Save comparison results to CSV
    output_file = os.path.join(OUTPUT_PATH, 'ploidy_comparison_with_ambiguous_gate.csv')
    df.to_csv(output_file, index=False)
    print(f"Comparison results saved to {output_file}")
    
    # Create a summary file
    summary_file = os.path.join(OUTPUT_PATH, 'ploidy_comparison_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("Ploidy Comparison Summary\n")
        f.write("========================\n\n")
        f.write(f"Total wells analyzed: {summary['total_wells']}\n")
        f.write(f"Wells with both methods: {summary['wells_with_both_methods']}\n")
        f.write(f"Wells with only benomyl: {summary['wells_with_only_benomyl']}\n")
        f.write(f"Wells with only flow cytometry: {summary['wells_with_only_flow']}\n\n")
        
        f.write("Method Agreement:\n")
        f.write(f"  Agreements: {summary['agreements']} ({summary['percent_agreement']:.2f}%)\n")
        f.write(f"  Disagreements: {summary['disagreements']} ({summary['percent_disagreement']:.2f}%)\n")
        f.write(f"  Borderline cases: {summary['borderline_cases']} ({summary['percent_borderline']:.2f}%)\n\n")
        
        f.write("Ambiguous Gate Statistics:\n")
        f.write(f"  Mean percentage: {summary['ambiguous_gate_stats']['mean']:.2f}%\n")
        f.write(f"  Median percentage: {summary['ambiguous_gate_stats']['median']:.2f}%\n")
        f.write(f"  Min percentage: {summary['ambiguous_gate_stats']['min']:.2f}%\n")
        f.write(f"  Max percentage: {summary['ambiguous_gate_stats']['max']:.2f}%\n")
        f.write(f"  Standard deviation: {summary['ambiguous_gate_stats']['std']:.2f}%\n\n")
        
        f.write("Results by Plate:\n")
        for plate, stats in summary['by_plate'].items():
            agreement_pct = stats['agreements'] / stats['total'] * 100 if stats['total'] > 0 else 0
            f.write(f"  Plate {plate}:\n")
            f.write(f"    Total comparisons: {stats['total']}\n")
            f.write(f"    Agreements: {stats['agreements']} ({agreement_pct:.2f}%)\n")
            f.write(f"    Disagreements: {stats['disagreements']} ({(100-agreement_pct):.2f}%)\n")
            f.write(f"    Borderline cases: {stats['borderline']}\n\n")
    
    print(f"Summary saved to {summary_file}")
    
    # Create a flagged results file for disagreements
    disagreement_df = df[(df['agreement'] == False) & (df['in_benomyl']) & (df['in_flow'])]
    flagged_file = os.path.join(OUTPUT_PATH, 'ploidy_comparison_disagreements.csv')
    disagreement_df.to_csv(flagged_file, index=False)
    print(f"Flagged disagreements saved to {flagged_file}")
    
    # Create a file with flow cytometry results including ambiguous gate
    flow_df = df[df['in_flow']].copy()
    flow_columns = ['well_id', 'plate', 'row', 'column', 
                   'flow_says_diploid', 'flow_haploid_count', 'flow_diploid_count', 
                   'flow_ambiguous_count', 'flow_singlets_count', 'flow_all_events_count',
                   'flow_haploid_percent', 'flow_diploid_percent', 'flow_ambiguous_percent']
    
    flow_output = os.path.join(OUTPUT_PATH, 'flow_cytometry_with_ambiguous_gate.csv')
    flow_df[flow_columns].to_csv(flow_output, index=False)
    print(f"Flow cytometry results with ambiguous gate saved to {flow_output}")
    
    print("Ploidy comparison workflow completed successfully!")

if __name__ == "__main__":
    main()