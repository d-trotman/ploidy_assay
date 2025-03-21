#!/usr/bin/env python3
"""
Ploidy Method Discrepancy Analysis

This script analyzes discrepancies between flow cytometry and benomyl assay results
for ploidy determination in S. cerevisiae. It calculates the ambiguous gate in flow 
cytometry data and conducts targeted analyses on disagreement cases.

Usage:
    python ploidy_discrepancy_analysis.py

Output:
    CSV files with comparison results and analysis of discrepancies
"""

import os
import pandas as pd
import numpy as np
import openpyxl
import re
import matplotlib.pyplot as plt
from scipy import stats

# Define paths
BENOMYL_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/raw_data/benamil_assay'
FLOW_CYTOMETRY_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/raw_data/flow_cytometry'
OUTPUT_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/analysis_output'
FIGURES_PATH = os.path.join(OUTPUT_PATH, 'figures')
CODE_PATH = '/Users/dawsontrotman/Documents/GitHub/ploidy_assay/code'

# Create output directories if they don't exist
os.makedirs(OUTPUT_PATH, exist_ok=True)
os.makedirs(FIGURES_PATH, exist_ok=True)

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
        
        # Determine if this is a "pure" population or mixed/ambiguous
        # High percentages in one gate suggest pure populations
        dominant_gate_percent = max(haploid_percent, diploid_percent)
        mixed_population = (dominant_gate_percent < 70) or (ambiguous_percent > 30)
        
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
            'mixed_population': mixed_population,
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
    """Compare benomyl and flow cytometry results with discrepancy analysis"""
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
            
            # Classify the type of disagreement if there is one
            disagreement_type = None
            if not agreement:
                if benomyl_says_diploid and not flow_says_diploid:
                    disagreement_type = "benomyl_diploid_flow_haploid"
                elif not benomyl_says_diploid and flow_says_diploid:
                    disagreement_type = "benomyl_haploid_flow_diploid"
            
            # Calculate confidence scores for each method
            # Benomyl confidence: lower if borderline
            benomyl_confidence = 0.5 if benomyl_info['is_borderline'] else 1.0
            
            # Flow confidence: based on how clear the dominant population is
            dominant_percent = max(flow_info['haploid_percent'], flow_info['diploid_percent'])
            flow_confidence = dominant_percent / 100
            
            # Flag potential issues
            mixed_population = flow_info['mixed_population']
            borderline_benomyl = benomyl_info['is_borderline']
            
            # Potential explanation for disagreement
            disagreement_explanation = None
            if not agreement:
                if mixed_population and borderline_benomyl:
                    disagreement_explanation = "Both methods show uncertainty"
                elif mixed_population:
                    disagreement_explanation = "Mixed or unclear population in flow cytometry"
                elif borderline_benomyl:
                    disagreement_explanation = "Borderline benomyl result"
                else:
                    disagreement_explanation = "Clear disagreement with high confidence both methods"
            
            comparison_results.append({
                'well_id': well_id,
                'plate': plate,
                'row': row,
                'column': column,
                'benomyl_says_diploid': benomyl_says_diploid,
                'flow_says_diploid': flow_says_diploid,
                'benomyl_is_borderline': benomyl_info['is_borderline'],
                'agreement': agreement,
                'disagreement_type': disagreement_type,
                'disagreement_explanation': disagreement_explanation,
                'benomyl_confidence': benomyl_confidence,
                'flow_confidence': flow_confidence,
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
                'flow_mixed_population': flow_info['mixed_population'],
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
                'disagreement_type': None,
                'disagreement_explanation': None,
                'benomyl_confidence': 0.5 if benomyl_info['is_borderline'] else 1.0,
                'flow_confidence': None,
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
                'flow_mixed_population': None,
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
                'disagreement_type': None,
                'disagreement_explanation': None,
                'benomyl_confidence': None,
                'flow_confidence': max(flow_info['haploid_percent'], flow_info['diploid_percent']) / 100,
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
                'flow_mixed_population': flow_info['mixed_population'],
                'flow_source_file': flow_info['source_file'],
                'in_benomyl': False,
                'in_flow': True
            })
    
    print(f"  Completed {len(comparison_results)} comparisons")
    return comparison_results

def analyze_discrepancies(comparison_results):
    """Analyze the patterns and potential causes of discrepancies between methods"""
    print("Analyzing discrepancies between methods...")
    
    # Convert to DataFrame for easier analysis
    df = pd.DataFrame(comparison_results)
    
    # Only consider wells with data from both methods
    both_methods_df = df[(df['in_benomyl']) & (df['in_flow'])]
    
    if len(both_methods_df) == 0:
        print("  No wells with data from both methods found")
        return None
    
    # Split into agreement and disagreement groups
    agreement_df = both_methods_df[both_methods_df['agreement'] == True]
    disagreement_df = both_methods_df[both_methods_df['agreement'] == False]
    
    # Further split disagreements by type
    benomyl_diploid_flow_haploid = disagreement_df[disagreement_df['disagreement_type'] == 'benomyl_diploid_flow_haploid']
    benomyl_haploid_flow_diploid = disagreement_df[disagreement_df['disagreement_type'] == 'benomyl_haploid_flow_diploid']
    
    # Calculate how many disagreements involve borderline benomyl results
    borderline_disagreements = disagreement_df[disagreement_df['benomyl_is_borderline'] == True]
    
    # Calculate how many disagreements involve mixed flow cytometry populations
    mixed_flow_disagreements = disagreement_df[disagreement_df['flow_mixed_population'] == True]
    
    # Calculate how many disagreements are "clean" (high confidence in both methods)
    clean_disagreements = disagreement_df[
        (disagreement_df['benomyl_is_borderline'] == False) & 
        (disagreement_df['flow_mixed_population'] == False)
    ]
    
    # Calculate statistics
    discrepancy_analysis = {
        'total_comparisons': len(both_methods_df),
        'agreements': len(agreement_df),
        'disagreements': len(disagreement_df),
        'benomyl_diploid_flow_haploid': len(benomyl_diploid_flow_haploid),
        'benomyl_haploid_flow_diploid': len(benomyl_haploid_flow_diploid),
        'borderline_benomyl_disagreements': len(borderline_disagreements),
        'mixed_flow_disagreements': len(mixed_flow_disagreements),
        'clean_disagreements': len(clean_disagreements),
        'percent_agreement': len(agreement_df) / len(both_methods_df) * 100 if len(both_methods_df) > 0 else 0,
        'percent_borderline_in_disagreements': len(borderline_disagreements) / len(disagreement_df) * 100 if len(disagreement_df) > 0 else 0,
        'percent_mixed_flow_in_disagreements': len(mixed_flow_disagreements) / len(disagreement_df) * 100 if len(disagreement_df) > 0 else 0,
        'percent_clean_disagreements': len(clean_disagreements) / len(disagreement_df) * 100 if len(disagreement_df) > 0 else 0
    }
    
    # For the clean disagreements, identify if there's a pattern by plate or row/column
    if len(clean_disagreements) > 0:
        clean_by_plate = clean_disagreements.groupby('plate').size().to_dict()
        clean_by_row = clean_disagreements.groupby('row').size().to_dict()
        clean_by_column = clean_disagreements.groupby('column').size().to_dict()
        
        discrepancy_analysis['clean_by_plate'] = clean_by_plate
        discrepancy_analysis['clean_by_row'] = clean_by_row
        discrepancy_analysis['clean_by_column'] = clean_by_column
    
    # Create some visualizations to help understand the patterns
    create_discrepancy_visualizations(df, both_methods_df, agreement_df, disagreement_df)
    
    print("  Completed discrepancy analysis")
    return discrepancy_analysis

def create_discrepancy_visualizations(df, both_methods_df, agreement_df, disagreement_df):
    """Create visualizations to help understand the patterns of disagreement"""
    print("Creating visualizations for discrepancy analysis...")
    
    # 1. Bar chart showing agreements vs. disagreements by plate
    plate_agreements = both_methods_df.groupby(['plate', 'agreement']).size().unstack().fillna(0)
    
    if True in plate_agreements.columns and False in plate_agreements.columns:
        plt.figure(figsize=(12, 6))
        plate_agreements.plot(kind='bar', figsize=(12, 6))
        plt.title('Method Agreement by Plate')
        plt.xlabel('Plate')
        plt.ylabel('Count')
        plt.legend(['Disagree', 'Agree'])
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.savefig(os.path.join(FIGURES_PATH, 'agreement_by_plate.png'))
        plt.close()
    
    # 2. Grouped bar chart of ambiguous percentages in agreement vs disagreement cases
    plt.figure(figsize=(10, 6))
    datasets = [
        agreement_df['flow_ambiguous_percent'], 
        disagreement_df['flow_ambiguous_percent']
    ]
    
    plt.boxplot(datasets, labels=['Methods Agree', 'Methods Disagree'])
    plt.ylabel('Ambiguous Gate Percentage (%)')
    plt.title('Ambiguous Gate Percentage by Method Agreement')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(FIGURES_PATH, 'ambiguous_by_agreement.png'))
    plt.close()
    
    # 3. Pie chart showing breakdown of disagreement types
    if len(disagreement_df) > 0:
        disagreement_types = disagreement_df['disagreement_type'].value_counts()
        
        plt.figure(figsize=(10, 7))
        plt.pie(disagreement_types, labels=disagreement_types.index, autopct='%1.1f%%')
        plt.title('Types of Method Disagreement')
        plt.savefig(os.path.join(FIGURES_PATH, 'disagreement_types.png'))
        plt.close()
    
    # 4. Pie chart showing breakdown of disagreement explanations
    if len(disagreement_df) > 0:
        explanation_counts = disagreement_df['disagreement_explanation'].value_counts()
        
        plt.figure(figsize=(10, 7))
        plt.pie(explanation_counts, labels=explanation_counts.index, autopct='%1.1f%%')
        plt.title('Explanations for Method Disagreements')
        plt.savefig(os.path.join(FIGURES_PATH, 'disagreement_explanations.png'))
        plt.close()
    
    # 5. Scatter plot of haploid vs diploid percentages, colored by agreement
    plt.figure(figsize=(10, 8))
    plt.scatter(
        agreement_df['flow_haploid_percent'], 
        agreement_df['flow_diploid_percent'], 
        alpha=0.7, label='Methods Agree', color='blue'
    )
    plt.scatter(
        disagreement_df['flow_haploid_percent'], 
        disagreement_df['flow_diploid_percent'], 
        alpha=0.7, label='Methods Disagree', color='red'
    )
    plt.xlabel('Haploid Percentage (%)')
    plt.ylabel('Diploid Percentage (%)')
    plt.title('Flow Cytometry Cell Distribution by Method Agreement')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig(os.path.join(FIGURES_PATH, 'haploid_vs_diploid_by_agreement.png'))
    plt.close()
    
    print(f"  Visualizations saved to {FIGURES_PATH}")

def generate_summary_statistics(comparison_results, discrepancy_analysis):
    """Generate summary statistics including discrepancy analysis"""
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
        },
        'discrepancy_analysis': discrepancy_analysis
    }
    
    print("  Summary statistics generated")
    return summary

def main():
    """Main function to execute the discrepancy analysis workflow"""
    print("Starting ploidy method discrepancy analysis...")
    
    # Load data
    benomyl_data = load_benomyl_data()
    flow_data = load_flow_cytometry_data()
    
    # Compare methods
    comparison_results = compare_methods(benomyl_data, flow_data)
    
    # Analyze discrepancies
    discrepancy_analysis = analyze_discrepancies(comparison_results)
    
    # Generate summary statistics
    summary = generate_summary_statistics(comparison_results, discrepancy_analysis)
    
    # Convert comparison results to DataFrame
    df = pd.DataFrame(comparison_results)
    
    # Sort by plate, row, column
    df.sort_values(['plate', 'row', 'column'], inplace=True)
    
    # Save comparison results to CSV
    output_file = os.path.join(OUTPUT_PATH, 'ploidy_method_comparison.csv')
    df.to_csv(output_file, index=False)
    print(f"Comparison results saved to {output_file}")
    
    # Create a summary file
    summary_file = os.path.join(OUTPUT_PATH, 'ploidy_method_discrepancy_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("Ploidy Method Discrepancy Analysis\n")
        f.write("=================================\n\n")
        
        f.write("OVERVIEW\n")
        f.write("--------\n")
        f.write(f"Total wells analyzed: {summary['total_wells']}\n")
        f.write(f"Wells with both methods: {summary['wells_with_both_methods']}\n")
        f.write(f"Wells with only benomyl: {summary['wells_with_only_benomyl']}\n")
        f.write(f"Wells with only flow cytometry: {summary['wells_with_only_flow']}\n\n")
        
        f.write("METHOD AGREEMENT\n")
        f.write("--------------\n")
        f.write(f"Agreements: {summary['agreements']} ({summary['percent_agreement']:.2f}%)\n")
        f.write(f"Disagreements: {summary['disagreements']} ({summary['percent_disagreement']:.2f}%)\n")
        f.write(f"Borderline benomyl cases: {summary['borderline_cases']} ({summary['percent_borderline']:.2f}%)\n\n")
        
        if discrepancy_analysis:
            f.write("DISCREPANCY ANALYSIS\n")
            f.write("-------------------\n")
            f.write("Types of disagreements:\n")
            f.write(f"  Benomyl says diploid, Flow says haploid: {discrepancy_analysis['benomyl_diploid_flow_haploid']} cases\n")
            f.write(f"  Benomyl says haploid, Flow says diploid: {discrepancy_analysis['benomyl_haploid_flow_diploid']} cases\n\n")
            
            f.write("Potential explanations for disagreements:\n")
            f.write(f"  Borderline benomyl results: {discrepancy_analysis['borderline_benomyl_disagreements']} cases ")
            f.write(f"({discrepancy_analysis['percent_borderline_in_disagreements']:.2f}% of disagreements)\n")
            f.write(f"  Mixed flow cytometry populations: {discrepancy_analysis['mixed_flow_disagreements']} cases ")
            f.write(f"({discrepancy_analysis['percent_mixed_flow_in_disagreements']:.2f}% of disagreements)\n")
            f.write(f"  Clean disagreements: {discrepancy_analysis['clean_disagreements']} cases ")
            f.write(f"({discrepancy_analysis['percent_clean_disagreements']:.2f}% of disagreements)\n\n")
            
            if discrepancy_analysis.get('clean_by_plate'):
                f.write("Distribution of clean disagreements:\n")
                f.write("  By plate:\n")
                for plate, count in discrepancy_analysis['clean_by_plate'].items():
                    f.write(f"    Plate {plate}: {count} cases\n")
                
                f.write("  By row:\n")
                for row, count in discrepancy_analysis['clean_by_row'].items():
                    f.write(f"    Row {row}: {count} cases\n")
                
                f.write("  By column:\n")
                for col, count in discrepancy_analysis['clean_by_column'].items():
                    f.write(f"    Column {col}: {count} cases\n\n")
        
        f.write("AMBIGUOUS GATE STATISTICS\n")
        f.write("------------------------\n")
        f.write(f"Mean ambiguous gate percentage: {summary['ambiguous_gate_stats']['mean']:.2f}%\n")
        f.write(f"Median ambiguous gate percentage: {summary['ambiguous_gate_stats']['median']:.2f}%\n")
        f.write(f"Min ambiguous gate percentage: {summary['ambiguous_gate_stats']['min']:.2f}%\n")
        f.write(f"Max ambiguous gate percentage: {summary['ambiguous_gate_stats']['max']:.2f}%\n")
        f.write(f"Standard deviation: {summary['ambiguous_gate_stats']['std']:.2f}%\n\n")
        
        f.write("RESULTS BY PLATE\n")
        f.write("---------------\n")
        for plate, stats in summary['by_plate'].items():
            agreement_pct = stats['agreements'] / stats['total'] * 100 if stats['total'] > 0 else 0
            f.write(f"Plate {plate}:\n")
            f.write(f"  Total comparisons: {stats['total']}\n")
            f.write(f"  Agreements: {stats['agreements']} ({agreement_pct:.2f}%)\n")
            f.write(f"  Disagreements: {stats['disagreements']} ({(100-agreement_pct):.2f}%)\n")
            f.write(f"  Borderline cases: {stats['borderline']}\n\n")
        
        f.write("CONCLUSIONS\n")
        f.write("-----------\n")
        f.write("1. Overall agreement between methods: ")
        if summary['percent_agreement'] > 75:
            f.write("Strong agreement (>75%)\n")
        elif summary['percent_agreement'] > 50:
            f.write("Moderate agreement (50-75%)\n")
        else:
            f.write("Poor agreement (<50%)\n")
        
        if discrepancy_analysis:
            f.write("2. Primary source of disagreement: ")
            max_pct = max(
                discrepancy_analysis['percent_borderline_in_disagreements'],
                discrepancy_analysis['percent_mixed_flow_in_disagreements'],
                discrepancy_analysis['percent_clean_disagreements']
            )
            
            if max_pct == discrepancy_analysis['percent_borderline_in_disagreements']:
                f.write("Borderline benomyl results\n")
            elif max_pct == discrepancy_analysis['percent_mixed_flow_in_disagreements']:
                f.write("Mixed populations in flow cytometry\n")
            else:
                f.write("Clean disagreements (likely biological differences in method sensitivity)\n")
            
            f.write("3. Predominant disagreement type: ")
            if discrepancy_analysis['benomyl_diploid_flow_haploid'] > discrepancy_analysis['benomyl_haploid_flow_diploid']:
                f.write("Benomyl tends to classify as diploid when flow cytometry indicates haploid\n")
            elif discrepancy_analysis['benomyl_diploid_flow_haploid'] < discrepancy_analysis['benomyl_haploid_flow_diploid']:
                f.write("Benomyl tends to classify as haploid when flow cytometry indicates diploid\n")
            else:
                f.write("No clear pattern in disagreement types\n")
    
    print(f"Analysis summary saved to {summary_file}")
    
    # Create specialized output files for different types of disagreements
    # File 1: Cases with borderline benomyl results
    borderline_df = df[(df['agreement'] == False) & (df['benomyl_is_borderline'] == True)]
    borderline_file = os.path.join(OUTPUT_PATH, 'discrepancy_borderline_benomyl.csv')
    borderline_df.to_csv(borderline_file, index=False)
    print(f"Borderline benomyl disagreements saved to {borderline_file}")
    
    # File 2: Cases with mixed flow cytometry populations
    mixed_flow_df = df[(df['agreement'] == False) & (df['flow_mixed_population'] == True)]
    mixed_flow_file = os.path.join(OUTPUT_PATH, 'discrepancy_mixed_flow.csv')
    mixed_flow_df.to_csv(mixed_flow_file, index=False)
    print(f"Mixed flow population disagreements saved to {mixed_flow_file}")
    
    # File 3: Clean disagreements (high confidence from both methods)
    clean_df = df[(df['agreement'] == False) & 
                 (df['benomyl_is_borderline'] == False) & 
                 (df['flow_mixed_population'] == False)]
    clean_file = os.path.join(OUTPUT_PATH, 'discrepancy_clean.csv')
    clean_df.to_csv(clean_file, index=False)
    print(f"Clean disagreements saved to {clean_file}")
    
    print("Ploidy method discrepancy analysis completed successfully!")

if __name__ == "__main__":
    main()