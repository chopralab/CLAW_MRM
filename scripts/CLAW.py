"""
Lipidomics Data Analysis Pipeline - MRM Only
Processes mzML files for lipid identification and quantification using MRM transitions
"""

import os
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import pymzml
from matplotlib import pyplot as plt

# Suppress warnings
warnings.filterwarnings('ignore')


# ============================================================================
# GLOBAL DATAFRAME - Used across parsing functions
# ============================================================================
master_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Sample_ID'])


# ============================================================================
# PROJECT SETUP AND CONFIGURATION
# ============================================================================

def setup_project_directories(database_path, project_path, project_name, 
                              data_folder, results_folder, output_filename, 
                              tolerance, remove_standards, save_data):
    """
    Setup and verify project directory structure.
    
    Args:
        database_path: Path to MRM database file
        project_path: Root project directory
        project_name: Name of the project
        data_folder: Folder containing mzML data files
        results_folder: Folder for output results
        output_filename: Base name for output files
        tolerance: Mass tolerance for ion matching (Da)
        remove_standards: Whether to filter to specific lipid classes
        save_data: Whether to save results to CSV
    
    Returns:
        tuple: Configuration parameters for downstream processing
    """
    # Create directories if they don't exist
    os.makedirs(os.path.dirname(database_path), exist_ok=True)
    os.makedirs(project_path, exist_ok=True)
    os.makedirs(data_folder, exist_ok=True)
    os.makedirs(results_folder, exist_ok=True)

    # Log configuration
    config = {
        "database_path": database_path,
        "project_path": project_path,
        "project_name": project_name,
        "data_folder": data_folder,
        "results_folder": results_folder,
        "output_filename": output_filename,
        "tolerance": tolerance,
        "remove_standards": remove_standards,
        "save_data": save_data
    }
    
    for key, value in config.items():
        print(f"{key}: {value}")
    
    return database_path, data_folder, results_folder, output_filename, tolerance, remove_standards, save_data


# ============================================================================
# MRM DATABASE HANDLING
# ============================================================================

def load_mrm_database(filepath, remove_standards=True, deuterated=False):
    """
    Load and process Multiple Reaction Monitoring (MRM) lipid database.
    
    Args:
        filepath: Path to Excel file containing MRM transitions
        remove_standards: If True, keep only common lipid classes
        deuterated: If True, adjust m/z values for deuterated lipids
    
    Returns:
        pd.DataFrame: Processed MRM database with columns:
            - Lipid: Compound name
            - Parent_Ion: Precursor m/z
            - Product_Ion: Fragment m/z
            - Class: Lipid class
            - Transition: Formatted transition string
    """
    # Read all sheets and concatenate
    raw_data = pd.read_excel(filepath, sheet_name=None)
    combined_data = pd.concat(raw_data.values(), ignore_index=True)

    # Extract and rename columns
    mrm_data = combined_data[['Compound Name', 'Parent Ion', 'Product Ion', 'Class']].copy()
    mrm_data.columns = ['Lipid', 'Parent_Ion', 'Product_Ion', 'Class']
    
    # Round m/z values to 1 decimal place
    mrm_data['Parent_Ion'] = mrm_data['Parent_Ion'].round(1)
    mrm_data['Product_Ion'] = mrm_data['Product_Ion'].round(1)
    
    # Adjust for deuterated lipids if needed
    if deuterated:
        mrm_data['Parent_Ion'] += 1
        mrm_data['Product_Ion'] += 1
    
    # Create transition string
    mrm_data['Transition'] = (mrm_data['Parent_Ion'].astype(str) + ' -> ' + 
                              mrm_data['Product_Ion'].astype(str))

    # Filter to major lipid classes if requested
    if remove_standards:
        major_classes = ['PS', 'PG', 'CE', 'PC', 'DAG', 'PE', 'TAG', 'FA', 
                        'Cer', 'CAR', 'PI', 'SM']
        mrm_data = mrm_data[mrm_data['Class'].isin(major_classes)]
    
    return mrm_data


def build_ion_lookup_dict(mrm_database):
    """
    Create a lookup dictionary for fast ion matching.
    
    Args:
        mrm_database: DataFrame with Parent_Ion, Product_Ion, Lipid, Class columns
    
    Returns:
        dict: Keys are (parent_mz, product_mz) tuples, 
              values are lists of (lipid_name, lipid_class) tuples
    """
    ion_dict = defaultdict(list)
    for _, row in mrm_database.iterrows():
        key = (row['Parent_Ion'], row['Product_Ion'])
        value = (row['Lipid'], row['Class'])
        ion_dict[key].append(value)
    return ion_dict


# ============================================================================
# MZML FILE PARSING
# ============================================================================

def parse_mzml_file(filepath, plot_chromatogram=False, target_transition=None):
    """
    Parse a single mzML file and extract MRM transition data.
    
    Args:
        filepath: Path to mzML file
        plot_chromatogram: Whether to plot extracted ion chromatogram
        target_transition: Optional tuple (parent_mz, product_mz) to plot
    
    Side Effects:
        Updates global master_df
    """
    global master_df
    
    transition_data = []
    
    run = pymzml.run.Reader(filepath, skip_chromatogram=False)
    sample_id = os.path.basename(filepath).replace('.mzML', '')
    
    for spectrum in run:
        # Extract Q1 (parent) and Q3 (product) m/z values from spectrum ID
        parent_mz = product_mz = 0
        
        for element in spectrum.ID.split(' '):
            if 'Q1=' in element:
                parent_mz = round(float(element.split('=')[1]), 1)
            if 'Q3=' in element:
                product_mz = round(float(element.split('=')[1]), 1)
        
        if parent_mz == 0 or product_mz == 0:
            continue
        
        # Plot chromatogram if requested for specific transition
        if plot_chromatogram and target_transition:
            target_parent, target_product = target_transition
            if (abs(parent_mz - target_parent) <= 0.1 and 
                abs(product_mz - target_product) <= 0.1):
                times, intensities = zip(*spectrum.peaks())
                plt.plot(times, intensities)
                plt.xlabel('Retention Time (min)')
                plt.ylabel('Intensity')
                plt.title(f'Chromatogram: {parent_mz} -> {product_mz}')
                plt.show()
        
        # Calculate total intensity for transition
        intensities = np.array([intensity for _, intensity in spectrum.peaks()])
        total_intensity = np.sum(intensities)
        
        transition_str = f"{parent_mz} -> {product_mz}"
        
        # Store summarized transition data
        transition_data.append({
            'Parent_Ion': parent_mz,
            'Product_Ion': product_mz,
            'Intensity': total_intensity,
            'Transition': transition_str,
            'Sample_ID': sample_id
        })
    
    # Append to global dataframe
    master_df = pd.concat([master_df, pd.DataFrame(transition_data)], ignore_index=True)
    
    print(f'Parsed: {filepath}')


def parse_mzml_batch(folder_path, plot_chromatogram=False, target_transition=None):
    """
    Parse all mzML files in a directory.
    
    Args:
        folder_path: Directory containing mzML files
        plot_chromatogram: Whether to plot chromatograms
        target_transition: Optional tuple (parent_mz, product_mz) to plot
    """
    files = sorted([f for f in os.listdir(folder_path) if f.endswith('.mzML')])
    
    for filename in files:
        filepath = os.path.join(folder_path, filename)
        parse_mzml_file(filepath, plot_chromatogram, target_transition)
    
    print(f'\nFinished parsing {len(files)} mzML files')


# ============================================================================
# LIPID IDENTIFICATION
# ============================================================================

def is_within_tolerance(value1, value2, tolerance=0.3):
    """Check if two m/z values are within specified tolerance."""
    return abs(value1 - value2) <= tolerance


def match_row_to_database(row, ion_lookup, tolerance=0.3):
    """
    Match a single row's ions to the database.
    
    Args:
        row: DataFrame row with Parent_Ion and Product_Ion
        ion_lookup: Dictionary from build_ion_lookup_dict()
        tolerance: Mass tolerance in Da
    
    Returns:
        pd.Series: Original row with added Lipid and Class columns
    """
    parent_mz = row['Parent_Ion']
    product_mz = row['Product_Ion']
    
    matched_lipids = []
    matched_classes = []
    
    # Check all database entries for matches within tolerance
    for (db_parent, db_product), annotations in ion_lookup.items():
        if (is_within_tolerance(parent_mz, db_parent, tolerance) and 
            is_within_tolerance(product_mz, db_product, tolerance)):
            for lipid_name, lipid_class in annotations:
                matched_lipids.append(lipid_name)
                matched_classes.append(lipid_class)
    
    # Add matched data to row
    if matched_lipids:
        row['Lipid'] = ' | '.join(matched_lipids)
        row['Class'] = ' | '.join(matched_classes)
    
    return row


def match_lipids_to_database(mrm_database, transition_df, tolerance=0.3):
    """
    Match all transitions in dataframe to lipid database.
    
    Args:
        mrm_database: MRM database from load_mrm_database()
        transition_df: DataFrame with transition data
        tolerance: Mass tolerance for matching (Da)
    
    Returns:
        pd.DataFrame: Input dataframe with added Lipid and Class columns
    """
    ion_lookup = build_ion_lookup_dict(mrm_database)
    matched_df = transition_df.apply(
        lambda row: match_row_to_database(row, ion_lookup, tolerance), 
        axis=1
    )
    return matched_df


# ============================================================================
# DATA EXPORT
# ============================================================================

def save_results_to_csv(dataframe, results_folder, filename, max_attempts=5):
    """
    Save DataFrame to CSV with automatic versioning if file exists.
    
    Args:
        dataframe: DataFrame to save
        results_folder: Output directory (absolute or relative path)
        filename: Base filename (without .csv extension)
        max_attempts: Maximum number of filename variations to try
    """
    output_dir = results_folder
    os.makedirs(output_dir, exist_ok=True)
    
    for attempt in range(max_attempts):
        if attempt == 0:
            filepath = os.path.join(output_dir, f'{filename}.csv')
        else:
            filepath = os.path.join(output_dir, f'{filename}_{attempt}.csv')
        
        if not os.path.exists(filepath):
            dataframe.to_csv(filepath, index=False)
            print(f"Results saved to: {filepath}")
            return filepath
    
    print(f"Warning: Could not save file after {max_attempts} attempts")
    return None


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def run_full_analysis(database_path, data_folder, results_folder, output_filename,
                     tolerance=0.3, remove_standards=True, save_data=False, 
                     batch_processing=True, plot_chromatogram=False, 
                     target_transition=None):
    """
    Execute complete lipidomics MRM analysis pipeline.
    
    Args:
        database_path: Path to MRM database Excel file
        data_folder: Directory with mzML files (or single file path if batch_processing=False)
        results_folder: Directory for output files
        output_filename: Base name for output CSV
        tolerance: Mass tolerance for lipid matching (Da)
        remove_standards: Whether to filter to major lipid classes
        save_data: Whether to save results to CSV
        batch_processing: If True, process all files in folder; if False, process single file
        plot_chromatogram: Whether to plot extracted ion chromatograms
        target_transition: Optional (parent_mz, product_mz) tuple for plotting
    
    Returns:
        pd.DataFrame: Dataframe containing identified lipids with intensities
    """
    global master_df
    
    # Reset global dataframe
    master_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Sample_ID'])
    
    # Load database
    print("Loading MRM database...")
    mrm_database = load_mrm_database(database_path, remove_standards=remove_standards)
    print(f"Loaded {len(mrm_database)} transitions\n")
    
    # Parse mzML files
    print("Parsing mzML files...")
    if batch_processing:
        parse_mzml_batch(data_folder, plot_chromatogram, target_transition)
    else:
        parse_mzml_file(data_folder, plot_chromatogram, target_transition)
    
    # Match lipids
    print("\nMatching transitions to lipid database...")
    matched_df = match_lipids_to_database(mrm_database, master_df, tolerance)
    print(f"Identified {matched_df['Lipid'].notna().sum()} transitions\n")
    
    # Save results
    if save_data:
        save_results_to_csv(matched_df, results_folder, output_filename)
    
    return matched_df


if __name__ == "__main__":
    print("Lipidomics MRM Analysis Module Loaded")
    print("Use run_full_analysis() to process mzML files")