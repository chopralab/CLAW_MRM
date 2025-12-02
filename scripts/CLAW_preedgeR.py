#!/usr/bin/env python3
"""
Create pre-EdgeR formatted CSV from intensity data and sample labels.

This script:
1. Takes intensity data (with Lipid, Class, Sample_ID, Intensity columns)
2. Takes a labels file mapping samples to groups (e.g., Genotype: WT vs 5XFAD)
3. Outputs a wide-format CSV suitable for EdgeR analysis

Example usage:
  python CLAW_preedgeR.py \\
      --data results.csv \\
      --labels labels.csv \\
      --sample-col "Sample Name" \\
      --factor-col Genotype \\
      --group1 5XFAD \\
      --group2 WT \\
      --out-prefix preedger_output

Output format:
  Lipid, Class, <sample1>, <sample2>, ..., Title1, Title2, Title, length1, length2, Blank_name
"""

import argparse
import pandas as pd
from pathlib import Path
import re


def extract_sample_name_from_id(sample_id, labels_sample_names):
    """
    Extract the sample name from a Sample_ID by matching against known labels.
    
    Sample_ID format example: PC_FAD199_m1_66_B_N1_loading_size_2_13_25
    Sample name format example: FAD199_m1_66_B
    """
    sample_id_str = str(sample_id)
    
    # Try to find a matching sample name from labels
    for label_name in labels_sample_names:
        if label_name in sample_id_str:
            return label_name
    
    return None


def build_preedger_matrix(data_df, labels_df, sample_col, factor_col, group1, group2):
    """
    Build a pre-EdgeR formatted matrix with intensity values.
    
    Args:
        data_df: DataFrame with Lipid, Class, Sample_ID, Intensity columns
        labels_df: DataFrame mapping sample names to groups
        sample_col: Column name in labels_df containing sample names
        factor_col: Column name in labels_df for grouping (e.g., Genotype)
        group1: First group (treatment), e.g., 5XFAD
        group2: Second group (control), e.g., WT
    
    Returns:
        DataFrame in pre-EdgeR format
    """
    # Get sample names from labels that belong to either group
    labels_df = labels_df.copy()
    labels_df[factor_col] = labels_df[factor_col].astype(str)
    
    mask = labels_df[factor_col].isin([str(group1), str(group2)])
    group_labels = labels_df.loc[mask].copy()
    
    if group_labels.empty:
        raise ValueError(
            f"No samples found where {factor_col} is {group1!r} or {group2!r}. "
            f"Available values: {labels_df[factor_col].unique().tolist()}"
        )
    
    # Get sample names for each group
    group1_samples = group_labels.loc[
        group_labels[factor_col] == str(group1), sample_col
    ].tolist()
    group2_samples = group_labels.loc[
        group_labels[factor_col] == str(group2), sample_col
    ].tolist()
    
    # Find blank sample if exists
    blank_mask = labels_df[factor_col].str.lower().str.contains('blank', na=False)
    blank_samples = labels_df.loc[blank_mask, sample_col].tolist() if blank_mask.any() else []
    blank_name = blank_samples[0] if blank_samples else ""
    
    print(f"  Group1 ({group1}) samples: {group1_samples}")
    print(f"  Group2 ({group2}) samples: {group2_samples}")
    print(f"  Blank sample: {blank_name}")
    
    # Map Sample_ID to sample name from labels
    all_label_names = labels_df[sample_col].tolist()
    data_df = data_df.copy()
    data_df['matched_sample'] = data_df['Sample_ID'].apply(
        lambda x: extract_sample_name_from_id(x, all_label_names)
    )
    
    # Filter to only samples in our groups (or blank)
    valid_samples = group1_samples + group2_samples + blank_samples
    data_filtered = data_df[data_df['matched_sample'].isin(valid_samples)].copy()
    
    if data_filtered.empty:
        raise ValueError(
            "No data rows matched the sample names from labels. "
            "Check that Sample_ID values contain the sample names from labels."
        )
    
    # Pivot: rows = (Lipid, Class), columns = matched_sample, values = Intensity
    # Average intensity if there are multiple entries per lipid per sample
    pivot_df = data_filtered.pivot_table(
        index=['Lipid', 'Class'],
        columns='matched_sample',
        values='Intensity',
        aggfunc='mean'
    ).reset_index()
    
    # Reorder columns: Lipid, Class, group1_samples, group2_samples, blank
    ordered_cols = ['Lipid', 'Class']
    for s in group1_samples:
        if s in pivot_df.columns:
            ordered_cols.append(s)
    for s in group2_samples:
        if s in pivot_df.columns:
            ordered_cols.append(s)
    if blank_name and blank_name in pivot_df.columns:
        ordered_cols.append(blank_name)
    
    # Keep only existing columns
    ordered_cols = [c for c in ordered_cols if c in pivot_df.columns]
    pivot_df = pivot_df[ordered_cols]
    
    # Count samples in each group that are actually in the data
    n_group1 = len([s for s in group1_samples if s in pivot_df.columns])
    n_group2 = len([s for s in group2_samples if s in pivot_df.columns])
    
    # Add metadata columns
    pivot_df['Title1'] = f"Group: {group1}"
    pivot_df['Title2'] = f"Group: {group2}"
    pivot_df['Title'] = f"Group: {group1} vs Group: {group2}"
    pivot_df['length1'] = n_group1
    pivot_df['length2'] = n_group2
    pivot_df['Blank_name'] = blank_name
    
    return pivot_df


def main():
    parser = argparse.ArgumentParser(
        description="Create pre-EdgeR formatted CSV from intensity data and sample labels."
    )
    parser.add_argument(
        "--data",
        required=True,
        help="Path to intensity data CSV (must contain Lipid, Class, Sample_ID, Intensity columns).",
    )
    parser.add_argument(
        "--labels",
        required=True,
        help="Path to labels CSV mapping samples to groups.",
    )
    parser.add_argument(
        "--sample-col",
        default="Sample Name",
        help="Column name in labels file with sample identifiers (default: 'Sample Name').",
    )
    parser.add_argument(
        "--factor-col",
        required=True,
        help="Column name in labels file for grouping (e.g., Genotype, Group).",
    )
    parser.add_argument(
        "--group1",
        required=True,
        help="First group level (treatment), e.g., 5XFAD.",
    )
    parser.add_argument(
        "--group2",
        required=True,
        help="Second group level (control), e.g., WT.",
    )
    parser.add_argument(
        "--out-prefix",
        default=None,
        help="Prefix/path for output file (default: <factor>_<group1>_vs_<group2>).",
    )

    args = parser.parse_args()

    # Validate input files
    data_path = Path(args.data)
    labels_path = Path(args.labels)
    
    if not data_path.exists():
        raise SystemExit(f"Data file not found: {data_path}")
    if not labels_path.exists():
        raise SystemExit(f"Labels file not found: {labels_path}")

    # Load data
    print(f"\nLoading data from: {data_path}")
    data_df = pd.read_csv(data_path)
    
    # Validate required columns in data
    required_data_cols = ['Lipid', 'Class', 'Sample_ID', 'Intensity']
    missing_data_cols = [c for c in required_data_cols if c not in data_df.columns]
    if missing_data_cols:
        raise SystemExit(
            f"Missing columns in data file: {missing_data_cols}. "
            f"Available columns: {list(data_df.columns)}"
        )
    
    # Load labels
    print(f"Loading labels from: {labels_path}")
    labels_df = pd.read_csv(labels_path)
    
    # Validate required columns in labels
    for col in [args.sample_col, args.factor_col]:
        if col not in labels_df.columns:
            raise SystemExit(
                f"Column {col!r} not found in labels file. "
                f"Available columns: {list(labels_df.columns)}"
            )

    factor = args.factor_col
    g1 = args.group1
    g2 = args.group2

    # Determine output path
    out_prefix = args.out_prefix
    if out_prefix is None:
        out_prefix = f"{factor}_{g1}_vs_{g2}"
    
    out_path = Path(f"{out_prefix}_preedger.csv")

    print(f"\nBuilding pre-EdgeR matrix...")
    print(f"  Factor column: {factor}")
    print(f"  Comparison: {g1} vs {g2}")
    
    # Build the pre-EdgeR matrix
    preedger_df = build_preedger_matrix(
        data_df=data_df,
        labels_df=labels_df,
        sample_col=args.sample_col,
        factor_col=factor,
        group1=g1,
        group2=g2,
    )

    # Save output
    preedger_df.to_csv(out_path, index=False)

    print("")
    print("=" * 70)
    print("Pre-EdgeR file created successfully!")
    print("=" * 70)
    print(f"  Data file:      {data_path}")
    print(f"  Labels file:    {labels_path}")
    print(f"  Factor column:  {factor}")
    print(f"  Comparison:     {g1} (treatment) vs {g2} (control)")
    print("")
    print(f"  Output file:    {out_path}")
    print(f"  Total lipids:   {len(preedger_df)}")
    print(f"  Group 1 count:  {preedger_df['length1'].iloc[0] if len(preedger_df) > 0 else 0}")
    print(f"  Group 2 count:  {preedger_df['length2'].iloc[0] if len(preedger_df) > 0 else 0}")
    print("=" * 70)
    print("")


if __name__ == "__main__":
    main()
