# script to format data as csv files to input into biopan
#!/usr/bin/env python3
"""
Build a grouped counts matrix for 2 groups (e.g. 5xFAD vs WT).

Inputs
------
1) labels CSV (metadata) with:
   - a sample ID column (e.g. Sample_ID or Sample Name)
   - a factor column (e.g. Genotype) with values like 5xFAD / WT

2) long-format intensity CSV with:
   - sample ID column (same name as in labels)
   - lipid ID column (e.g. Lipid)
   - intensity column (e.g. Intensity)

Output
------
CSV with:
  Lipid,<Group1_col>...,<Group2_col>...

For example:

  Lipid,Group_B,Group_B,Group_B,Group_B,Group__A,Group__A,Group__A,Group__A
  1-O-behenoyl-Cer(...),1179...,1124..., ...

By default:
  - group1 (e.g. 5xFAD) -> columns named "Group_B"
  - group2 (e.g. WT)   -> columns named "Group__A"

You can override these names with --group1-col-name / --group2-col-name.
"""

import argparse
from pathlib import Path

import pandas as pd


def build_group_matrix(
    labels_df: pd.DataFrame,
    intensities_df: pd.DataFrame,
    sample_col: str,
    factor_col: str,
    group1: str,
    group2: str,
    lipid_col: str = "Lipid",
    intensity_col: str = "Intensity",
    group1_col_name: str = "Group_B",
    group2_col_name: str = "Group__A",
) -> pd.DataFrame:
    """
    Build a grouped counts matrix.

    labels_df:
        Metadata with at least [sample_col, factor_col].
    intensities_df:
        Long-format data with at least [sample_col, lipid_col, intensity_col].
    sample_col:
        Column name for sample IDs (must exist in BOTH dataframes).
    factor_col:
        Column in labels_df defining the group (e.g. Genotype).
    group1, group2:
        Two levels of factor_col (e.g. '5xFAD', 'WT').
    group1_col_name, group2_col_name:
        Column names to use for group1 and group2 replicates, respectively.

    Returns
    -------
    out_df : pd.DataFrame
        Columns: Lipid, <group1_col_name>..., <group2_col_name>...
    """
    # ------------------------------------------------------------------
    # 1. Filter labels to the two groups and get sample lists
    # ------------------------------------------------------------------
    mask = labels_df[factor_col].astype(str).isin([str(group1), str(group2)])
    meta = labels_df.loc[mask].copy()

    if meta.empty:
        raise ValueError(
            f"No rows in labels where {factor_col} is {group1!r} or {group2!r}."
        )

    if sample_col not in meta.columns:
        raise ValueError(
            f"sample_col={sample_col!r} not found in labels. "
            f"Available columns: {list(labels_df.columns)}"
        )

    # Drop NA sample IDs
    meta = meta.dropna(subset=[sample_col])
    if meta.empty:
        raise ValueError(
            f"After dropping NA in {sample_col!r}, no samples remain in labels."
        )

    # Sample lists for each group (preserve order in labels file)
    group1_samples = (
        meta.loc[meta[factor_col].astype(str) == str(group1), sample_col]
        .astype(str)
        .tolist()
    )
    group2_samples = (
        meta.loc[meta[factor_col].astype(str) == str(group2), sample_col]
        .astype(str)
        .tolist()
    )

    if len(group1_samples) == 0 or len(group2_samples) == 0:
        raise ValueError(
            f"Need at least one sample in each group. "
            f"Found {len(group1_samples)} in {group1!r}, "
            f"{len(group2_samples)} in {group2!r}."
        )

    # ------------------------------------------------------------------
    # 2. Prepare intensities: pivot to Lipid x Sample matrix
    # ------------------------------------------------------------------
    for col in [sample_col, lipid_col, intensity_col]:
        if col not in intensities_df.columns:
            raise ValueError(
                f"Column {col!r} not found in intensities file. "
                f"Available columns: {list(intensities_df.columns)}"
            )

    df = intensities_df.copy()
    df[sample_col] = df[sample_col].astype(str)

    # Keep only samples that are present in the labels for these groups
    allowed_samples = set(group1_samples + group2_samples)
    df = df[df[sample_col].isin(allowed_samples)].copy()
    if df.empty:
        raise ValueError(
            "No rows in intensity data match the group samples from labels."
        )

    # Aggregate if multiple rows per Lipid / Sample
    grouped = (
        df.groupby([lipid_col, sample_col])[intensity_col]
        .sum()
        .reset_index()
    )

    # Pivot: rows = Lipid, columns = sample IDs
    pivot = grouped.pivot(
        index=lipid_col,
        columns=sample_col,
        values=intensity_col,
    )

    # Make sure columns name is dropped
    pivot.columns.name = None

    # ------------------------------------------------------------------
    # 3. Build the grouped output: Lipid, [Group_B...], [Group__A...]
    # ------------------------------------------------------------------
    # Keep only the samples that exist in the pivot
    group1_cols = [s for s in group1_samples if s in pivot.columns]
    group2_cols = [s for s in group2_samples if s in pivot.columns]

    if len(group1_cols) == 0 or len(group2_cols) == 0:
        raise ValueError(
            "After intersecting with intensity data, one of the groups "
            "has no samples left. Check sample IDs between labels and data."
        )

    # Extract values in the desired sample order
    ordered_cols = group1_cols + group2_cols
    values = pivot[ordered_cols].to_numpy()

    # Create repeated group column names
    new_cols = (
        [group1_col_name] * len(group1_cols)
        + [group2_col_name] * len(group2_cols)
    )

    out_df = pd.DataFrame(values, index=pivot.index, columns=new_cols)

    # Lipid as a proper column
    out_df.insert(0, "Lipid", out_df.index)

    # Replace any NaNs with 0 (optional; change if you prefer NaN)
    out_df = out_df.reset_index(drop=True).fillna(0)

    return out_df


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create a grouped counts matrix (Lipid + Group_B/Group__A columns) "
            "from labels and long-format intensity data."
        )
    )

    parser.add_argument(
        "--labels",
        required=True,
        help="Path to labels/metadata CSV (must contain sample + grouping columns).",
    )
    parser.add_argument(
        "--intensities",
        required=True,
        help="Path to long-format intensity CSV.",
    )
    parser.add_argument(
        "--sample-col",
        default="Sample_ID",
        help=(
            "Column with sample IDs (must exist in BOTH labels and intensities). "
            "Default: Sample_ID"
        ),
    )
    parser.add_argument(
        "--factor-col",
        required=True,
        help="Column in labels defining the group, e.g. Genotype.",
    )
    parser.add_argument(
        "--group1",
        required=True,
        help="First group (e.g. 5xFAD) -> becomes Group_B by default.",
    )
    parser.add_argument(
        "--group2",
        required=True,
        help="Second group (e.g. WT) -> becomes Group__A by default.",
    )

    parser.add_argument(
        "--lipid-col",
        default="Lipid",
        help="Column name for lipid ID in intensities file (default: Lipid).",
    )
    parser.add_argument(
        "--intensity-col",
        default="Intensity",
        help="Column name for numeric intensity in intensities file (default: Intensity).",
    )

    parser.add_argument(
        "--group1-col-name",
        default="Group_B",
        help="Column name to use for group1 replicates (default: Group_B).",
    )
    parser.add_argument(
        "--group2-col-name",
        default="Group__A",
        help="Column name to use for group2 replicates (default: Group__A).",
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Path to output CSV (e.g. .../Genotype_5XFAD_vs_WT_GroupMatrix.csv).",
    )

    args = parser.parse_args()

    labels_path = Path(args.labels)
    intensities_path = Path(args.intensities)
    out_path = Path(args.out)

    if not labels_path.exists():
        raise SystemExit(f"Labels file not found: {labels_path}")
    if not intensities_path.exists():
        raise SystemExit(f"Intensity file not found: {intensities_path}")

    labels_df = pd.read_csv(labels_path)
    intensities_df = pd.read_csv(intensities_path)

    out_df = build_group_matrix(
        labels_df=labels_df,
        intensities_df=intensities_df,
        sample_col=args.sample_col,
        factor_col=args.factor_col,
        group1=args.group1,
        group2=args.group2,
        lipid_col=args.lipid_col,
        intensity_col=args.intensity_col,
        group1_col_name=args.group1_col_name,
        group2_col_name=args.group2_col_name,
    )

    out_df.to_csv(out_path, index=False)

    print("")
    print("=" * 70)
    print("Grouped counts matrix created")
    print("=" * 70)
    print(f"  Labels file:         {labels_path}")
    print(f"  Intensity file:      {intensities_path}")
    print(f"  Sample column:       {args.sample_col}")
    print(f"  Factor column:       {args.factor_col}")
    print(f"  Group1 (-> {args.group1_col_name}): {args.group1}")
    print(f"  Group2 (-> {args.group2_col_name}): {args.group2}")
    print("")
    print(f"  Output:              {out_path}")
    print("=" * 70)
    print("")


if __name__ == "__main__":
    main()
