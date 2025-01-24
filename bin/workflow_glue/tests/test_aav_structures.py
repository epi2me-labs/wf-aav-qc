"""Test the AAV structures data generation code.

The AAV structures process is formed of two main parts
1: Assigning an AlnType to each alignment based on its start and end positions in
relation to the ITR sequences
2: Assigning a ReadType by looking at the individual AlnTypes of the read
"""
from pathlib import Path

import pandas as pd
import polars as pl
from polars.testing import assert_frame_equal
from workflow_glue.aav_structures import (
    annotate_reads,
    assign_genome_types_to_alignments,
    get_subtype_definitions,
)
import yaml


def build_bam_info():
    """
    Build a DataFrame of alignment summaries.

    The bam alignment summaries from seqkit are used in the workflow to
    determine the genome structures present. The function builds a test DataFrame
    with only the required columns present.

    The alignment DataFrame is populated from entries in tests/data/structures.yaml,
    where examples of each aln_type are given.

    This function takes each example and uses it to populate a pandas dataframe
    similar to what would be expected by the workflow.
    """
    test_data_dr = Path(__file__).resolve().parent / 'data'
    structures_config_file = test_data_dr / 'structures.yaml'
    with open(structures_config_file, 'r') as fh:
        cfg = yaml.load(fh, Loader=yaml.FullLoader)

    records = []
    read_id = 0
    for gt_data in cfg['aln_types']:
        for aln in gt_data['examples']:
            records.append(
                [read_id] + aln.split() + [gt_data['aln_type']], )
            read_id += 1

    df = pd.DataFrame.from_records(records)
    header = ['Read', 'Ref', 'Pos', 'EndPos', 'Strand', 'IsSec', 'expected_aln_type']
    df.columns = header
    df = df.astype({
        'Read': str,
        'Ref': str,
        'Pos': int,
        'EndPos': int,
        'IsSec': int
    })
    dfp = pl.from_pandas(df)

    return cfg, dfp


def build_genome_type_df():
    """Build a test DataFrame containing alignment assignments.

    With columns
    - Read
    - Ref
    - Pos
    - EndPos
    - strand
    - IsSec
    - AlnType # This is the expected G
    """
    test_data_dir = Path(__file__).resolve().parent / 'data'
    # Each AAV genome type has one or more example reads in structures.yaml
    # Load these into a dataframe
    structures_config_file = test_data_dir / 'structures.yaml'
    with open(structures_config_file, 'r') as fh:
        cfg = yaml.load(fh, Loader=yaml.FullLoader)

    records = []
    read_to_type = []

    for read_type in cfg['genome_types']:
        for example in read_type['examples']:
            for aln, aln_type in example:
                aln_info = aln.split()
                read_id = aln_info[0]
                records.append(
                    [read_id] + aln_info[1:] + [aln_type])
                read_to_type.append([
                    read_id, read_type['genome_subtype'], read_type['genome_type']])

    # Build the DataFrame that the test alignments will be generated from
    df_alignments = pd.DataFrame.from_records(records)
    header = [
        'Read', 'Ref', 'Pos', 'EndPos', 'Strand',
        'IsSec', 'aln_type']
    df_alignments.columns = header
    df_alignments = df_alignments.astype({
        'Read': str,
        'Ref': str,
        'Pos': int,
        'EndPos': int,
        'IsSec': int,
        'aln_type': str
    })

    df_read_subtype_map = (
        pd.DataFrame.from_records(
            read_to_type,
            columns=['Read', 'Assigned_genome_subtype', 'Assigned_genome_type'])
    )
    return pl.from_pandas(df_alignments), df_read_subtype_map


def test_assign_genome_types_to_alignments():
    """Test the initial classification of alignments."""
    itr1_start = 11
    itr1_end = 156
    itr2_start = 2324
    itr2_end = 2454
    itr_fl_threshold = 100
    itr_backbone_threshold = 20

    _, bam_info_df = build_bam_info()

    result_df = assign_genome_types_to_alignments(
        bam_info_df,
        itr1_start,
        itr1_end,
        itr2_start,
        itr2_end,
        itr_fl_threshold,
        itr_backbone_threshold
    )

    actual = result_df.select(['Read', 'aln_type'])
    expected = (
        bam_info_df.select(['Read', 'expected_aln_type'])
        .rename({'expected_aln_type': 'aln_type'})
    )

    assert_frame_equal(
        expected, actual,
        check_dtype=False, check_row_order=False, check_column_order=False
    )


def test_annotate_reads():
    """Test assigning Assigned_genome_subtype to reads."""
    df_alns, read_to_subtype_map_df = build_genome_type_df()

    actual_per_read_subtypes, actual_summary = annotate_reads(
        'sample', df_alns, get_subtype_definitions(), symmetry_threshold=10)

    expected_per_read_info = (
        pl.from_pandas(read_to_subtype_map_df)
        .unique(subset='Read')
        .sort('Read'))

    assert_frame_equal(
        expected_per_read_info,
        actual_per_read_subtypes.select(
            ['Read', 'Assigned_genome_subtype', 'Assigned_genome_type']),
        check_dtype=False, check_row_order=False, check_column_order=False
    )

    expected_summary = pl.from_pandas(pd.DataFrame({
        "Assigned_genome_type": [
            "Full_ssAAV",
            "Full_scAAV",
            "Partial ssAAV",
            "Partial scAAV",
            "ITR region only",
            "Complex",
            "Backbone contamination"
        ],
        "count": [1, 2, 10, 19, 3, 2, 3],
        "percentage": [2.5, 5.0, 25.0, 47.5, 7.5, 5.0, 7.5],
        "sample_id": ["sample"] * 7,
    }))

    assert_frame_equal(
        expected_summary, actual_summary,
        check_dtype=False, check_row_order=False, check_column_order=False
    )
