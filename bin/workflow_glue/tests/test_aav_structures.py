"""Test the AAV structures data generation code.

The AAV structures process is formed of two main parts
1: Assigning an AlnType to each alignment based on its start and end positions in
relation to the ITR sequences
2: Assigning a ReadType by looking at the individual AlnTypes of the read
"""
from collections import defaultdict
from pathlib import Path
import sys

import pandas as pd
import polars as pl
from polars.testing import assert_frame_equal
import pytest
from workflow_glue.aav_structures import (
    annotate_reads,
    assign_genome_types_to_alignments,
    get_subtype_definitions,
)
import yaml


def assert_equal_view_diffs(expected, actual, **kwargs):
    """Assert two dataframes are equal, printing differences if not."""
    defaults = {
        "check_dtype": False,
        "check_row_order": False,
        "check_column_order": False
    }
    options = {**defaults, **kwargs}
    pl.Config.set_tbl_rows(100)
    try:
        assert_frame_equal(expected, actual, **options)
    except AssertionError:
        diff = expected.join(actual, on='Read', how="anti")
        sys.stdout.write("In expected, not in actual:\n")
        sys.stdout.write(str(diff) + '\n')

        diff2 = actual.join(expected, on='Read', how="anti")
        sys.stdout.write("In actual, not in expected:\n")
        sys.stdout.write(str(diff2))

        sys.stdout.write("\nFull expected:\n")
        sys.stdout.write(str(expected) + '\n')

        sys.stdout.write("\nFull actual:\n")
        sys.stdout.write(str(actual) + '\n')

        raise


@pl.api.register_dataframe_namespace("pd")
class PolarsPandasDebug:
    """For debugging purposes."""

    def __init__(self, df):
        """Set a pandas version of a polars dataframe for debugging."""
        self.df = df.to_pandas()


@pytest.fixture
def test_data(request):
    """Define data location fixture."""
    return Path(request.config.getoption("--test_data")) / "workflow_glue"


@pytest.fixture
def bam_info_df(test_data):
    """
    Build a DataFrame of alignment summaries.

    The bam alignment summaries from seqkit are used in the workflow to
    determine the genome structures present. The function builds a test DataFrame
    with only the required columns present.

    The alignment DataFrame is populated from entries in structures.yaml,
    where examples of each aln_type are given.

    This function takes each example and uses it to populate a pandas dataframe
    similar to what would be expected by the workflow.
    """
    structures_config_file = test_data / 'structures.yaml'
    with open(structures_config_file, 'r') as fh:
        cfg = yaml.load(fh, Loader=yaml.FullLoader)

    records = []
    read_id = 0
    for alignment_type in cfg['aln_types']:
        for example in alignment_type['examples']:

            records.append(
                [
                    read_id, example.get('reference', 'flip_flop'), example['start'],
                    example['end'], example['strand'], alignment_type['aln_type']
                ])
            read_id += 1

    df = pd.DataFrame.from_records(records)
    header = ['Read', 'Ref', 'Pos', 'EndPos', 'Strand', 'expected_aln_type']
    df.columns = header
    df = df.astype({
        'Read': str,
        'Ref': str,
        'Pos': int,
        'EndPos': int
    })
    dfp = pl.from_pandas(df)

    return dfp


@pytest.fixture
def genome_type_df(test_data):
    """Build a test DataFrame containing alignment assignments.

    With columns
    - Read
    - Ref
    - Pos
    - EndPos
    - strand
    - AlnType # This is the expected G
    """
    # Each AAV genome type has one or more example reads in structures.yaml
    # Load these into a dataframe
    structures_config_file = test_data / 'structures.yaml'
    with open(structures_config_file, 'r') as fh:
        cfg = yaml.load(fh, Loader=yaml.FullLoader)

    records = []
    read_to_type = []

    for genome_type in cfg['genome_types']:
        for example in genome_type['examples']:
            for aln in example:
                records.append(
                    [
                        aln['read_id'], aln.get('reference', 'flip_flop'), aln['start'],
                        aln['end'], aln['strand'], aln['aln_type']
                    ])
                read_to_type.append([
                    aln['read_id'], genome_type['genome_subtype'],
                    genome_type['genome_type']])

    # Build the DataFrame that the test alignments will be generated from
    df_alignments = pd.DataFrame.from_records(records)
    header = [
        'Read', 'Ref', 'Pos', 'EndPos', 'Strand', 'aln_type']
    df_alignments.columns = header
    df_alignments = df_alignments.astype({
        'Read': str,
        'Ref': str,
        'Pos': int,
        'EndPos': int,
        'aln_type': str
    })

    df_read_subtype_map = (
        pd.DataFrame.from_records(
            read_to_type,
            columns=['Read', 'Assigned_genome_subtype', 'Assigned_genome_type'])
    )
    return pl.from_pandas(df_alignments), pl.from_pandas(df_read_subtype_map), cfg


def test_assign_genome_types_to_alignments(bam_info_df):
    """Test the initial classification of alignments."""
    itr1_start = 11
    itr1_end = 156
    itr2_start = 2324
    itr2_end = 2454
    itr_fl_threshold = 100
    itr_backbone_threshold = 20

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


def test_annotate_reads(genome_type_df):
    """Test assigning Assigned_genome_subtype to reads."""
    df_alns, read_to_subtype_map_df, cfg = genome_type_df

    actual_per_read_subtypes, actual_summary = annotate_reads(
        'sample', df_alns, get_subtype_definitions(), symmetry_threshold=10)
    actual_per_read_subtypes = actual_per_read_subtypes.sort('Read')

    expected_per_read_info = (
        read_to_subtype_map_df
        .unique(subset='Read')
        .sort('Read'))

    assert_equal_view_diffs(
        expected_per_read_info,
        actual_per_read_subtypes.select(
            ['Read', 'Assigned_genome_subtype', 'Assigned_genome_type'])
    )

    actual_summary = actual_summary.sort('Assigned_genome_type')
    expected_summary = defaultdict(list)
    for type_ in cfg['genome_types']:
        expected_summary['Assigned_genome_type'].append(type_['genome_type'])
        expected_summary['count'].append(len(type_['examples']))

    expected_summary_df = pl.from_dict(expected_summary)
    expected_summary_df = (
        expected_summary_df
        .group_by('Assigned_genome_type')
        .sum()
        .with_columns(
            (pl.col("count") / pl.sum("count") * 100)
            .round(2).alias("percentage"),
            pl.lit("sample").alias("sample_id"))
    ).sort('Assigned_genome_type')

    assert_equal_view_diffs(expected_summary_df, actual_summary)
