"""Test the AAV structures data generation code.

The workflow can be broken down into two main parts
1: Assigning a Genome_type to each alignment based on its start and end positions in
relation to the ITR sequences
2: Assignng Assigned_genome_subtype which is derived from the Genome_subtypes of the
read's alignments ....
"""

from pathlib import Path

import pandas as pd
import polars as pl
from workflow_glue.aav_structures import (
    annotate_reads,
    assign_genome_types_to_alignments
)
import yaml


def build_bam_info():
    """
    Build alignment summary DatFrame.

    The bam info alignment results (seqkit stats) are used in the workflow to
    determine the genome structures present.

    In test/data/structures.yaml, each genome type is defined and examples are given.
    This function takes each example and generates a pandas dataframe similar to
    what would be expected by the workflow.
    """
    test_data_dr = Path(__file__).resolve().parent / 'data'
    structures_config_file = test_data_dr / 'structures.yaml'
    with open(structures_config_file, 'r') as fh:
        cfg = yaml.load(fh, Loader=yaml.FullLoader)

    records = []
    read_id = 0
    for gt_data in cfg['genome_types']:
        for aln in gt_data['examples']:
            records.append(
                [read_id] + aln.split() + [gt_data['genome_type']], )
            read_id += 1

    df = pd.DataFrame.from_records(records)
    header = ['Read', 'Ref', 'Pos', 'EndPos', 'Strand', 'IsSec', 'expected_genome_type']
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
    """Build a dataframe containing alignments that have been asigned a Genome_type.

    With columns
    - Read
    - Ref
    - Pos
    - EndPos
    - Strand
    - GenomeType
    """
    test_data_dir = Path(__file__).resolve().parent / 'data'
    structures_config_file = test_data_dir / 'structures.yaml'
    with open(structures_config_file, 'r') as fh:
        cfg = yaml.load(fh, Loader=yaml.FullLoader)

    records = []
    read_subtype_map = {}  # Reads to expected subtype mapping
    read_id = 0

    for subtype in cfg['assigned_genome_subtypes']:
        read_subtype_map[f'read_{read_id}'] = subtype['subtype']
        for aln, genome_type in subtype['alignments']:
            records.append(
                [f'read_{read_id}'] + aln.split() + [genome_type])
            read_id += 1

    df = pd.DataFrame.from_records(records)
    header = [
        'Read', 'Ref', 'Pos', 'EndPos', 'Strand',
        'IsSec', 'Genome_type']
    df.columns = header
    df = df.astype({
        'Read': str,
        'Ref': str,
        'Pos': int,
        'EndPos': int,
        'IsSec': int,
        'Genome_type': str
    })

    df_read_subtype_map = (
        pd.DataFrame.from_dict(read_subtype_map, orient='index')
        .reset_index(drop=False)
    )
    df_read_subtype_map.columns = ['Read', 'Assigned_genome_subtype']
    dfp = pl.from_pandas(df)
    return dfp, df_read_subtype_map


def test_assign_genome_types_to_alignments():
    """Test the initial classification of alignments."""
    itr1_start = 11
    itr1_end = 156
    itr2_start = 2324
    itr2_end = 2454
    padding = 10

    cfg, bam_info_df = build_bam_info()

    result_df = assign_genome_types_to_alignments(
        bam_info_df,
        itr1_start,
        itr1_end,
        itr2_start,
        itr2_end,
        padding
    )

    result_df = result_df.to_pandas()[['Read', 'Genome_type']]
    test_df = bam_info_df.to_pandas()[['Read', 'expected_genome_type']]
    test_df.rename(columns={'expected_genome_type': 'Genome_type'}, inplace=True)

    pd.testing.assert_frame_equal(test_df, result_df)


def test_annotate_reads():
    """Test annotationg reads.

    # Cols that we need in trans_df
    - Read Ref Pos EndPos Strand GenomeType

    TODO:
        - Add more all the assigned_genome_type entries into the structures.yaml
        - Test assigned_genome_types_summary
    """
    df_subtype, read_to_subtype_map_df = build_genome_type_df()
    # For this test we convert the expected_result column to 'GenomeType' we know
    # GenomeType is correct as the `test_assign_genome_types_to_alignments` will
    # Have passed to get to this point

    assigned_genome_types_summary, per_read_output = annotate_reads(df_subtype)
    pd.testing.assert_frame_equal(
        read_to_subtype_map_df,
        per_read_output.to_pandas()[['Read', 'Assigned_genome_subtype']]
    )
