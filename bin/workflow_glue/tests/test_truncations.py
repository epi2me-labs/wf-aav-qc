"""Test the truncations data generation code."""

import tempfile
from unittest.mock import Mock

import pandas as pd
from workflow_glue.truncations import main


def test_main():
    """Test main."""
    itr1_start = 11
    itr2_end = 2454

    # Note: bam_info will contain more columns, but only three are used in the tested
    # script.
    bam_info_entries = [
        ['Read',  'Ref',  'Pos',  'EndPos'],
        # Contained fully in ITR-ITR
        ['read1', 'aav8', '20',   '2000'],
        # Not fully in ITR-ITR, but within padded region
        ['read2', 'aav8', '10',   '2000'],
        # Not fully in ITR-ITR
        ['read3', 'aav8', '2500', '4000'],
        # A non-transgene plasmid alignment
        ['read4', 'cell_line', '300', '3000']
    ]

    bam_info = tempfile.NamedTemporaryFile()
    with open(bam_info.name, 'w') as fh_b:
        for line in bam_info_entries:
            fh_b.write("\t".join(line) + '\n')

    outfile = tempfile.NamedTemporaryFile(mode='w')
    args = Mock()
    args.itr_range = [itr1_start, itr2_end]
    args.bam_info = bam_info.name
    args.sample_id = "test_id"
    args.outfile = outfile.name
    args.transgene_plasmid_name = 'aav8'
    main(args)

    result_df = pd.read_csv(args.outfile, sep='\t')

    # Two results where fully in ITR-ITR should be returned
    assert len(result_df) == 2
    # First results is from read1
    assert result_df.iloc[0].values.tolist() == [20, 2000, 'test_id']
    # Second results from read2
    assert result_df.iloc[1].values.tolist() == [10, 2000, 'test_id']
