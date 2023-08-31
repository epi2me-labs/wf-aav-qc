"""Test the truncations data generation code."""

import tempfile
from unittest.mock import Mock

import pandas as pd
from workflow_glue.truncations import main


def test_main():
    """Test main."""
    annotation_str =  \
        """aav8\t11\t 156\tITR1\n
aav8\t175\t891\tcmv_enhancer\n
aav8\t897\t1099\tcmv_promoter\n
aav8\t1223\t1941\tGFP\n
aav8\t2324\t2454\tITR2\n
aav8\t3252\t3347\tAmpR_promoter\n
aav8\t3350\t4206\tAmpR\n
aav8\t4383\t4920\tori
        """

    annotation = tempfile.NamedTemporaryFile()
    with open(annotation.name, 'w') as fh_a:
        fh_a.write(annotation_str)

    # Note: bam_info will contain more columns, but only three are used.
    bam_info_entries = [
        ['Read',  'Ref',       'Pos',  'EndPos'],
        # Contained fully in ITR-ITR
        ['read1', 'flip_flop', '20',   '2000'],
        # Not fully in ITR-ITR, but within padded region
        ['read2', 'flop_flop', '10',   '2000'],
        # Not fully in ITR-ITR
        ['read3', 'flop_flop', '2500', '4000'],
    ]

    bam_info = tempfile.NamedTemporaryFile()
    with open(bam_info.name, 'w') as fh_b:
        for line in bam_info_entries:
            fh_b.write("\t".join(line) + '\n')

    outfile = tempfile.NamedTemporaryFile(mode='w')
    args = Mock()
    args.annotation = annotation.name
    args.bam_info = bam_info.name
    args.sample_id = "test_id"
    args.outfile = outfile.name
    main(args)

    result_df = pd.read_csv(args.outfile, sep='\t')

    # Two results where fully in ITR-ITR should be returned
    assert len(result_df) == 2
    # First results is from read1
    assert result_df.iloc[0].values.tolist() == [20, 2000, 'test_id']
    # Second resultsi s from read2
    assert result_df.iloc[1].values.tolist() == [10, 2000, 'test_id']
