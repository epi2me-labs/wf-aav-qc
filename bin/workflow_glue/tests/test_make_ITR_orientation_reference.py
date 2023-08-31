"""Test ITR orientation generation code."""

from pathlib import Path
import tempfile

import pysam
from workflow_glue import make_ITR_orientation_reference


def test_make_flip_flop_aavs():
    """Test the plasmid oirientation generation."""
    test_data_dr = Path(__file__).resolve().parent / 'data'

    # Input is a bed file containing transgene plasmid features and AAV8 reference fasta
    annotations = test_data_dr / 'AAV8_new.bed'
    input_aav_fasta = test_data_dr / 'AAV8_transgene_fixed_ITR1_no_3bp_del.fasta'

    # Expected output is a fasta file containing the 4 different possible configurations
    # flip-flop flip-flip etc.
    # Note: the expected fasta file is simply the output from the original Snakemake
    # pipeline. At the moment we're just checkig that the nextflow workflow does not
    # deviate.
    expected = test_data_dr / 'AAV8_rep2_transgene_plasmid_ITR-ITR_orientations.fasta'

    outpath = tempfile.NamedTemporaryFile(suffix='.fa')
    make_ITR_orientation_reference.make_flip_flop_aavs(
        input_aav_fasta, annotations, outpath.name)

    for exp, tst in zip(
        pysam.FastxFile(
            str(expected)), pysam.FastxFile(
            outpath.name)):
        assert (exp.name == tst.name)
        assert (exp.sequence == tst.sequence)
