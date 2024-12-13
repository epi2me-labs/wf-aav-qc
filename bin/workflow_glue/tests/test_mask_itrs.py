"""Test ITR orientation generation code."""

import tempfile

from workflow_glue.mask_itrs import (
    get_variable_itr_regions, mask_itrs
)


def test_get_variable_itr_regions():
    """Test getting the locations of the variable positions in a an ITR sequence."""
    # Original ITR orientation is Flip, alternate is flop (varaible in lower case)
    original_itr = "CGccCGggC-AAA-GccCGggCG-T-CGggCGacC-TTT-GgtCGccCG"
    # altern_itr = "CGggCGacC-AAA-GgtCGccCG-A-CGccCGggC-TTT-GccCGggCG"
    # 0-based
    variable_locs = [2, 3, 6, 7, 13, 14, 17, 18, 21, 24, 25, 28, 29, 35, 36, 39, 40]
    # Make the ITR - add 10nt prefix as real sequence will contain ITR upstream bases
    prefix_len = 10
    itr_to_test = f"{'G' * prefix_len}{original_itr.replace('-', '')}"
    # Add the prefix len to the expected positions
    variable_locs = [x + prefix_len for x in variable_locs]

    found_locations = get_variable_itr_regions(itr_to_test)
    assert found_locations == variable_locs


def test_make_flip_flop_aavs():
    """Test the whole of the plasmid masking script."""
    # A contrived ITR sequence B-C region where the end of each element is variable
    flip = 'CCCCCCCCg-AAA-cGGGGGGGG-T-CCCCCCCCc-TTT-gGGGGGGGG'.replace('-', '')
    flop = 'CCCCCCCCc-AAA-gGGGGGGGG-A-CCCCCCCCg-TTT-cGGGGGGGG'.replace('-', '')
    # Variable locations within the ITR variable region  (5 positions)
    itr_b_c_var_locations = [8, 12, 21, 30, 34]

    itr_padding = 10
    midsection_len = 10

    flip = f"{'G' * itr_padding}{flip}{'G' * itr_padding}"
    flop = f"{'G' * itr_padding}{flop}{'G' * itr_padding}"

    flip_flop_plasmid_seq = (
        f"{flip}{'C' * midsection_len}{flop}"
    )

    # Create flip-flop transgene cassette with some N padding
    flip_flop_reference = (
        f">ref_name\n{flip_flop_plasmid_seq}\n")

    with tempfile.NamedTemporaryFile(
            mode='w+', suffix='.fasta', delete=False) as fh_fa:
        fh_fa.write(flip_flop_reference)

    # Define the ITR positions
    itr1_start = itr_padding
    itr1_end = itr1_start + len(flip)
    itr2_start = itr1_end + midsection_len
    itr2_end = itr2_start + len(flop)

    masked_seq = mask_itrs(
        fh_fa.name,
        [itr1_start, itr1_end, itr2_start, itr2_end]
    )

    # For ITR1, get the known variable regions.
    # Extract the corresponding bases from the masked sequence bases.
    # The 5 bases should be 'N'
    itr1_mask_locations = [x + itr1_start for x in itr_b_c_var_locations]
    itr1_mask_bases = ''.join(
        [masked_seq[x] for x in itr1_mask_locations])
    assert itr1_mask_bases == 'NNNNN'

    # Do the same for ITR2
    itr2_mask_locations = [x + itr2_start for x in itr_b_c_var_locations]
    itr2_mask_bases = ''.join(
        [masked_seq[x] for x in itr2_mask_locations])
    assert itr2_mask_bases == 'NNNNN'

    # Check that none of the other bases are masked
    # There should be 10 masked bases; 5 for each ITR
    assert masked_seq.count('N') == 10
