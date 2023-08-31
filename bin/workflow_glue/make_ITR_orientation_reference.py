#!/usr/bin/env python
"""Create reference with the various flip-flop orientations."""
import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("make_ITR_orientation_reference")

    parser.add_argument(
        '--fasta', help="Fasta file of the transgene plasmid")
    parser.add_argument(
        '--bed', help="Bed file for transgene plasmid")
    parser.add_argument(
        '--outfile', help="Path for output fasta file")

    return parser


def make_flip_flop_aavs(fasta, bed, outfile):
    """Make a fasta file of all 4 possible transgene plasmid orientations."""
    # Read in bed file
    logger = get_named_logger('ITR_ori')

    with open(bed, 'r') as infile:
        for line in infile:
            cols = line.strip("\n").split("\t")
            ref = cols[0]
            if "ITR1" in cols[3]:
                itr1_start_pos = int(cols[1])
                itr1_end_pos = int(cols[2])
            if "ITR2" in cols[3]:
                itr2_start_pos = int(cols[1])
                itr2_end_pos = int(cols[2])

    # Read in sequence
    fasta = pysam.FastaFile(fasta)
    original_length = len(fasta.fetch(ref))

    # Extract ITR sequences
    itr1_seq = fasta.fetch(ref, itr1_start_pos, itr1_end_pos).upper()
    itr2_seq = fasta.fetch(ref, itr2_start_pos, itr2_end_pos).upper()

    # Identify parts of the ITR
    # 4 parts = B1, B2, C1 and C2 (B1 and B2 bind together and same for C1 and
    # C2)

    def get_itr_regions(itr_seq, name):
        index_aaa = itr_seq.find('AAA')
        index_ttt = itr_seq.find('TTT')
        if index_aaa < index_ttt:
            itr_first = 'AAA'
            itr_second = 'TTT'
            first_idx = index_aaa
            second_idx = index_ttt
        else:
            itr_first = 'TTT'
            itr_second = 'AAA'
            first_idx = index_ttt
            second_idx = index_aaa

        itr_end = itr_seq[second_idx + 12: len(itr_seq) + 1]
        itr_between_b_and_c = itr_seq[first_idx + 12]
        itr_start = itr_seq[0: first_idx - 9]  # What is this?

        itr_b1 = itr_seq[index_aaa - 9: index_aaa]
        itr_b2 = itr_seq[index_aaa + 3: index_aaa + 12]
        itr_c1 = itr_seq[index_ttt - 9: index_ttt]
        itr_c2 = itr_seq[index_ttt + 3: index_ttt + 12]

        if itr_between_b_and_c == "A":
            itr_orientation = "flop"
        elif itr_between_b_and_c == "T":
            itr_orientation = "flip"

        logger.info(name)
        logger.info(itr_seq)
        logger.info(f"B1: {itr_b1}")
        logger.info(f"B2: {itr_b2}")
        logger.info(f"between_B_and_C: {itr_between_b_and_c}")
        logger.info(f"C1: {itr_c1}")
        logger.info(f"C2: {itr_c2}")

        return (
            itr_first, itr_second, itr_start, itr_end,
            itr_between_b_and_c, itr_b1, itr_b2, itr_c1, itr_c2, itr_orientation)

    (
        itr1_first, itr1_second, itr1_start, itr1_end,
        itr1_between_b_and_c, itr1_b1, itr1_b2, itr1_c1, itr1_c2, itr1_orientation
    ) = get_itr_regions(itr1_seq, 'ITR2')

    (
        itr2_first, itr2_second, itr2_start, itr2_end,
        itr2_between_b_and_c, itr2_b1, itr2_b2, itr2_c1, itr2_c2, itr2_orientation
    ) = get_itr_regions(itr2_seq, 'ITR2')

    # Checks
    if itr1_between_b_and_c not in ["A", "T"]:
        raise ValueError("Error with ITR1 structure")
    if itr2_between_b_and_c not in ["A", "T"]:
        raise ValueError("Error with ITR2 structure")
    # Add in something if AAA and/or TTT not found

    # Generate flip and flop alignments
    # Note: Do not know which orientation original reference is in
    # so can make the four alignments but not entirely sure which is which.
    # Flip has T between B and C sections, flop has an A between B and C
    # sections

    # Create different orientations
    if itr1_orientation == "flip":
        itr1_flip = itr1_seq
        itr1_flop = (
            f"{itr1_start}{itr1_c1}{itr1_first}{itr1_c2}A"
            f"{itr1_b1}{itr1_second}{itr1_b2}{itr1_end}")

    elif itr1_orientation == "flop":
        # swap A's and T's
        itr1_flip = (
            f"{itr1_start}{itr1_b1}{itr1_first}{itr1_b2}T"
            f"{itr1_c1}{itr1_second}{itr1_c2}{itr1_end}")
        itr1_flop = itr1_seq

    if itr2_orientation == "flip":
        itr2_flip = itr2_seq
        itr2_flop = (
            f"{itr2_start}{itr2_c1}{itr2_first}{itr2_c2}A"
            f"{itr2_b1}{itr2_second}{itr2_b2}{itr2_end}")

    elif itr2_orientation == "flop":
        itr2_flip = (
            f"{itr2_start}{itr2_b1}{itr2_first}{itr2_b2}T"
            f"{itr2_c1}{itr2_second}{itr2_c2}{itr2_end}")
        itr2_flop = itr2_seq

    logger.info(f"itr1_orientation: {itr1_orientation}")
    logger.info(f"itr2_orientation: {itr2_orientation}")
    logger.info(f"itr1_flip:{itr1_flip}, {len(itr1_flip)}")
    logger.info(f"itr1_flop:{itr1_flop}, {len(itr1_flop)}")
    logger.info(f"itr2_flip:{itr2_flip}, {len(itr2_flip)}")
    logger.info(f"itr2_flop:{itr2_flop}, {len(itr2_flop)}")

    # Check ITR lengths are the same
    if len(itr1_flip) != len(itr1_flop):
        raise ValueError("Issue with length of ITR1 orientations")
    if len(itr2_flip) != len(itr2_flop):
        raise ValueError("Issue with length of ITR2 orientations")

    # Make complete plasmid with 4 ITR-ITR orientations
    pre_itr1 = fasta.fetch(ref, 0, itr1_start_pos)
    inter_itr_region = fasta.fetch(ref, itr1_end_pos, itr2_start_pos)
    post_itr2 = fasta.fetch(reference=ref, start=itr2_end_pos)
    flip_flip = f"{pre_itr1}{itr1_flip}{inter_itr_region}{itr2_flip}{post_itr2}"
    flip_flop = f"{pre_itr1}{itr1_flip}{inter_itr_region}{itr2_flop}{post_itr2}"
    flop_flip = f"{pre_itr1}{itr1_flop}{inter_itr_region}{itr2_flip}{post_itr2}"
    flop_flop = f"{pre_itr1}{itr1_flop}{inter_itr_region}{itr2_flop}{post_itr2}"

    # Check that length of orientation sequences is same as original
    for orientation in [flip_flip, flip_flop, flop_flip, flop_flop]:
        if len(orientation) != original_length:
            logger.error(
                f"Incorrect length of orientation plasmid\n"
                f"Orientation length: {len(orientation)}, "
                f"Original length: {original_length}")
            raise ValueError(
                "Error: Issue with length of sequence (orientation vs original)")

    with open(outfile, 'w') as outfile:
        outfile.write(">flip_flip\n")
        outfile.write(f"{flip_flip}\n")
        outfile.write(">flip_flop\n")
        outfile.write(f"{flip_flop}\n")
        outfile.write(">flop_flip\n")
        outfile.write(f"{flop_flip}\n")
        outfile.write(">flop_flop\n")
        outfile.write(f"{flop_flop}\n")


def main(args):
    """Run main entry point."""
    make_flip_flop_aavs(args.fasta, args.bed, args.outfile)
