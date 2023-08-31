"""Create workflow report."""
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "--bed_in", help="Transgene annotation file")
    parser.add_argument(
        "--bed_out", help="Path for modified transgene annotation file")

    return parser


def main(args):
    """Add flaking bases to ITR region in bed file."""
    with open(args.bed_in, 'r') as infile, open(args.bed_out, 'w') as outfile:
        for line in infile:
            cols = line.strip("\n").split("\t")
            if "ITR1" in cols[3]:
                outfile.write(
                    f"{cols[0]}\t{cols[1]}\t{str(int(cols[2]) + 100)}\t{cols[3]}\n")
            elif "ITR2" in cols[3]:
                outfile.write(
                    f"{cols[0]}\t{str(int(cols[1]) - 100)}\t{cols[2]}\t{cols[3]}\n")
            else:
                outfile.write(f"{cols[0]}\t{cols[1]}\t{cols[2]}\t{cols[3]}\n")
