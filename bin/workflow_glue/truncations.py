"""
Identify truncation hotspots.

Create a CSV file with start and end locations for alignments fully contained within
the ITR-ITR region.
"""

from pathlib import Path

import pandas as pd

from .aav_util import trans_plasmid_ref_names  # noqa: ABS101
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("truncations")

    parser.add_argument(
        '--bam_info',
        help="The output of seqkit bam",
        type=Path)
    parser.add_argument(
        '--annotation',
        help="BED file with plasmid annotations",
        type=Path)
    parser.add_argument(
        '--sample_id',
        help="sample ID")
    parser.add_argument(
        '--outfile',
        help="Path to output",
        type=Path)

    return parser


def main(args):
    """Run main entry point."""
    df_ann = pd.read_csv(
        args.annotation, sep='\t',
        names=['serotype', 'start', 'end', 'feature']
    ).set_index('feature')

    df_bam = pd.read_csv(
        args.bam_info, sep='\t', usecols=['Read', 'Ref', 'Pos', 'EndPos'])

    # Get all the reads that map to one of the trans plasmids
    df_bam = df_bam.loc[df_bam.Ref.isin(trans_plasmid_ref_names)]

    # Filter for reads that start and end within the ITR-ITR region
    df_bam = df_bam.loc[
        (df_bam.Pos > df_ann.loc['ITR1', 'start'] - 3) &
        (df_bam.EndPos < df_ann.loc['ITR2', 'end'] + 3)
    ]

    # Write csv of start and end positions for plotting
    df_bam = df_bam[['Pos', 'EndPos']]
    df_bam.rename(columns={'Pos': 'Read Start', 'EndPos': 'Read end'})
    df_bam['sample_id'] = args.sample_id
    df_bam.to_csv(args.outfile, sep='\t', index=False)
