"""
Identify truncation hotspots.

Create a CSV file containing start and end locations for alignments fully contained
within the ITR-ITR region.
"""

from pathlib import Path

import pandas as pd

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("truncations")

    parser.add_argument(
        '--bam_info',
        help="The output of seqkit bam",
        type=Path)
    parser.add_argument(
        '--itr_range', help="[itr1_start, itr_2_end]",
        nargs='*', type=int)
    parser.add_argument(
        '--transgene_plasmid_name',
        help="Name of transgene plasmid reference")
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
    # Get the ITR locations
    itr1_start_pos, itr2_end_pos = args.itr_range

    loaded_dataframes = []
    # Get all the reads that map to one of the trans plasmids
    with pd.read_csv(
        args.bam_info,
        sep='\t',
        usecols=['Read', 'Ref', 'Pos', 'EndPos'],
        chunksize=50000
    ) as reader:
        for df_bam in reader:
            df_bam = df_bam.loc[df_bam.Ref.isin([args.transgene_plasmid_name])]

            # Filter for alignments that start and end within the ITR-ITR region
            df_bam = df_bam.loc[
                (df_bam.Pos > itr1_start_pos - 3) &
                (df_bam.EndPos < itr2_end_pos + 3)
            ]
            loaded_dataframes.append(df_bam[['Pos', 'EndPos']])

    # Write csv of start and end positions for plotting
    df_bam = (
        pd.concat(loaded_dataframes)
        .rename(columns={'Pos': 'Read Start', 'EndPos': 'Read end'})
    )
    df_bam['sample_id'] = args.sample_id
    df_bam.to_csv(args.outfile, sep='\t', index=False)
