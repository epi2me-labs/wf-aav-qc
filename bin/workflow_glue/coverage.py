"""Trim the depth file to positions within the ITR-ITR region."""
from pathlib import Path

import pandas as pd
import polars as pl

from .aav_util import load_annotation  # noqa: ABS101
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("coverage")

    parser.add_argument(
        '--annotation',
        help="transgene plasmidf annotation file",
        type=Path)

    parser.add_argument(
        '--itr_coverage',
        help="The output of samtools depth",
        type=Path)

    parser.add_argument(
        '--output',
        help="The output path",
        type=Path)

    return parser


def main(args):
    """Run main entry point."""
    ann_df = load_annotation(args.annotation)
    cov_df = pd.read_csv(args.itr_coverage, sep='\t')

    itr1_start = ann_df.filter(
        pl.col('feature') == 'ITR1').select(pl.col('start')).item()
    itr2_end = ann_df.filter(
        pl.col('feature') == 'ITR2').select(pl.col('end')).item()

    cov_df = cov_df.loc[
        (cov_df['pos'] > itr1_start - 3) & (cov_df['pos'] < itr2_end + 3)]

    cov_df.to_csv(args.output, sep='\t', index=False)
