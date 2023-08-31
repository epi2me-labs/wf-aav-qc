#!/usr/bin/env python
"""Create reference with the various flip-flop orientations."""

from pathlib import Path
import subprocess

import pandas as pd

from .aav_util import trans_plasmid_ref_names   # noqa: ABS101
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("truncations")

    parser.add_argument(
        '--bam_info',
        help="The output of seqkit bam",
        type=Path)
    parser.add_argument(
        '--sample_id',
        help="sample ID")
    parser.add_argument(
        '--helper_fasta',
        help="The helper fast sequence ID")
    parser.add_argument(
        '--rep_cap_fasta',
        help="The rep_cap fasta sequence ID")
    parser.add_argument(
        '--read_ids',
        help="TSV with containing all read IDs")
    parser.add_argument(
        '--contam_class_counts', help="Contamination counts output", type=Path)
    parser.add_argument(
        '--contam_read_lengths', help="Contamination lengths output", type=Path)

    return parser


def main(args):
    """Run main entry point."""
    helper_name = subprocess.check_output(
        ['sed', '-n', 's/^>//p', f"{args.helper_fasta}"], encoding='UTF-8').rstrip()

    rep_cap_name = subprocess.check_output(
        ['sed', '-n', 's/^>//p', f"{args.rep_cap_fasta}"], encoding='UTF-8').rstrip()

    df_bam = pd.read_csv(args.bam_info, sep='\t')
    # Assign reference category to reads
    df_bam['contam_class'] = None
    df_bam.loc[df_bam.Ref.isin(trans_plasmid_ref_names), 'contam_class'] = 'Transgene'
    df_bam.loc[df_bam.Ref.str.startswith('chr'), 'contam_class'] = 'Human'
    df_bam.loc[df_bam.Ref == helper_name, 'contam_class'] = 'Helper'
    df_bam.loc[df_bam.Ref == rep_cap_name, 'contam_class'] = 'Rep-cap'

    all_read_ids = pd.read_csv(args.read_ids, header=None)[0]
    unmapped = len(all_read_ids[~all_read_ids.isin(df_bam.Read)])
    unmapped_pct = 100 / len(df_bam) * unmapped
    df_contam_class = pd.DataFrame(df_bam['contam_class'].value_counts())
    df_contam_class = df_contam_class.rename(columns={'count': 'Number of alignments'})
    df_contam_class['Percentage of alignments'] = (
            100 / len(df_bam) * df_contam_class['Number of alignments']).round(2).T
    df_contam_class.index.name = 'Reference'
    df_contam_class.loc['Unmapped'] = [unmapped, unmapped_pct]
    df_contam_class['sample_id'] = args.sample_id

    df_lengths = df_bam[['ReadLen', 'contam_class']]
    df_lengths['sample_id'] = args.sample_id

    df_contam_class.to_csv(args.contam_class_counts, sep='\t')
    df_lengths.to_csv(args.contam_read_lengths, sep='\t', index=False)
