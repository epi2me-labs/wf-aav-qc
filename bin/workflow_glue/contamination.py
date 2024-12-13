#!/usr/bin/env python
"""Create reference with the various input sequences."""

import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("contam")

    parser.add_argument(
        '--bam_info',
        help="The output of seqkit bam",
        type=Path)
    parser.add_argument(
        '--sample_id',
        help="sample ID")
    parser.add_argument(
        '--transgene_fasta',
        help="The transgene plasmid fasta sequence")
    parser.add_argument(
        '--ref_ids',
        help="Json file with reference:contig mappings")
    parser.add_argument(
        '--n_reads',
        help="Total number of input reads", type=int)
    parser.add_argument(
        '--contam_class_counts', help="Contamination counts output", type=Path)

    return parser


def main(args):
    """Run main entry point."""
    # Extract the reference sequence names from the FATA files,
    # seqkit adds a newline to the output, so remove that with `rstrip`
    transgene_plasmid_name = subprocess.check_output(
        ['seqkit', 'seq', '--name', '--only-id', args.transgene_fasta],
        encoding='UTF-8').rstrip()

    with open(args.ref_ids, 'r') as json_fh:
        ref_ids = json.load(json_fh)

    # Read the per-alignment read summaries
    df_bam = pd.read_csv(
        args.bam_info,
        sep='\t',
        usecols=['Read', 'Ref', 'ReadLen'],
        dtype={
            'Read': str,
            'Ref': str,
            'ReadLen': np.uint32
            }
    )
    # Assign reference category to alignments
    df_bam['contam_class'] = None
    for ref, contigs in ref_ids.items():
        df_bam.loc[df_bam.Ref.isin(contigs), 'contam_class'] = ref
    df_bam.loc[df_bam.Ref == transgene_plasmid_name, 'contam_class'] = 'Transgene'

    # Calculate mapped/unmapped as a percentage of reads not alignments.
    # Note `seqkit bam` does not return info for unammoed reads, so we need to get the
    # total number of input reads separately.
    n_mapped_reads = df_bam.Read.nunique()
    n_input_reads = args.n_reads
    n_unmapped_reads = n_input_reads - n_mapped_reads
    unmapped_pct = 100 / n_input_reads * n_unmapped_reads
    mapped_pct = 100 - unmapped_pct

    df_contam_class = pd.DataFrame(df_bam['contam_class'].value_counts())
    df_contam_class = df_contam_class.rename(columns={'count': 'Number of alignments'})
    df_contam_class['Percentage of alignments'] = (
            100 / len(df_bam) * df_contam_class['Number of alignments']).round(2).T
    df_contam_class.index.name = 'Reference'
    df_contam_class.loc['Unmapped'] = [n_unmapped_reads, unmapped_pct]
    df_contam_class.loc['Mapped'] = [n_mapped_reads, mapped_pct]
    df_contam_class['sample_id'] = args.sample_id

    df_lengths = df_bam[['ReadLen', 'contam_class']]
    df_lengths['sample_id'] = args.sample_id

    df_contam_class.to_csv(args.contam_class_counts, sep='\t')
