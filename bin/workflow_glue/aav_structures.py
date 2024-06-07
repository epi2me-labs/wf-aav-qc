"""Categorise reads by gene type."""

from pathlib import Path

import polars as pl

from .structure_definitions import (  # noqa: ABS101
    ReadType, AlnType, get_subtype_definitions
)
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("structures")

    parser.add_argument(
        '--bam_info',
        help="File containing the output of seqkit bam",
        type=Path)
    parser.add_argument(
        '--itr_locations', help="[itr1_start, itr1_end, itr2_start, itr_2_end]",
        nargs='*', type=int)
    parser.add_argument(
        '--sample_id',
        help="'sample_id' is added as a column to both output files")
    parser.add_argument(
        '--itr_fl_threshold',
        help="Max missing bases to assign an ITR as full length",
        type=int)
    parser.add_argument(
        '--itr_backbone_threshold',
        help=("Alignments starting this many bases or more before ITR1 or after ITR2 "
              "will be clased as `backbone`"),
        type=int)
    parser.add_argument(
        '--transgene_plasmid_name',
        help="Name of transgene plasmid reference")
    parser.add_argument(
        '--output_plot_data',
        help="Path to output TSV used for plotting",
        type=Path)
    parser.add_argument(
        '--symmetry_threshold',
        help=("Threshold to determine whether alignments on opposite strand are "
              "symmetric or asymmetric."),
        type=int)
    parser.add_argument(
        '--output_per_read',
        help="Path to output TSV for inclusion in output directory",
        type=Path)

    return parser


def assign_genome_types_to_alignments(
        trans_df, itr1_start, itr1_end, itr2_start, itr2_end,
        itr_fl_threshold, itr_backbone_threshold):
    """Assign a type to each alignment.

    Assign an AlnType based on the positions of the start and end of the alignment
    relative to the ITRs.

    :param: trans_df
        polars df containing `seqkit bam` output of reads mapping to transgene plasmid
    :param: itr1_start: ITR1 start position
    :param: itr1_end: ITR1 end position
    :param: itr2_start: ITR2 start position
    :param: itr2_end ITR2 end position
    :param: itr_fl_threshold: Max allowed missing bases from ITR to class as full length
    :param: itr_backbone_threshold: Max allow extra bases at ITR ends to not be classed
        as backbone
    :return: polars Dataframe. Input dataframe with additional al_type column
    """
    # Make these polars column expressions up-front to make the following more concise
    start = pl.col('Pos')
    end = pl.col('EndPos')

    itr1_full = itr1_start + itr_fl_threshold
    itr2_full = itr2_end - itr_fl_threshold

    # Let's build some expressions up-front as they're used multiple times
    full_5prime = (start >= itr1_start - itr_backbone_threshold) & (start < itr1_full)
    full_3prime = (end > itr2_full) & (end <= itr2_end + itr_backbone_threshold)

    partial_5prime = (start >= itr1_full) & (start <= itr1_end)
    partial_3prime = (end >= itr2_start) & (end <= itr2_full)

    start_in_midsection = (start > itr1_end) & (start <= itr2_start)
    end_in_midsection = (end > itr1_end) & (end < itr2_start)

    trans_df = trans_df.with_columns(
        pl.when(full_5prime & full_3prime).then(AlnType.almost)

        .when(partial_5prime & full_3prime).then(AlnType.par5_full3)

        .when(partial_5prime & partial_3prime).then(AlnType.par5_par3)

        .when(full_5prime & partial_3prime).then(AlnType.full5_par3)

        .when((start > itr1_end) & (end <= itr2_start)).then(AlnType.par_no_itr)

        .when(full_5prime & end_in_midsection).then(AlnType.full5_par_mid)

        .when(start_in_midsection & full_3prime).then(AlnType.par_mid_full3)

        .when(partial_5prime & end_in_midsection).then(AlnType.par5_par_mid)

        .when(start_in_midsection & partial_3prime).then(AlnType.par_mid_par3)

        .when(
            (start >= itr1_start - itr_backbone_threshold) &
            (end <= itr1_end)).then(AlnType.itr5_only)

        .when(
            (start >= itr2_start) &
            (end <= itr2_end + itr_backbone_threshold)).then(AlnType.itr3_only)

        .when(
            (start < itr1_start - itr_backbone_threshold) &
            (end > itr2_end + itr_backbone_threshold)).then(AlnType.ext_itr)

        # Transgene plasmid backbone alignments
        .when(
            (start < itr1_start - itr_backbone_threshold) &
            (end >= itr1_start)).then(AlnType.vec_bb_5)

        .when(
            (start <= itr2_end) &
            (end > itr2_end + itr_backbone_threshold)).then(AlnType.vec_bb_3)

        .when((start < itr1_start) & (end < itr1_start)).then(AlnType.bb)

        .when((start > itr2_end) & (end > itr2_end)).then(AlnType.bb)

        # Unknown
        .otherwise(AlnType.unknown)

        .alias('aln_type')
    )

    return trans_df


def annotate_reads(aln_df, type_definitions, symmetry_threshold):
    """Assign each read an AAV ReadType annotation.

    Each read assigned to both of these categories:
    `Assigned_genome_type`, a coarse category used for in summary plot
    `Assigned_geome_subtype`, a more granular category for writing to TSV

    :param: aln_df
        polars df containing `seqkit bam output` of reads mapping to transgene plasmid
        with added column of `aln_type`, one row per alignment.
    :param: type_definitions: List[ReadTypeDefinition]
    :param: symmetry_threshold: min distance between start or end positions on opposite
        strands when testing for symmetry.
    :return: polars.DataFrame with columns: [Read, Assigned_genome_subtype, extra_info]
    """
    # Aggregate Alignments on Read. Apply some read summary info back to the alignment
    # DataFrame
    aln_df = aln_df.join(
        other=aln_df.groupby('Read')
        .agg(
            (pl.col('Strand').unique().count().alias('n_strands')),
            (pl.count('aln_type').alias('n_alignments')),
            (pl.col('aln_type').unique().count().alias('n_aln_types'))
        ).select(pl.col(['Read', 'n_strands', 'n_aln_types', 'n_alignments'])),
        on='Read',
        how='left'
    )

    # Calculate symmetry status at the left and right positions
    # where there are two alignments on both strands
    sym_results = (
        aln_df.filter(
            (pl.col('n_strands') == 2) &
            (pl.col('n_alignments') == 2)
        )
        .groupby('Read')
        .agg(
            ((pl.col('Pos').take(0) - pl.col('Pos').take(1)).abs()
             <= symmetry_threshold).alias('5prime_sym'),
            ((pl.col('EndPos').take(0) - pl.col('EndPos').take(1)).abs()
             <= symmetry_threshold).alias('3prime_sym')
        )
    )
    aln_df = aln_df.join(
        on='Read',
        other=sym_results,
        how='left'
    )

    read_subtype_dfs = []  # Results

    for subtype_def in type_definitions:
        # Add filter expressions associated with each definition to a list
        filter_expressions = []
        if subtype_def.n_alignments != 0:
            filter_expressions.append(
                pl.col('n_alignments') == subtype_def.n_alignments)
        if subtype_def.n_strands != 0:
            filter_expressions.append(
                pl.col('n_strands') == subtype_def.n_strands)
        if subtype_def.n_aln_types != 0:
            filter_expressions.append(
                pl.col('n_aln_types') == subtype_def.n_aln_types)
        # Subset the alignments based on the filtering criteria
        df_subtype = aln_df.filter(pl.all_horizontal(filter_expressions))

        if df_subtype.is_empty():  # All alignments filtered
            continue

        if subtype_def.symmetry:
            # Apply any required symmetry filter
            if subtype_def.symmetry == 'end_sym':
                df_subtype = df_subtype.filter(
                    pl.col('3prime_sym'))
            elif subtype_def.symmetry == 'end_asym':
                df_subtype = df_subtype.filter(
                    ~pl.col('3prime_sym'))
            elif subtype_def.symmetry == 'start_sym':
                df_subtype = df_subtype.filter(
                    pl.col('5prime_sym'))
            elif subtype_def.symmetry == 'start_asym':
                df_subtype = df_subtype.filter(
                    ~pl.col('5prime_sym'))

        # Generate polars expressions with filtering on AlnTypes
        aln_type_expressions = []
        for condition in subtype_def.aln_type_combinations:

            if condition == 'no_overlap':
                aln_type_expressions.append(
                    (pl.col('Pos').max() - pl.col('EndPos').min() > 0))

            if condition == 'overlap':
                aln_type_expressions.append(
                    (pl.col('Pos').max() - pl.col('EndPos').min() <= 0))

            elif len(condition) == 1:
                # Only one aln_type defines this ReadType. Check if it is present.
                exp = pl.any_horizontal(
                    pl.lit(condition[0].value).is_in(pl.col("aln_type")))
                aln_type_expressions.append(exp)

            elif len(condition) == 2:
                # Two aln_types define this condition. Check if both present.
                aln_type_expressions.append(
                    pl.lit(condition[0]).is_in(pl.col("aln_type")) &
                    pl.lit(condition[1]).is_in(pl.col("aln_type"))

                )
        # Apply the aln_type expressions
        assigned_reads_df = (
            df_subtype.groupby('Read')
            .agg(
                pl.any_horizontal(aln_type_expressions)
                .alias('Assigned_genome_subtype'))
            # Keep only reads that have been assigned (True)
            .filter(pl.col('Assigned_genome_subtype'))
            # Set the name of the Assigned_genome_type to the current definition name
            .with_columns(pl.lit(subtype_def.type_).alias('Assigned_genome_subtype'))
        )

        # Make extra_info columns. Only required for a subset of types
        extra_info = '-'
        if subtype_def.symmetry:
            if subtype_def.symmetry in ['start_sym', 'end_sym']:
                extra_info = 'symmetrical'
            if subtype_def.symmetry in ['start_asym', 'end_asym']:
                extra_info = 'asymmetrical'
        elif subtype_def.type_ == ReadType.itr1_cat:
            extra_info = ReadType.itr1_cat
        elif subtype_def.type_ == ReadType.itr2_cat:
            extra_info = ReadType.itr2_cat
        elif subtype_def.type_ == ReadType.itr1_2_cat:
            extra_info = ReadType.itr1_2_cat
        assigned_reads_df = assigned_reads_df.with_columns(
            pl.lit(extra_info).alias('extra_info'))

        read_subtype_dfs.append(assigned_reads_df)

    # Assign complex (> 2 alignments)
    complex_ = (
        aln_df.filter(pl.col('n_alignments') > 2)
        .with_columns(pl.lit('Complex').alias('Assigned_genome_subtype'))
        .select(['Read', 'Assigned_genome_subtype'])
    ).unique()
    read_subtype_dfs.append(complex_)
    all_reads_assigned_df = pl.concat(read_subtype_dfs, how="diagonal")

    # Assign 'Unknown'
    unknown_df = (
        (
          aln_df.filter(
              pl.col('Read').is_in(all_reads_assigned_df['Read']).is_not())[
              ['Read']]
        ).unique()
        .with_columns(
            pl.lit(ReadType.unknown)
            .alias('Assigned_genome_subtype')
         )
        .select(['Read', 'Assigned_genome_subtype'])
    )
    merged_df = pl.concat([all_reads_assigned_df, unknown_df], how='diagonal')

    # Backbone fix.
    # It's possible to get duplicate read assignments with backbone and another type
    # For example complex (with 2 strands with overlap) can also have backbone in it
    # If any read is assigned backbone, remove any other reads assigned to other
    # categories
    dups = merged_df.filter(pl.col("Read").is_duplicated())

    # Get read IDs that have backbone assignment
    bb_ids = dups.filter(
        pl.col('Assigned_genome_subtype').is_in(pl.lit(ReadType.bb_contam)))['Read']

    # Drop reads that also have the same read id but not backbone
    dup_df = dups.filter(
        (pl.col('Read').is_in(bb_ids))
        & (pl.col('Assigned_genome_subtype') != pl.lit(ReadType.bb_contam)))

    final_df = merged_df.join(
        other=dup_df[['Read', 'Assigned_genome_subtype']],
        on=['Read', 'Assigned_genome_subtype'],
        how="anti"
    )
    # f = final_df.to_pandas()
    # f1 = f[f.Read.duplicated(keep=False)].sort_values('Read')
    return final_df


def main(args):
    """Entry point."""
    # Load the BAM info file
    schema = {
        'Ref':  pl.Utf8,
        'Read': pl.Utf8,
        'Pos': pl.UInt32,
        'EndPos': pl.UInt32,
        'ReadLen': pl.UInt32,
        'Strand': pl.UInt8,
        'IsSec': pl.UInt8,
        'IsSup': pl.UInt8
    }

    df_bam = (pl.read_csv(
        source=args.bam_info,
        separator='\t',
        columns=list(schema.keys()),
        dtypes=list(schema.values())
        )
        .with_columns([
            pl.col('Strand').cast(pl.Boolean),
            pl.col('IsSec').cast(pl.Boolean),
            pl.col('IsSup').cast(pl.Boolean)
        ])
    )

    # Get the ITR locations
    itr1_start, itr1_end, itr2_start, itr2_end = args.itr_locations

    # Get reads that map to transgene plasmid
    non_plasmid_read_ids = (
        df_bam.filter(~pl.col('Ref').is_in([args.transgene_plasmid_name]))
        .select(pl.col('Read')).to_series()
    )
    trans_df = df_bam.filter(~pl.col('Read').is_in(non_plasmid_read_ids))

    # Assign types to each alignment
    df_assigned_alns = assign_genome_types_to_alignments(
        trans_df,
        itr1_start,
        itr1_end,
        itr2_start,
        itr2_end,
        args.itr_fl_threshold,
        args.itr_backbone_threshold
    )

    # Assign final types to each read
    df_assigned_reads = annotate_reads(
        df_assigned_alns, get_subtype_definitions(), args.symmetry_threshold)

    # Map SBG with symmetric subtypes up to SBG
    # This is here as the tests are expecting the SBG subcategories.
    df_assigned_reads = df_assigned_reads.with_columns(
        pl.when(pl.col('Assigned_genome_subtype').is_in(
            [ReadType.sbg5_sym, ReadType.sbg5_asym]))
        .then(ReadType.sbg5)
        .when(pl.col('Assigned_genome_subtype').is_in(
            [ReadType.sbg3_sym, ReadType.sbg3_asym]))
        .then(ReadType.sbg3)
        .when(pl.col('Assigned_genome_subtype').is_in(
            [ReadType.sbg5_incom_sym, ReadType.sbg5_incom_asym]))
        .then(ReadType.sbg5_incomp)
        .when(pl.col('Assigned_genome_subtype').is_in(
            [ReadType.sbg3_incom_sym, ReadType.sbg3.sbg3_incom_asym]))
        .then(ReadType.sbg3_incomp)
        .otherwise(pl.col('Assigned_genome_subtype'))
        .alias('Assigned_genome_subtype')
    )

    # Add sample_id column
    per_read_summary = (
        df_assigned_reads.with_columns(pl.lit(args.sample_id).alias('sample_id'))
    )
    # Write the per read outputs with an assigned genome subtype.
    per_read_summary.write_csv(args.output_per_read, separator='\t')
    # Get the more coarse Assigned_genome_type for the plots.
    per_read_summary = per_read_summary.with_columns(
        pl.when(pl.col('Assigned_genome_subtype') == ReadType.bb_contam)
        .then(ReadType.bb_contam)
        .when(pl.col('Assigned_genome_subtype') == ReadType.full_ss)
        .then(ReadType.full_ss)
        .when(pl.col('Assigned_genome_subtype').is_in([
            ReadType.icg5, ReadType.icg3, ReadType.par_icg_incom_itrs,
            ReadType.par_icg_no_itrs, ReadType.gdm]))
        .then(ReadType.par_ss)
        .when(pl.col('Assigned_genome_subtype') == ReadType.full_sc)
        .then(ReadType.full_sc)
        .when(pl.col('Assigned_genome_subtype').is_in([
            ReadType.itr1_cat, ReadType.itr2_cat, ReadType.itr1_2_cat,
            ReadType.itr3_only, ReadType.itr5_only, ReadType.single_itr]))
        .then(ReadType.itr_only)
        .when(pl.col('Assigned_genome_subtype') == ReadType.complex_)
        .then(ReadType.complex_)
        .when(pl.col('Assigned_genome_subtype').is_in([
            ReadType.sbg5, ReadType.sbg3, ReadType.sbg5_incomp,
            ReadType.sbg3_incomp, ReadType.sbg_unresolved]))
        .then(ReadType.par_sc)
        .otherwise(ReadType.unknown)
        .alias('Assigned_genome_type')
    )

    # create summary
    assigned_genome_types_summary = (
        per_read_summary.groupby('Assigned_genome_type').count()
        .with_columns(
            (pl.col('count') / pl.col('count').sum() * 100).round(2)
            .alias('percentage'))
    )
    data_for_plot = (
        assigned_genome_types_summary
        .with_columns(pl.lit(args.sample_id).alias('sample_id'))
    )
    # Write the summary
    data_for_plot.write_csv(args.output_plot_data, separator='\t')
