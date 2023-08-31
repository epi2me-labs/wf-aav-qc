#!/usr/bin/env python
"""
Get the number of reads that are assigned to the various possible transgene structures.

Note: In the original snakemake workflowm this was called
plot_AAV_structures_ssAAV_scAAV.r
"""
from enum import Enum
from pathlib import Path

import polars as pl

from .aav_util import trans_plasmid_ref_names  # noqa: ABS101
from .util import wf_parser  # noqa: ABS101


class GenomeType(str, Enum):
    """Enum for GenomeType categories.

    Subclassing str allows us to access the values as strings and not do .value all over
    the place.
    StrEnum available in Python 3.11
    """

    # Main types
    full_ss = 'Full_ssAAV'
    full_sc = 'Full_scAAV'
    par3_ss = '3` partial ssAAV'
    par5_ss = '5` partial ssAAV'
    bb_int = 'Backbone integration'
    vbb_int = 'Vector-backbone integration'
    par_sc = 'Partial scAAV'
    par_ss = 'Partial ssAAV'

    # Subtypes. Put in separate ENUM?
    almost = 'Almost complete'
    full5_par3 = 'Full 5` ITR and partial 3` ITR'
    par5_full3 = 'Partial 5` ITR and full 3` ITR'
    par_no_itr = 'Partial - no ITRs'
    par5_par3 = 'Partial 5` ITR and partial 3` ITR'
    full5_par_mid = 'Full 5` ITR and partial mid section'
    par_mid_full3 = 'Partial mid section and full 3` ITR'
    par5_par_mid = 'Partial 5` ITR and partial mid section'
    par_mid_par3 = 'Partial mid section and partial 3` ITR'
    itr5 = '5` ITR'
    itr3 = '3` ITR'
    vec_bb_5 = 'Vector backbone - 5` end'
    vec_bb_3 = 'Vector backbone - 3` end'
    bb = 'Backbone'
    ext_itr = 'Extended ITR-ITR region'
    unknown = 'Unknown'

    symetric = 'Symetric'
    asymetric = 'Asymetric'

    sbg_unresolved = 'SBG (unresolved)'
    sbg3 = '3` SBG'
    sbg3_sym = 'SBG 3 symmetric'
    sbg3_asym = 'SBG 3 asymmetric'
    sbg3_incomp = "3' SBG (incomplete)"
    sbg3_incom_asym = 'SBG 3` incomplete ITR asymmetric'
    sbg3_incom_sym = 'SBG 3` incomplete ITR symmetric'

    sbg5 = "5` SBG"
    sbg5_sym = 'SBG 5` symmetric'
    sbg5_asym = 'SBG 5` symmetric'
    sbg5_incomp = "5` SBG (incomplete)"
    sbg5_incom_asym = 'SBG 5` incomplete ITR asymmetric'
    sbg5_incom_sym = 'SBG 5` incomplete ITR symmetric'

    itr_cat = 'ITR concatemer'
    itr1_cat = 'ITR1 concatemer'
    itr2_cat = 'ITR2 concatemer'
    itr1_2_cat = 'ITR1-ITR2 concatemer'

    bb_contam = 'Backbone contamination'
    complex_ = 'Complex'

    par_icg = "Partial ICG - part of both ITRs"  # rename this var
    icg5 = '5` ICG'
    icg3 = '3` ICG'
    par_icg_incom_itrs = 'Partial ICG - incomplete ITRs'
    par_icg_no_itrs = 'Partial ICG - no ITRs'

    gdm = 'GDM'


# this ENUM is used a lot, so let's give it a smaller variable name
gt = GenomeType


def argparser():
    """Create argument parser."""
    parser = wf_parser("structures")

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
        '--itr_padding',
        help="Padding added to ITR regions for checking inclusion in ITR regions",
        type=int,
        default=10)
    parser.add_argument(
        '--output_plot_data', help="Path to output", type=Path)
    parser.add_argument(
        '--output_per_read', help="Path to output", type=Path)

    return parser


def assign_genome_types_to_alignments(trans_df, itr1s, itr1e, itr2s, itr2e, padding):
    """Assign genome_subtypes and make sumamries.

    Based on the relationship of the start and end of the alignemnt to the
    ITR positions, assign a GenomeType to each.

    :param: trans_df
        polars df containing `seqkit bam output`
    :param itr1s: ITR1 start
    :param: itr1e: ITR1 end
    :param: itr2s: ITR2 start
    :param: itr2e ITR2 end
    :param: padding: Padding for inclusion into ITR region
    """
    # Make the column expressions up-front to make the following more concise
    st = pl.col('Pos')
    en = pl.col('EndPos')
    p = padding

    trans_df = trans_df.with_columns(
        pl.when((st >= itr1s) & (st < itr1s + p) & (en > itr2e - p) & (en <= itr2e))
        .then(gt.almost)
        .when((st >= itr1s + p) & (st <= itr1e) & (en > itr2e - p) & (en <= itr2e))
        .then(gt.par5_full3)
        .when((st >= itr1s + p) & (st <= itr1e) & (en >= itr2s) & (en <= itr2e - p))
        .then(gt.par5_par3)
        .when((st >= itr1s) & (st < itr1s + p) & (en >= itr2s) & (en <= itr2e - p))
        .then(gt.full5_par3)
        .when((st >= itr1e) & (en <= itr2s))  #
        .then(gt.par_no_itr)
        .when((st >= itr1s) & (st < itr1s + p) & (en > itr1e) & (en < itr2s))
        .then(gt.full5_par_mid)
        .when((st >= itr1e) & (st <= itr2s) & (en >= itr2e - p) & (en <= itr2e))
        .then(gt.par_mid_full3)
        .when((st >= itr1s + p) & (st < itr1e) & (en > itr1e) & (en < itr2s))
        .then(gt.par5_par_mid)
        .when((st > itr1e) & (st < itr2s) & (en > itr2s) & (en <= itr2e - p))  #
        .then(gt.par_mid_par3)
        .when((st >= itr1s) & (en <= itr1e))
        .then(gt.itr5)
        .when((st >= itr2s) & (en <= itr2e))
        .then(gt.itr3)
        # Check these 2
        .when((st < itr1s) & (en > itr2e))
        .then(gt.ext_itr)
        .when((st < itr1s) & (en >= itr1s))
        .then(gt.vec_bb_5)
        .when((st <= itr2e) & (en > itr2e))
        .then(gt.vec_bb_3)
        .when((st < itr1s) & (en < itr1s))
        .then(gt.bb)
        .when((st > itr2e) & (en > itr2e))
        .then(gt.bb)
        .alias('Genome_type')
    )

    return trans_df


def annotate_reads(trans_df):
    """Assign each read an AAV genome type annotation.

    Use the GenomeTypes applied to the alignemnts for each read to assign
    coarse `Assigned_genome_type` used for summarising in a plot, and more granular
    `Assigned_geome_subtype` for writing to CSV
    """
    # We are accessing a lot of columns, so shorten pl.col
    c = pl.col

    # Get aggregate information for each read from its alignemnts
    # Can we rename Genome_type to alignment_class?
    types_df = trans_df.groupby('Read').agg(
        (c('Ref').unique().str.concat(',')).alias('Reference'),
        # (c('ReadLen').first()).alias('Read_length'),
        # (c('IsSup').sum()).alias('Sup_reads'),
        ((c('Strand').unique().count()-1).cast(pl.Boolean).alias('both_strands')),
        (pl.count().alias('n_alignments')),
        (c('Genome_type').unique().count().alias('n_genome_types')),
        (c('Genome_type').str.concat(',')).alias('Assigned_alignments'),
        (c('Genome_type').count().alias('n_Assigned_alignments')),
    ).sort(by='Read')

    # Apply per-read info back onto the per alignment DF
    # Note this might not be necesarry soon
    trans_df = trans_df.join(
        on='Read',
        other=types_df.select(
            c(['Read', 'both_strands', 'n_genome_types', 'n_alignments'])))

    symmetric_threshold = 10
    # Filter the transdf for alignemnts that may be GDMs then assign reads

    # For a read to be classed as GDM, it must satisfy any of the following conditions
    # Alignemnts are grouoped by Read then if each alignemnt has a match for in
    # the Genome_types specified, set the column 'gdm' to True
    gdm_df = (
        trans_df
        # GDM (genome deletion mutant) -
        # requires 2 alignments both aligning to the same strand
        .filter(c('n_alignments') == 2)
        .filter(~c('both_strands'))
    )
    if not gdm_df.is_empty():
        # Could this be swithched to just check for par_mid Genome_types?
        gdm_df = (
            gdm_df.groupby('Read')
            .agg((
                pl.lit([gt.itr5, gt.par_mid_full3]).is_in(c('Genome_type')).all() |
                pl.lit([gt.itr5, gt.par_mid_par3]).is_in(c('Genome_type')).all() |
                pl.lit([gt.full5_par_mid, gt.itr3]).is_in(c('Genome_type')).all() |
                pl.lit([gt.par5_par_mid, gt.itr3]).is_in(c('Genome_type')).all() |
                pl.lit([gt.full5_par_mid, gt.par_mid_full3]).is_in(
                    c('Genome_type')).all() |
                pl.lit([gt.full5_par_mid, gt.par_mid_par3]).is_in(
                    c('Genome_type')).all() |
                pl.lit([gt.par5_par_mid, gt.par_mid_full3]).is_in(
                    c('Genome_type')).all() |
                pl.lit([gt.par5_par_mid, gt.par_mid_par3]).is_in(
                    c('Genome_type')).all()
            ).alias('gdm'))
        )
        gdm_reads = gdm_df.filter(c('gdm')).select(c('Read')).to_series()
    else:
        gdm_reads = pl.Series([])

    sc_df = (
        # identify scAAV reads and assign subtypes
        trans_df
        # scAAV (self-complementary AAV) - requires 2 alignments aligning to
        # different strands
        .filter(c('n_alignments') == 2)
        .filter(c('both_strands'))
        .filter(c('n_genome_types') == 1)
        .groupby('Read')
        .agg((
            pl.when(
                (c('Genome_type').first() == gt.almost) |
                (c('Genome_type').first() == gt.par5_full3) |
                (c('Genome_type').first() == gt.full5_par3) |
                (c('Genome_type').first() == gt.par5_par3))
            .then(gt.full_sc)
            .when((c('Genome_type').first() == gt.par_mid_full3) &
                  ((c('Pos').take(0) - c('Pos').take(1)).abs()
                  < symmetric_threshold))
            .then(gt.sbg3_sym)

            .when((c('Genome_type').first() == gt.par_mid_full3) &
                  ((c('Pos').take(0) - c('Pos').take(1)).abs()
                   >= symmetric_threshold))
            .then(gt.sbg3_asym)

            .when((c('Genome_type').first() == gt.full5_par_mid) &
                  ((c('EndPos').take(0) - c('EndPos').take(1)).abs()
                  < symmetric_threshold))
            .then(gt.sbg5_sym)

            .when((c('Genome_type').first() == gt.full5_par_mid) &
                  ((c('EndPos').take(0) - c('EndPos').take(1)).abs()
                  >= symmetric_threshold))
            .then(gt.sbg5_asym)

            .when(c('Genome_type').first() == gt.itr5)
            .then(gt.itr1_cat)

            .when(c('Genome_type').first() == gt.itr3)
            .then(gt.itr2_cat)
            .alias('subtype')
        ))
    )

    sc_df2 = (
        trans_df
        # scAAV (self-complementary AAV) - requires 2 alignments aligning to
        # different strands
        .filter(c('n_alignments') == 2)
        .filter(c('both_strands'))
        .filter(c('n_genome_types') == 2)
        .groupby('Read')
        .agg(
            (pl.when(
                (c('Genome_type').is_in([gt.almost, gt.par5_full3]).all()) |
                (c('Genome_type').is_in(
                    [gt.almost, gt.full5_par3]).all()) |
                (c('Genome_type').is_in(
                    [gt.almost, gt.par5_par3]).all()) |
                (c('Genome_type').is_in([
                    gt.full5_par3, gt.par5_par3]).all()) |
                (c('Genome_type').is_in([
                    gt.par5_full3, gt.par5_par3
                ]).all()) |
                (c('Genome_type').is_in([
                    gt.full5_par3, gt.par5_full3
                ]).all())
            )
             .then(gt.full_sc)

             # 3' SBG
             .when(
                (c('Genome_type').is_in([
                    gt.almost, gt.par_mid_full3
                ]).all()) |
                (c('Genome_type').is_in([
                    gt.par_mid_full3, gt.par5_full3
                ]).all())
            )
             .then(gt.sbg3_asym)

             # Incomplete 3' SBG
             .when((c('Genome_type').is_in([gt.almost, gt.par_mid_par3]).all()) |
                   (c('Genome_type').is_in([gt.full5_par3, gt.par_mid_full3]).all()) |
                   (c('Genome_type').is_in([gt.par_mid_par3, gt.par5_full3]).all()) |
                   (c('Genome_type').is_in([gt.par_mid_full3, gt.par5_par3]).all()) |
                   (c('Genome_type').is_in([gt.par_mid_full3, gt.par_mid_par3]).all() &
                    ((c('Pos').take(0) - c('Pos').take(1)).abs()
                        >= symmetric_threshold)) |
                   (c('Genome_type').is_in([gt.full5_par3, gt.par_mid_par3]).all()) |
                   (c('Genome_type').is_in([gt.par_mid_par3, gt.par5_par3]).all()))
             .then(gt.sbg3_incom_asym)
             .when(
                (c('Genome_type').is_in([gt.par_mid_full3, gt.par_mid_par3]).all() &
                 ((c('Pos').take(0) - c('Pos').take(1)).abs()
                  < symmetric_threshold)))
             .then(gt.sbg3_incom_sym)  # Here

             # 5` SBG
             .when(
                (c('Genome_type').is_in([gt.full5_par_mid, gt.full5_par3]).all()) |
                (c('Genome_type').is_in([gt.almost, gt.full5_par_mid]).all()) |
                (c('Genome_type').is_in([gt.almost, gt.par_no_itr]).all()))
             .then(gt.sbg5_asym)

             # Incomplete 5' SBG.
             # SBG_5_incomplete_ITR_asymmetric
             .when(
                (c('Genome_type').is_in([gt.almost, gt.par5_par_mid]).all()) |
                (c('Genome_type').is_in([gt.par5_full3, gt.full5_par_mid]).all()) |
                ((c('Genome_type').is_in(
                    [gt.full5_par_mid, gt.par5_par_mid]).all())
                 & (c('EndPos').take(0) - c('EndPos').take(1)
                    >= symmetric_threshold)) |
                (c('Genome_type').is_in([gt.full5_par_mid, gt.par5_par3]).all()) |
                (c('Genome_type').is_in([gt.par5_par_mid, gt.full5_par3]).all()) |
                (c('Genome_type').is_in([gt.par5_par_mid, gt.par5_par3]).all()) |
                (c('Genome_type').is_in([gt.par5_full3, gt.par5_par_mid]).all())
            )
             .then(gt.sbg5_incom_asym)

             # SBG_5_incomplete_ITR_symmetric
             .when(
                ((c('Genome_type').is_in([gt.full5_par_mid, gt.par5_par_mid])
                  .all())
                 & (c('EndPos').take(0) - c('EndPos').take(1)
                    < symmetric_threshold))
            )
             .then(gt.sbg5_incom_sym)

             # SBG - unresolved
             .when(
                (c('Genome_type').is_in([gt.par5_par3, gt.par_no_itr]).all()) |
                (c('Genome_type').is_in([gt.par5_par_mid, gt.par_mid_par3]).all()) |
                (c('Genome_type').is_in([gt.full5_par_mid, gt.par_mid_par3]).all()) |
                (c('Genome_type').is_in([gt.par5_par_mid, gt.par_mid_full3]).all()) |
                (c('Genome_type').is_in([gt.full5_par_mid, gt.itr3]).all()) |
                (c('Genome_type').is_in([gt.par_mid_full3, gt.itr5]).all()) |
                (c('Genome_type').is_in([gt.par5_par_mid, gt.itr3]).all()) |
                (c('Genome_type').is_in([gt.par_mid_par3, gt.itr5]).all()) |
                (c('Genome_type').is_in([gt.itr5, gt.par_no_itr]).all()) |
                (c('Genome_type').is_in([gt.itr3, gt.par_no_itr]).all()) |
                (c('Genome_type').is_in([gt.full5_par_mid, gt.par_mid_full3]).all()) |
                (c('Genome_type').is_in([gt.almost, gt.itr3]).all()) |
                (c('Genome_type').is_in([gt.par_mid_full3, gt.itr3]).all()) |
                (c('Genome_type').is_in([gt.par5_full3, gt.itr3]).all()) |
                (c('Genome_type').is_in([gt.par_mid_full3, gt.par_no_itr]).all()) |
                (c('Genome_type').is_in([gt.par5_full3, gt.par_no_itr]).all()) |
                (c('Genome_type').is_in([gt.par_mid_par3, gt.par_no_itr]).all()) |
                (c('Genome_type').is_in([gt.par_mid_par3, gt.itr3]).all()) |
                (c('Genome_type').is_in([gt.full5_par3, gt.itr3]).all()) |
                (c('Genome_type').is_in([gt.par5_par3, gt.itr3]).all()) |
                (c('Genome_type').is_in([gt.par5_par3, gt.itr5]).all()) |
                (c('Genome_type').is_in([gt.full5_par3, gt.itr5]).all()) |
                (c('Genome_type').is_in([gt.full5_par_mid, gt.par_no_itr]).all()) |
                (c('Genome_type').is_in([gt.full5_par3, gt.par_no_itr]).all()) |
                (c('Genome_type').is_in([gt.par5_par_mid, gt.par_no_itr]).all()) |
                (c('Genome_type').is_in([gt.par5_par_mid, gt.itr5]).all()) |
                (c('Genome_type').is_in([gt.full5_par_mid, gt.itr5]).all()) |
                (c('Genome_type').is_in([gt.par5_full3, gt.itr5]).all()) |
                (c('Genome_type').is_in([gt.par5_par3, gt.itr5]).all()) |
                (c('Genome_type').is_in([gt.almost, gt.par_no_itr]).all()))
             .then(gt.sbg_unresolved)
             .when((c('Genome_type').is_in([gt.almost, gt.itr5]).all()))
             .then('SBG_5_asymmetric')
             .when((c('Genome_type').is_in([gt.itr5, gt.itr3]).all()))
             .then('ITR1_ITR2_concatemer')

             .when(
                (c('Genome_type') == gt.vec_bb_5).any() |
                (c('Genome_type') == gt.vec_bb_3).any() |
                (c('Genome_type') == gt.bb).any()
            )
             .then('Other')
             .alias('subtype'))  # Assign the result column

        ))

    sc_df = sc_df.vstack(sc_df2)

    # Assign scAAV Assigned_genome_subtypes
    sc_df = sc_df.with_columns(
        pl.when(c('subtype').is_in([gt.sbg5_sym, gt.sbg5_asym]))
        .then(gt.sbg5)
        .when(c('subtype').is_in([gt.sbg3_sym, gt.sbg3_asym]))
        .then(gt.sbg3)
        .when(c('subtype').is_in([gt.sbg5_incom_sym, gt.sbg5_incom_asym]))
        .then(gt.sbg5_incomp)
        .when(c('subtype').is_in([gt.sbg3_incom_sym, gt.sbg3.sbg3_incom_asym]))
        .then(gt.sbg3_incomp)
        .when(c('subtype') == gt.sbg_unresolved)
        .then(gt.sbg_unresolved)
        .when(c('subtype').is_in([gt.itr1_cat, gt.itr2_cat, gt.itr1_2_cat]))
        .then(gt.itr_cat)
        .when(c('subtype') == gt.full_sc)
        .then(gt.full_sc)
        .alias('Assigned_genome_subtype')
    )
    # 3` sbg incomplete fine in sc_df
    types_df = types_df.with_columns(
        pl.when(
            c('Assigned_alignments').is_in([gt.vec_bb_5, gt.vec_bb_3, gt.bb]))
        .then(gt.bb_contam)
        .when(c('n_Assigned_alignments') > 2)
        .then(gt.complex_)
        .alias('Assigned_genome_subtype')
    )

    types_df2 = (types_df.filter(c('n_Assigned_alignments') == 1))
    types_df2 = types_df2.with_columns(
        pl.when(c('Assigned_alignments') == pl.lit(gt.almost))
        .then(gt.full_ss)
        .when(c('Assigned_alignments') == pl.lit(gt.par5_par3))
        .then(pl.lit(gt.par_icg))  # Incomplete genome
        .when((c('Assigned_alignments') == pl.lit(gt.full5_par3)) |
              (c('Assigned_alignments') == pl.lit(gt.full5_par_mid)))
        .then(gt.icg5)
        .when((c('Assigned_alignments') == pl.lit(gt.par5_full3)) |
              (c('Assigned_alignments') == pl.lit(gt.par_mid_full3)))
        .then(gt.icg3)
        .when((c('Assigned_alignments') == pl.lit(gt.par5_par_mid)) |
              (c('Assigned_alignments') == pl.lit(gt.par_mid_par3)))
        .then(gt.par_icg_incom_itrs)
        .when(c('Assigned_alignments') == pl.lit(gt.par_no_itr))
        .then(gt.par_icg_no_itrs)
        .when(c('Read').is_in(gdm_reads))
        .then(gt.gdm)
        .alias('Assigned_genome_subtype')
     )
    # Join the new Asssigned_genotype_type (not working)
    types_df = (types_df.join(
        types_df2.select(['Read', 'Assigned_genome_subtype']),
        on='Read',
        how='outer')
        .with_columns(
        pl.when(c('Assigned_genome_subtype') == None)  # noqa: E711
        .then(c('Assigned_genome_subtype_right'))
        .otherwise(c('Assigned_genome_subtype'))
        .alias('Assigned_genome_subtype')).drop('Assigned_genome_subtype_right')
    )

    # Merge the scAAV df to the types df
    # Some duplication of code here - fix
    types_df = (types_df.join(
        sc_df, on='Read', how='left')
        .with_columns(
            pl.when(c('Assigned_genome_subtype') == None)  # noqa: E711
            .then(c('Assigned_genome_subtype_right'))
            .otherwise(c('Assigned_genome_subtype'))
            .alias('Assigned_genome_subtype')
        )
        .drop('Assigned_genome_subtype_right')
    )
    # Add extra info
    types_df = types_df.with_columns(
        pl.when(c('subtype').is_in([
            gt.sbg3_sym, gt.sbg5_incom_sym, gt.sbg3_incom_sym]))
        .then(gt.symetric)
        .when(c('subtype').is_in([
            gt.sbg3_asym, gt.sbg5_incom_asym, gt.sbg3_incom_asym]))
        .then(gt.asymetric)
        .when(c('subtype') == gt.itr1_cat)
        .then(gt.itr1_cat)
        .when(c('subtype') == gt.itr2_cat)
        .then(gt.itr2_cat)
        .when(c('subtype') == gt.itr1_2_cat)
        .then(gt.itr1_2_cat)
        .alias('Assigned_extra_info')
    )

    # Assign assigned genome types
    # Why are we reassigning subtypes here
    types_df = types_df.with_columns(
        pl.when(c('Assigned_genome_subtype') == gt.bb_contam)
        .then(gt.bb_contam)
        .when(c('Assigned_genome_subtype') == gt.full_ss)
        .then(gt.full_ss)
        .when(c('Assigned_genome_subtype').is_in([
            gt.icg5, gt.icg3, gt.par_icg_incom_itrs, gt.par_icg_no_itrs, gt.gdm]))
        .then(gt.par_ss)
        .when(c('Assigned_genome_subtype') == gt.full_sc)
        .then(gt.full_sc)
        .when(c('Assigned_genome_subtype') == gt.itr_cat)
        .then(gt.itr_cat)
        .when(c('Assigned_genome_subtype') == gt.complex_)
        .then(gt.complex_)
        .otherwise(gt.unknown)
        .alias('Assigned_genome_type')
    )

    types_df = types_df.with_columns(
        pl.when(c('Assigned_genome_subtype') == None)  # noqa: E711
        .then(gt.unknown)
        .otherwise(c('Assigned_genome_subtype'))
        .alias('Assigned_genome_subtype'))

    # create summary
    assigned_genome_types_summary = (
        types_df.groupby('Assigned_genome_type').count()
        .with_columns(
            (c('count') / c('count').sum() * 100).round(2)
            .alias('percentage'))
    )

    per_read_output = types_df.select([
        'Read',
        'Assigned_genome_subtype',
        'Assigned_extra_info'
    ])
    return assigned_genome_types_summary, per_read_output


def main(args):
    """Entry point."""
    c = pl.col
    # Load the annotation BED file
    ann_schema = {
        'ref': str,
        'start': pl.UInt32,
        'end': pl.UInt32,
        'feature': str
    }
    df_ann = pl.read_csv(
        source=args.annotation,
        has_header=False,
        separator='\t',
        dtypes=ann_schema
    )

    # Load the BAM info file
    df_bam = pl.read_csv(
        source=args.bam_info,
        separator='\t'
    )

    # Get rid of secondary alignemnts
    df_bam = df_bam.filter(c('IsSec') == 0)
    df_bam = df_bam.drop(columns=['IsSec'])

    # Get ITR start and stop positions
    itr_bed = df_ann.filter(c('feature').is_in(['ITR1', 'ITR2']))
    if itr_bed.shape[0] != 2:
        raise ValueError('Bed file should contain both ITR1 and ITR2 entries')
    itr1_start = itr_bed.item(0, 'start')
    itr1_end = itr_bed.item(0, 'end')
    itr2_start = itr_bed.item(1, 'start')
    itr2_end = itr_bed.item(1, 'end')

    # Get reads that map to trans plasmid
    non_plasmid_read_ids = (
        df_bam.filter(~c('Ref').is_in(trans_plasmid_ref_names))
        .select(c('Read')).to_series()
    )

    trans_df = df_bam.filter(~c('Read').is_in(non_plasmid_read_ids))

    # Get dataframe containing alignemnt annotatoed with genome type
    df_assigned_alns = assign_genome_types_to_alignments(
        trans_df,
        itr1_start,
        itr1_end,
        itr2_start,
        itr2_end,
        padding=10  # Remove hard code
    )

    data_for_plot, per_read_summary = annotate_reads(
        df_assigned_alns)

    per_read_summary = (
        per_read_summary.with_columns(pl.lit(args.sample_id).alias('sample_id'))
    )
    data_for_plot = (
        data_for_plot.with_columns(pl.lit(args.sample_id).alias('sample_id'))
    )

    data_for_plot.write_csv(args.output_plot_data, separator='\t')
    per_read_summary.write_csv(args.output_per_read, separator='\t')
