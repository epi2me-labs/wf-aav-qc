"""Categorise reads by gene type.

The script classifies per-read AAV genome structures from alignment evidence on the
transgene plasmid. Conceptually it runs in three layers:

  1.	Per-alignment labeling → assign an AlnType to each alignment relative to
  ITRs and backbone.
  2.	Per-read rule evaluation → apply ordered rules (definitions) that combine
  alignment patterns (and symmetry/overlap constraints) into read subtypes;
  then map those subtypes to coarse genome types for reporting.
  3.	Outputs & summaries → produce a per-read table (subtype + coarse type) and
  an aggregate count/percent breakdown.
"""

import contextlib
from pathlib import Path

import polars as pl
import pysam

from .structure_definitions import (  # noqa: ABS101
    ReadType, AlnType, get_subtype_definitions
)
from .util import wf_parser, get_named_logger  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("structures")

    parser.add_argument(
        '--bam_in',
        help="Path to input BAM file to tag with assigned genome type")
    parser.add_argument(
        '--bam_out',
        help="Out path for tagged BAM")
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
              "will be classed as `backbone`"),
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


def _l(x):
    """Convert string to polars literal."""
    return pl.lit(x.value)


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
        pl.when(full_5prime & full_3prime)
        .then(_l(AlnType.almost))

        .when(partial_5prime & full_3prime)
        .then(_l(AlnType.par5_full3))

        .when(partial_5prime & partial_3prime)
        .then(_l(AlnType.par5_par3))

        .when(full_5prime & partial_3prime)
        .then(_l(AlnType.full5_par3))

        .when((start > itr1_end) & (end <= itr2_start))
        .then(_l(AlnType.par_no_itr))

        .when(full_5prime & end_in_midsection)
        .then(_l(AlnType.full5_par_mid))

        .when(start_in_midsection & full_3prime)
        .then(_l(AlnType.par_mid_full3))

        .when(partial_5prime & end_in_midsection)
        .then(_l(AlnType.par5_par_mid))

        .when(start_in_midsection & partial_3prime)
        .then(_l(AlnType.par_mid_par3))

        .when(
            (start >= itr1_start - itr_backbone_threshold) &
            (end <= itr1_end))
        .then(_l(AlnType.itr5_only))

        .when(
            (start >= itr2_start) &
            (end <= itr2_end + itr_backbone_threshold))
        .then(_l(AlnType.itr3_only))

        .when(
            (start < itr1_start - itr_backbone_threshold) &
            (end > itr2_end + itr_backbone_threshold))
        .then(_l(AlnType.ext_itr))

        # Transgene plasmid backbone alignments
        .when(
            (start < itr1_start - itr_backbone_threshold) &
            (end >= itr1_start))
        .then(_l(AlnType.vec_bb_5))

        .when(
            (start <= itr2_end) &
            (end > itr2_end + itr_backbone_threshold))
        .then(_l(AlnType.vec_bb_3))

        # Backbone only, no transgene cassette sequence
        .when((start < itr1_start) & (end < itr1_start))
        .then(_l(AlnType.bb))

        .when((start > itr2_end) & (end > itr2_end))
        .then(_l(AlnType.bb))

        # Unknown
        .otherwise(_l(AlnType.transgene_unclassified))

        .alias('aln_type')
    )

    return trans_df


def annotate_reads(sample_id, aln_df, type_definitions, symmetry_threshold):
    """Assign each read an AAV ReadType annotation.

    Each read assigned to both of these categories:
    `Assigned_genome_type`, a coarse category used for in summary plot
    `Assigned_geome_subtype`, a more granular category for writing to TSV

    :param: aln_df
        polars DF containing `seqkit bam output` of reads mapping to transgene plasmid
        with added column of `aln_type`, one row per alignment.
    :param: type_definitions: List[ReadTypeDefinition]
    :param: symmetry_threshold: min distance between start or end positions on opposite
        strands when testing for symmetry.
    :return: polars.DataFrame with columns: [Read, Assigned_genome_subtype, extra_info]
    """
    # Annotate each alignment row with per-read summary statistics.
    # Note that using `.over` expressions to calculate these stats is more efficient
    # than a groupby-agg-join approach.
    aln_df = aln_df.with_columns([
        pl.col('Strand').n_unique().over('Read').alias('n_strands'),
        pl.len().over('Read').alias('n_alignments'),
        pl.col('aln_type').n_unique().over('Read').alias('n_aln_types'),
    ])

    # Calculate symmetry status at the left and right positions
    # where there are two alignments on both strands
    sym_results = (
        aln_df.filter(
            (pl.col('n_strands') == 2) &
            (pl.col('n_alignments') == 2)
        )
        .group_by('Read')
        .agg(
            ((pl.col('Pos').get(0) - pl.col('Pos').get(1)).abs()
             <= symmetry_threshold).alias('5prime_sym'),
            ((pl.col('EndPos').get(0) - pl.col('EndPos').get(1)).abs()
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

            if condition == 'overlap':
                # The dataframes are sorted by Pos. Check if any alignemnt end
                # is greater than an adjacent alignment start (overlap)
                aln_type_expressions.append(
                    ((pl.col("EndPos") > pl.col("Pos").shift(-1)).any()))

            if len(condition) == 1:
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
            df_subtype.group_by('Read')
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

    # Assign 'transgene_unclassified' Check if any of the original reads are missing
    # from the assigned DF and set to 'transgene_unclassified' if so.
    unknown_df = (
        (
          aln_df.filter(
              ~pl.col('Read').is_in(all_reads_assigned_df['Read']))[
              ['Read']]
        ).unique()
        .with_columns(
            pl.lit(ReadType.transgene_unclassified)
            .alias('Assigned_genome_subtype')
         )
        .select(['Read', 'Assigned_genome_subtype'])
    )
    merged_df = pl.concat([all_reads_assigned_df, unknown_df], how='diagonal')

    # Reads may be assigned to both 'backbone' and 'complex' categories if, for example,
    # they have three or more alignments or overlapping alignments
    # (classified as 'complex'), but also match the criteria for 'backbone'.
    # In these cases, the 'complex' assignment takes priority.
    # Similarly, GDM-assigned reads with two 'partial - no ITR' alignments
    # can also be considered 'complex' if they overlap, so 'complex' should take
    # precedence.

    priority = (
        pl.when(pl.col("Assigned_genome_subtype") == pl.lit(ReadType.complex_))
        .then(0)
        .when(pl.col("Assigned_genome_subtype") == pl.lit(ReadType.bb_contam))
        .then(1)
        .otherwise(2)
    )

    final_df = (
        merged_df
        .with_columns(priority.alias("prior"))
        .sort(["Read", "prior"])
        .group_by("Read")
        .agg(pl.all().first())  # keep the highest-priority assignment
        .select(["Read", "Assigned_genome_subtype"])
    )

    # Assign subtypes to higher level Assigned_genome_type
    per_read_info = final_df.with_columns(
        pl.when(pl.col('Assigned_genome_subtype') == ReadType.bb_contam)
        .then(_l(ReadType.bb_contam))

        .when(pl.col('Assigned_genome_subtype') == ReadType.full_ss)
        .then(_l(ReadType.full_ss))

        .when(pl.col('Assigned_genome_subtype').is_in([
            ReadType.icg5, ReadType.icg3, ReadType.par_icg_incom_itrs,
            ReadType.par_icg_no_itrs, ReadType.par_icg, ReadType.gdm
        ]))
        .then(_l(ReadType.par_ss))

        .when(pl.col('Assigned_genome_subtype') == ReadType.full_sc)
        .then(_l(ReadType.full_sc))

        .when(pl.col('Assigned_genome_subtype').is_in([
            ReadType.itr1_cat, ReadType.itr2_cat, ReadType.itr1_2_cat,
            ReadType.itr3_only, ReadType.itr5_only, ReadType.single_itr]))
        .then(_l(ReadType.itr_only))

        .when(pl.col('Assigned_genome_subtype') == ReadType.complex_)
        .then(_l(ReadType.complex_))

        .when(pl.col('Assigned_genome_subtype').is_in([
            ReadType.sbg5, ReadType.sbg3, ReadType.sbg5_incomp,
            ReadType.sbg3_incomp, ReadType.sbg_unresolved,
            ReadType.sbg5_sym, ReadType.sbg5_asym,
            ReadType.sbg5_incom_sym, ReadType.sbg3_sym, ReadType.sbg3_asym,
            ReadType.sbg5_incom_asym, ReadType.sbg3_incom_sym, ReadType.sbg3_incom_asym
        ]))
        .then(_l(ReadType.par_sc))

        .when(pl.col('Assigned_genome_subtype') == ReadType.transgene_unclassified)
        .then(_l(ReadType.transgene_unclassified))

        .alias('Assigned_genome_type')
    )

    # create summary
    assigned_genome_types_summary = (
        per_read_info.group_by('Assigned_genome_type').len()
        .rename({'len': 'count'})
        .with_columns(
            (pl.col('count') / pl.col('count').sum() * 100).round(2)
            .alias('percentage'))
        .with_columns(pl.lit(sample_id).alias('sample_id'))
    )

    return per_read_info, assigned_genome_types_summary


class BamChunkReader:
    """Iterator that yields batches BAM alignments.

    Identically-named records are guaranteed to be in the same batch as long as
    the input BAM is processed with samtools sort -n, sort-N or collate.
    """

    def __init__(
        self,
        bam_path,
        batch_size=200_000,
        threads=2
    ):
        """Initialize the reader."""
        self.bam_path = bam_path
        self.batch_size = batch_size
        self.threads = threads
        self._overflow_records = []

        self.bam_file = pysam.AlignmentFile(self.bam_path, 'rb', threads=self.threads)
        self._bam_iter = self.bam_file.fetch(until_eof=True)
        self.all_processed = False

    def __iter__(self):
        """Return the iterator object itself."""
        return self

    def __next__(self):
        """Return the next batch of BAM alignments."""
        batch = self._load_next_complete_batch()
        # It's possible to have an empty batch returned if all were unmapped and
        # filtered. There may still be more to process.
        if not batch and self.all_processed:
            raise StopIteration
        return batch

    def _load_next_complete_batch(self):
        this_batch = self._overflow_records
        self._overflow_records = []

        # Fill up to batch_size
        while len(this_batch) < self.batch_size:
            try:
                record = next(self._bam_iter)
                if record.is_unmapped:
                    continue
                this_batch.append(record)
            except StopIteration:
                self.all_processed = True
                break

        if not this_batch:
            # No records mappped reads found.
            return []

        # Ensure full set of records for last read id in batch
        batch_last_id = this_batch[-1].query_name
        while True:
            try:
                peek = next(self._bam_iter)
                if peek.is_unmapped:
                    continue
            except StopIteration:
                break
            if peek.query_name == batch_last_id:
                this_batch.append(peek)
            else:
                self._overflow_records.append(peek)
                break
        return this_batch


def main(args):
    """Entry point."""
    logger = get_named_logger('aav_structures')
    # Get the ITR locations
    itr1_start, itr1_end, itr2_start, itr2_end = args.itr_locations

    schema = {
        'Ref': pl.Utf8,
        'Read': pl.Utf8,
        'Pos': pl.UInt32,
        'EndPos': pl.UInt32,
        'Strand': pl.Int8,
    }

    reader = BamChunkReader(bam_path=args.bam_in)

    # Accumulate summaries across batches
    all_summaries = []

    with contextlib.ExitStack() as stack:

        bam_out = stack.enter_context(pysam.AlignmentFile(
                args.bam_out, 'wb', template=reader.bam_file, threads=2))
        assigned_out = stack.enter_context(open(args.output_per_read, 'a'))

        first_batch = True
        logger.info("Starting processing BAM in chunks\n")
        for batch_n, bam_batch_iter in enumerate(reader, start=1):
            records = [
                [
                    r.reference_name, r.query_name, r.reference_start,
                    r.reference_end, r.is_reverse
                ] for r in bam_batch_iter
            ]

            batch_df = pl.from_records(records, schema)
            logger.info(f"Loaded {len(batch_df)} alignments in batch {batch_n}\n")
            # Get reads that map to transgene plasmid
            non_plasmid_read_ids = (
                batch_df.filter(~pl.col('Ref').is_in([args.transgene_plasmid_name]))
                .select(pl.col('Read')).to_series()
            )
            trans_df = batch_df.filter(~pl.col('Read').is_in(non_plasmid_read_ids))

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
            df_assigned_reads, type_summary = annotate_reads(
                args.sample_id, df_assigned_alns, get_subtype_definitions(),
                args.symmetry_threshold)

            all_summaries.append(type_summary)

            df_assigned_reads.write_csv(
                assigned_out, separator='\t',
                include_header=True if first_batch else False)
            first_batch = False

            # Make a lookup dictionary for read name → assigned genome type
            read_gtype = df_assigned_reads.with_columns(
                pl.col('Assigned_genome_type')
                .str.replace_all(' ', '_')
                .str.to_lowercase()
                .alias('Assigned_genome_type'))
            gtype_lookup = dict(
                zip(
                    read_gtype["Read"].to_list(),
                    read_gtype["Assigned_genome_type"].to_list()))

            # Tag the BAM chunk and append to file
            # This could be put on a separate thread if needed
            for aln in bam_batch_iter:
                gtype = gtype_lookup.get(aln.query_name, 'Non_transgene')
                aln.set_tag('AV', gtype, value_type="Z")
                bam_out.write(aln)

        # Write aggregated summary after all batches processed
        if all_summaries:
            final_summary = (
                pl.concat(all_summaries)
                .group_by(['sample_id', 'Assigned_genome_type'])
                .agg(pl.col('count').sum())
                .with_columns(
                    (pl.col('count') / pl.col('count').sum() * 100).round(2)
                    .alias('percentage'))
            )
            final_summary.write_csv(args.output_plot_data, separator='\t')
