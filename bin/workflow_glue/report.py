"""Create workflow report."""
import json

from dominate.tags import h5, p
import ezcharts as ezc
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.layout.snippets.table import DataTable
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def plot_trucations(report, truncations_file):
    """Make report section with truncation plots."""
    df = pd.read_csv(truncations_file, sep='\t')

    with report.add_section("Truncations", "Truncations"):
        p(
            "The plots show the frequency of start and end positions of reads "
            "that map completely within the ITR-ITR region of the transgene plasmid, "
            "helping to identify truncation hotspots"
        )
        tabs = Tabs()

        for sample, df_sample in df.groupby('sample_id'):
            with tabs.add_tab(sample):
                with Grid(columns=2):
                    df_sample.drop(columns=['sample_id'], inplace=True)
                    plt = ezc.histplot(data=df_sample, binwidth=5)
                    plt.title = dict(text='Start positions')
                    plt.xAxis = dict(name='Genomic position')
                    plt.yAxis = dict(name='Number of alignments')
                    plt.legend = dict(orient='horizontal', top=30)
                    EZChart(plt, theme='epi2melabs')


def plot_itr_coverage(report, coverage_file):
    """Make report section with ITR-ITR coverage."""
    df = pd.read_csv(coverage_file, sep=r"\s+")

    with report.add_section("ITR-ITR coverage", "Coverage"):
        p(
            "For each transgene reference the sequencing depth is calculated "
            "for both forward and reverse mapping reds. "
            "All primary and secondary reads "
            "Note: no MAPQ filtering at the moment."
        )
        tabs = Tabs()
        for sample, df_sample in df.groupby('sample_id'):
            with tabs.add_tab(sample):
                with Grid(columns=2):
                    for ref, df_ref in df_sample.groupby('ref'):
                        plt = ezc.lineplot(
                            data=df_ref, x='pos', y='depth', hue='strand')
                        plt.title = dict(text=ref)
                        plt.legend = dict(
                            orient='horizontal', top=30, icon='rect')
                        for s in plt.series:
                            s.showSymbol = False
                        EZChart(plt, theme='epi2melabs', height='300px')


def plot_contamination(report, class_counts, read_lengths):
    """Make report section with contamination plots."""
    df_class_counts = pd.read_csv(class_counts, sep='\t')
    # df_len = pd.read_csv(read_lengths, sep='\t')

    with report.add_section("Contamination", "Contam."):
        p(
            "This analysis details the references that each read align to:  "
            "human, helper plasmid, rep-cap plasmid and transgene plasmid "
            "(or unknown if do not map to any of the given references). "
            "For reads that map to multiple references, they will be counted once "
            "for each alignment as the counts are made per-alignment not per-read."
        )
        p(
            "Read depth distribution histograms are also displayed per reference."
        )

        tabs = Tabs()
        for sample, df_sample in df_class_counts.groupby('sample_id'):
            with tabs.add_tab(sample):
                # df_sample.drop(columns=['sample_id'], inplace=True)
                plt = ezc.barplot(df_sample[['Reference', 'Percentage of alignments']])
                h5('Reference mapping counts')
                EZChart(plt, theme='epi2melabs', height='400px')
                # read length plots in same view
                # df_sample_lengths = df_len[df_len.sample_id == sample]
                # h5('Contamination read lengths')
                # with Grid(columns=4):
                #     for cat, df_cat_len in df_sample_lengths.groupby(
                #             'contam_class_pct'):
                #         plt = ezc.histplot(data=df_cat_len['ReadLen'])
                #         plt.title = dict(text=cat)
                #         plt.xAxis = dict(name='Genomic position')
                #         plt.yAxis = dict(name='Percentage of alignments')
                #         EZChart(plt, theme='epi2melabs', height='300px')


def plot_aav_structures(report, structures_file):
    """Make report section barplots detailing the AAV structures found."""
    df = pd.read_csv(structures_file, sep='\t')

    with report.add_section("AAV Structures", "Structures"):
        p(
            "The different possible transgene plasmid structure are detailed in the "
            "following plots."
        )
        tabs = Tabs()
        for sample, df_sample in df.groupby('sample_id'):
            with tabs.add_tab(sample):
                df_sample = df_sample.sort_values('percentage', ascending=False)
                # Plot of main genome type counts
                plt = ezc.barplot(
                    df_sample,
                    x='Assigned_genome_type',
                    y='percentage')
                plt.title = dict(text='Genome types')
                plt.xAxis.axisLabel = dict(rotate=45)
                EZChart(plt, theme='epi2melabs')

                # Table with counts and percentages
                # (in lieu of being able to annotate bar plots in ezchrts)
                DataTable.from_pandas(df_sample, use_index=False)


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "AAV QC workflow report", "wf-aav-qc",
        args.params, args.versions)

    with open(args.metadata) as metadata:
        sample_details = sorted([
            {
                'sample': d['alias'],
                'type': d['type'],
                'barcode': d['barcode']
            } for d in json.load(metadata)
        ], key=lambda d: d["sample"])

    if args.stats:
        with report.add_section("Read summary", "Read summary"):
            # TODO fix this. Do we need o concat stats?
            fastcat.SeqSummary(args.stats[0])

    plot_contamination(
        report,
        args.contam_class_counts,
        args.contam_read_lengths)
    plot_trucations(report, args.truncations)
    plot_itr_coverage(report, args.itr_coverage)
    plot_aav_structures(report, args.aav_structures)

    with report.add_section("Metadata", "Metadata"):
        tabs = Tabs()
        for d in sample_details:
            with tabs.add_tab(d["sample"]):
                df = pd.DataFrame.from_dict(
                    d, orient="index", columns=["Value"])
                df.index.name = "Key"
                DataTable.from_pandas(df)

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "report", help="Report output file")
    parser.add_argument(
        "--stats", nargs='*', help="Fastcat per-read stats file(s).")
    parser.add_argument(
        "--truncations", help="TSV with ? for truncations.")
    parser.add_argument(
        "--itr_coverage", help="TSV with ? for truncations.")
    parser.add_argument(
        "--contam_class_counts", help="TSV with ? for contamination.")
    parser.add_argument(
        "--contam_read_lengths", help="TSV with ? for contamination.")
    parser.add_argument(
        "--aav_structures", help="TSV of reads with AAV structure assignment.")
    parser.add_argument(
        "--consensus_accuracy",
        help="TSV accuracies per transgene plasmid orientation.")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    return parser
