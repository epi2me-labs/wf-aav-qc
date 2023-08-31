"""Common code to be used acrosss workflow scripts."""

import polars as pl
# The refrence names of the trans plasmid references in different
# possible configurations
trans_plasmid_ref_names = ['flip_flip', 'flip_flop', 'flop_flip', 'flop_flop']


def load_annotation(annotation):
    """Load the annotation BED file."""
    ann_schema = {
        'ref': str,
        'start': pl.UInt32,
        'end': pl.UInt32,
        'feature': str
    }
    df_ann = pl.read_csv(
        source=annotation,
        has_header=False,
        separator='\t',
        dtypes=ann_schema
    )
    return df_ann
