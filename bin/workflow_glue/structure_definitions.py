"""This module contains definitions for AAV genome types.

This module declares the rules governing how AAV genome types are defined based on
the alignments that make up each read. The definitions are used in aav_structures.py
to assign genome types to reads.
The main components are:

  -	Enums for AlnType and ReadType.
  -	AlnType: labels to assign to individual alignments on the transgene plasmid
  -	ReadType: final, per-read categories the workflow reports
"""

from enum import Enum


class AlnType(str, Enum):
    """Enum for Assigning categories to alignments.

    An alignment category defines its ITR and midsection status as well as whether
    the alignment maps to the vector backbone.

    Subclassing str allows us to access the values as strings and not have to
    do .value all over the place.
    """

    # These are alignments that represent almost full AAV genomes. They have varying
    # amounts of ITR on both sides of the alignment and contain full mid-sections
    almost = 'Almost complete'
    full5_par3 = 'Full 5` ITR and partial 3` ITR'
    par5_full3 = 'Partial 5` ITR and full 3` ITR'
    par5_par3 = 'Partial 5` ITR and partial 3` ITR'

    # These alignments are truncated at the mod-section region but contain some
    # ITR region on one of the ends
    full5_par_mid = 'Full 5` ITR and partial mid section'
    par_mid_full3 = 'Partial mid section and full 3` ITR'
    par5_par_mid = 'Partial 5` ITR and partial mid section'
    par_mid_par3 = 'Partial mid section and partial 3` ITR'

    # Only midsection
    par_no_itr = 'Partial - no ITRs'

    # Alignment starts and ends within ITR
    itr5_only = '5` ITR'
    itr3_only = '3` ITR'

    # Transgene plamsid backbone alignments
    vec_bb_5 = 'Vector backbone - 5` end'
    vec_bb_3 = 'Vector backbone - 3` end'
    # Backbone only, no transgene cassette sequence
    bb = 'Backbone'

    ext_itr = 'Extended ITR-ITR region'
    transgene_unclassified = 'Transgene unclassified'
    non_transgene = 'Non transgene'


class ReadType(str, Enum):
    """Enum for ReadType categories.

    Read types are defined by the identity of the constituent alignment types.
    """

    full_ss = 'Full ssAAV'
    full_sc = 'Full scAAV'
    par3_ss = '3` partial ssAAV'
    par5_ss = '5` partial ssAAV'
    bb_int = 'Backbone integration'
    vbb_int = 'Vector-backbone integration'
    par_sc = 'Partial scAAV'
    par_ss = 'Partial ssAAV'
    transgene_unclassified = 'Transgene unclassified'
    non_transgene = 'Non transgene'

    par_icg = "Partial ICG - part of both ITRs"
    icg5 = '5` ICG'
    icg3 = '3` ICG'
    par_icg_incom_itrs = 'Partial ICG - incomplete ITRs'
    par_icg_no_itrs = 'Partial ICG - no ITRs'
    itr5_only = '5` ITR'
    itr3_only = '3` ITR'
    # All the ss and sc categories containing only ITR
    itr_only = 'ITR region only'

    symmetric = 'Symmetric'
    asymmetric = 'Asymmetric'

    sbg_unresolved = 'SBG (unresolved)'
    sbg3 = '3` SBG'
    sbg3_sym = 'SBG 3` symmetric'
    sbg3_asym = 'SBG 3` asymmetric'
    sbg3_incomp = "3` SBG (incomplete)"
    sbg3_incom_asym = 'SBG 3` incomplete ITR asymmetric'
    sbg3_incom_sym = 'SBG 3` incomplete ITR symmetric'

    sbg5 = "5` SBG"
    sbg5_sym = 'SBG 5` symmetric'
    sbg5_asym = 'SBG 5` asymmetric'
    sbg5_incomp = "5` SBG (incomplete)"
    sbg5_incom_asym = 'SBG 5` incomplete ITR asymmetric'
    sbg5_incom_sym = 'SBG 5` incomplete ITR symmetric'

    itr_cat = 'ITR concatemer'
    itr1_cat = 'ITR1 concatemer'
    itr2_cat = 'ITR2 concatemer'
    itr1_2_cat = 'ITR1-ITR2 concatemer'
    single_itr = "ITR single strand"
    potential_cat = 'Potential concatemer'

    bb_contam = 'Backbone contamination'
    complex_ = 'Complex'
    gdm = 'GDM'


class ReadTypeDefinition:
    """Define AAV read types."""

    def __init__(
            self, def_type, n_alignments, n_strands,
            n_aln_types, aln_type_combinations, symmetry=None):
        """
        Initialise a ReadTypeDefinition.

        The AAV type is defined based on the numbers and type of alignments assigned to
        each read. For example:


        >>> TypeDefinition(ReadType.sbg3_incom_sym, 2, 2, 2, [
        >>>        [AlnType.par_mid_full3, AlnType.par_mid_par3]], 'start_sym')


        The above definition defines the ReadType.sbg3_incom_sym (Snapback genome with
        full 3prime ITR and symmetric incomplete 5prime end) as having two
        alignments, alignments on both strands, and two genome types.

        `aln_type_combinations` is a list (a single entry in the example above instance)
        that define the possible additional combinations of AlnTypes that must be
        present for a read to be classed as this ReadType.

        `check_symmetry` defines whether this category includes a symmetry check on
        one side of the read (defaults to None, in which case no check is done) and the
        options are 'start_sym' or 'end_sym', which will check for symmetry/asymmetry of
        alignment between the two strand at the start of the read, and
        'end_sym' or 'end_asym', which specifies a requirement for symmetry/asymmetry
        at the end of the read.

        :param name: ReadType name
        :param n_alignments: Number of alignemnts from read
        :param n_strands: Number of strands reads are on (1/2)
        :param n_aln_types: Number of AlnTypes associated with read
        :param aln_type_combinations: A list of potential combinations of AlnTypes
            that specify a ReadType
        :param symmetry: How to check for symmetrical ends of the alignments
        """
        symmetry_allowed_opts = [None, 'start_sym', 'end_sym', 'start_asym', 'end_asym']
        if symmetry not in symmetry_allowed_opts:
            raise ValueError(
                f'`symmetry` must be one of {symmetry_allowed_opts}'
            )

        self.type_ = def_type
        self.n_alignments = n_alignments
        self.n_strands = n_strands
        self.n_aln_types = n_aln_types
        self.aln_type_combinations = aln_type_combinations
        self.symmetry = symmetry


def get_subtype_definitions():
    """Get TypeDefinitions used to define AAV genome subtypes."""
    definitions = [
        ReadTypeDefinition(ReadType.full_ss, 1, 1, 1, [[AlnType.almost]]),

        # Genome deletion mutants
        # can be different alignment types or the same with a gap
        ReadTypeDefinition(ReadType.gdm, 2, 1, 2, [
                [AlnType.itr5_only, AlnType.par_mid_full3],
                [AlnType.itr5_only, AlnType.par_mid_par3],
                [AlnType.full5_par_mid, AlnType.itr3_only],
                [AlnType.par5_par_mid, AlnType.itr3_only],
                [AlnType.full5_par_mid, AlnType.par_mid_full3],
                [AlnType.full5_par_mid, AlnType.par_mid_par3],
                [AlnType.par5_par_mid, AlnType.par_mid_full3],
                [AlnType.par5_par_mid, AlnType.par_mid_par3]
            ]),
        # Can have two alignments in mid-section that form GDM
        # TODO: we need to be able to specify no overalp here as we are getting complex
        # reads assigned to this type
        ReadTypeDefinition(ReadType.gdm, 2, 1, 1, [
                [AlnType.par_no_itr, AlnType.par_no_itr]
           ]),

        # Genome deletion mutant if there are two alignments, 1 strand, anny number
        # and types of genome type ([]). The condition here is simply no overlap
        # between primary and supplementary alignments.
        # ReadTypeDefinition(ReadType.gdm, 2, 1, 0, ['no_overlap']),

        # Incomplete genomes. ssAAV subtype with either one end missing or just
        # partial ITRs present on both sides. Single alignment
        ReadTypeDefinition(ReadType.icg5, 1, 1, 1, [
                [AlnType.full5_par3],
                [AlnType.full5_par_mid]
            ]),
        ReadTypeDefinition(ReadType.icg3, 1, 1, 1, [
                [AlnType.par5_full3],
                [AlnType.par_mid_full3]
            ]),
        ReadTypeDefinition(ReadType.par_icg_incom_itrs, 1, 1, 1, [
            [AlnType.par5_par_mid],
            [AlnType.par_mid_par3]
        ]),
        ReadTypeDefinition(ReadType.par_icg_no_itrs, 1, 1, 1, [
            [AlnType.par_no_itr]
        ]),

        ReadTypeDefinition(ReadType.par_icg, 1, 1, 1, [
            [AlnType.par5_par3]
        ]),

        # scAAV definitions
        # full_scAAV can have either 1 or two genome types
        #  so there's two definitions
        ReadTypeDefinition(ReadType.full_sc, 2, 2, 1, [
            [AlnType.almost, AlnType.almost],
            [AlnType.par5_full3, AlnType.par5_full3],
            [AlnType.full5_par3, AlnType.full5_par3],
            [AlnType.par5_par3, AlnType.par5_par3]
        ]),

        ReadTypeDefinition(ReadType.full_sc, 2, 2, 2, [
            [AlnType.almost, AlnType.par5_full3],
            [AlnType.almost, AlnType.full5_par3],
            [AlnType.almost, AlnType.par5_par3],
            [AlnType.full5_par3, AlnType.par5_par3],
            [AlnType.par5_full3, AlnType.par5_par3],
            [AlnType.full5_par3, AlnType.par5_full3]
        ]),
        ReadTypeDefinition(ReadType.sbg3_asym, 2, 2, 2, [
            [AlnType.almost, AlnType.par_mid_full3],
            [AlnType.par_mid_full3, AlnType.par5_full3]
        ]),
        ReadTypeDefinition(ReadType.sbg3_incom_asym, 2, 2, 2, [
            [AlnType.almost, AlnType.par_mid_par3],
            [AlnType.full5_par3, AlnType.par_mid_full3],
            [AlnType.par_mid_par3, AlnType.par5_full3],
            [AlnType.par_mid_full3, AlnType.par5_par3],
            [AlnType.full5_par3, AlnType.par_mid_par3],
            [AlnType.par_mid_par3, AlnType.par5_par3]
        ]),
        ReadTypeDefinition(ReadType.sbg3_incom_sym, 2, 2, 1, [
            [AlnType.par_mid_par3, AlnType.par_mid_par3]], symmetry='start_sym'
                           ),
        ReadTypeDefinition(ReadType.sbg3_incom_asym, 2, 2, 1, [
            [AlnType.par_mid_par3, AlnType.par_mid_par3]], symmetry='start_asym'
                           ),
        ReadTypeDefinition(ReadType.sbg3_incom_sym, 2, 2, 2, [
            [AlnType.par_mid_full3, AlnType.par_mid_par3]], symmetry='start_sym'
                           ),
        ReadTypeDefinition(ReadType.sbg3_incom_asym, 2, 2, 2, [
            [AlnType.par_mid_full3, AlnType.par_mid_par3]], symmetry='start_asym'
                           ),
        ReadTypeDefinition(ReadType.sbg5_sym, 2, 2, 1, [
            [AlnType.full5_par_mid, AlnType.full5_par_mid]], symmetry='end_sym'
                           ),
        ReadTypeDefinition(ReadType.sbg5_asym, 2, 2, 2, [
            [AlnType.full5_par_mid, AlnType.full5_par3],
            [AlnType.almost, AlnType.full5_par_mid],
            [AlnType.almost, AlnType.itr5_only]
        ]),
        ReadTypeDefinition(ReadType.sbg5_asym, 2, 2, 1, [
            [AlnType.full5_par_mid, AlnType.full5_par_mid]], symmetry='end_asym'
                           ),
        ReadTypeDefinition(ReadType.sbg5_incom_asym, 2, 2, 2, [
            [AlnType.almost, AlnType.par5_par_mid],
            [AlnType.par5_full3, AlnType.full5_par_mid],
            [AlnType.full5_par_mid, AlnType.par5_par3],
            [AlnType.par5_par_mid, AlnType.full5_par3],
            [AlnType.par5_par_mid, AlnType.par5_par3],
            [AlnType.par5_full3, AlnType.par5_par_mid]
        ]),
        ReadTypeDefinition(ReadType.sbg5_incom_asym, 2, 2, 1, [
            [AlnType.par5_par_mid, AlnType.par5_par_mid]], symmetry='end_asym'
                           ),
        ReadTypeDefinition(ReadType.sbg5_incom_sym, 2, 2, 1, [
            [AlnType.par5_par_mid, AlnType.par5_par_mid]], symmetry='end_sym'
                           ),
        # Same as above but differentiating symmetry
        ReadTypeDefinition(ReadType.sbg5_incom_asym, 2, 2, 2, [
            [AlnType.full5_par_mid, AlnType.par5_par_mid]], symmetry='end_asym'
                           ),
        ReadTypeDefinition(ReadType.sbg5_incom_sym, 2, 2, 2, [
            [AlnType.full5_par_mid, AlnType.par5_par_mid]], symmetry='end_sym'
                           ),
        ReadTypeDefinition(ReadType.sbg3_sym, 2, 2, 1, [
            [AlnType.par_mid_full3, AlnType.par_mid_full3]], symmetry='start_sym'
                           ),
        ReadTypeDefinition(ReadType.sbg3_asym, 2, 2, 1, [
            [AlnType.par_mid_full3, AlnType.par_mid_full3]], symmetry='start_asym'
                           ),
        # scAAV with ITRs on one side only
        ReadTypeDefinition(ReadType.sbg_unresolved, 2, 2, 2, [
            [AlnType.par5_par3, AlnType.par_no_itr],
            [AlnType.par5_par_mid, AlnType.par_mid_par3],
            [AlnType.full5_par_mid, AlnType.par_mid_par3],
            [AlnType.par5_par_mid, AlnType.par_mid_full3],
            [AlnType.full5_par_mid, AlnType.itr3_only],
            [AlnType.par_mid_full3, AlnType.itr5_only],
            [AlnType.par5_par_mid, AlnType.itr3_only],
            [AlnType.par_mid_par3, AlnType.itr5_only],
            [AlnType.itr5_only, AlnType.par_no_itr],
            [AlnType.itr3_only, AlnType.par_no_itr],
            [AlnType.full5_par_mid, AlnType.par_mid_full3],
            [AlnType.almost, AlnType.itr3_only],
            [AlnType.par_mid_full3, AlnType.itr3_only],
            [AlnType.par5_full3, AlnType.itr3_only],
            [AlnType.par_mid_full3, AlnType.par_no_itr],
            [AlnType.par5_full3, AlnType.par_no_itr],
            [AlnType.par_mid_par3, AlnType.par_no_itr],
            [AlnType.par_mid_par3, AlnType.itr3_only],
            [AlnType.full5_par3, AlnType.itr3_only],
            [AlnType.par5_par3, AlnType.itr3_only],
            [AlnType.par5_par3, AlnType.itr5_only],
            [AlnType.full5_par3, AlnType.itr5_only],
            [AlnType.full5_par_mid, AlnType.par_no_itr],
            [AlnType.full5_par3, AlnType.par_no_itr],
            [AlnType.par5_par_mid, AlnType.par_no_itr],
            [AlnType.par5_par_mid, AlnType.itr5_only],
            [AlnType.full5_par_mid, AlnType.itr5_only],
            [AlnType.par5_full3, AlnType.itr5_only],
            [AlnType.par5_par3, AlnType.itr5_only],
            [AlnType.almost, AlnType.par_no_itr]
        ]),

        ReadTypeDefinition(ReadType.sbg_unresolved, 2, 2, 1, [
            [AlnType.par_no_itr]
        ]),

        # ITR concatemers
        ReadTypeDefinition(ReadType.itr1_cat, 2, 2, 1, [
            [AlnType.itr5_only, AlnType.itr5_only]
        ]),
        ReadTypeDefinition(ReadType.itr2_cat, 2, 2, 1, [
            [AlnType.itr3_only, AlnType.itr3_only]
        ]),
        ReadTypeDefinition(ReadType.itr1_2_cat, 2, 2, 2, [
            [AlnType.itr5_only, AlnType.itr3_only]
        ]),

        # Single ITR
        ReadTypeDefinition(ReadType.single_itr, 1, 1, 1, [
            [AlnType.itr5_only],
            [AlnType.itr3_only]
        ]),

        # Backbone contamination: 1 alignment
        ReadTypeDefinition(ReadType.bb_contam, 1, 1, 1, [
            [AlnType.vec_bb_5],
            [AlnType.vec_bb_3],
            [AlnType.bb]
        ]),
        # Backbone contamination: 2 alignments, 1 or 2 strands, 1 or two AlnTypes
        ReadTypeDefinition(ReadType.bb_contam, 2, 0, 0, [
            [AlnType.vec_bb_5],
            [AlnType.vec_bb_3],
            [AlnType.bb]
        ]),
        # If we have primary and supplementary overlapping, no sure what's going on
        # so we'll class as complex
        # Also extra definition in aav_structures.py for when number of alignments > 2
        ReadTypeDefinition(ReadType.complex_, 2, 1, 0, [
            'overlap',
        ]),
    ]

    return definitions
