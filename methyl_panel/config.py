"""
config.py
=========
Comprehensive Primer3Plus-compatible configuration for the WBC methylation
panel pipeline.

Every setting from Primer3Plus (Main, General Settings, Advanced Settings,
Internal Oligo, Penalties, Advanced Seq.) is exposed as a configurable
parameter with the exact same name as in the Primer3Plus settings file.

Two salt-condition presets are provided:
  - "roche_dlc":  100 mM Na, 4.5 mM Mg, 1.2 mM dNTP, 250 nM primer (Roche DLC)
  - "primer3plus": 50 mM Na, 3 mM Mg, 1.2 mM dNTP, 250 nM primer (Primer3Plus default)

The Tm window is tightened to 58-62 C per the updated design rules.
"""

from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional


# ─────────────────────────────────────────────────────────────────────────────
# Universal TAG sequence (LRP5 105 bp) for QIAcuity dPCR TAG channel
# ─────────────────────────────────────────────────────────────────────────────
TAG_SEQUENCE = (
    "AGGTGCCCATTGGTATTGCGGCCGTGCGCCCGGCGGGCATGAATTAGCTGTGCCGCCTGG"
    "CTGCTGACGGGACGCCTCGCCTCGACTGAAAACTACCTGGAGCTG"
)

# T cell subtypes that require CD3T in background
T_CELL_SUBTYPES = {
    "Blood-T-CD3", "Blood-T-CD4", "Blood-T-CD8",
    "Blood-T-CenMem-CD4", "Blood-T-EffMem-CD4", "Blood-T-EffMem-CD8",
    "Blood-T-Eff-CD8", "Blood-T-Naive-CD4", "Blood-T-Naive-CD8",
    "CD4T", "CD8T", "CD3T", "CD4", "CD8",
}


# ─────────────────────────────────────────────────────────────────────────────
# Salt condition presets
# ─────────────────────────────────────────────────────────────────────────────
SALT_PRESETS = {
    "roche_dlc": {
        "PRIMER_SALT_MONOVALENT": 100.0,   # Na+ mM
        "PRIMER_SALT_DIVALENT": 4.5,       # Mg2+ mM
        "PRIMER_DNTP_CONC": 1.2,           # dNTPs mM
        "PRIMER_DNA_CONC": 250.0,          # primer nM
        "PRIMER_SALT_CORRECTIONS": 1,      # 1 = SantaLucia 1998 (Owczarzy not supported by design_primers)
        "PRIMER_TM_FORMULA": 1,            # 1 = SantaLucia 1998
    },
    "primer3plus": {
        "PRIMER_SALT_MONOVALENT": 50.0,    # Na+ mM
        "PRIMER_SALT_DIVALENT": 3.0,       # Mg2+ mM
        "PRIMER_DNTP_CONC": 1.2,           # dNTPs mM
        "PRIMER_DNA_CONC": 250.0,          # primer nM
        "PRIMER_SALT_CORRECTIONS": 1,      # 1 = SantaLucia 1998 salt correction
        "PRIMER_TM_FORMULA": 1,            # 1 = SantaLucia 1998
    },
}


@dataclass
class Primer3PlusConfig:
    """
    Complete Primer3Plus-compatible configuration.

    Every parameter from the Primer3Plus settings file is exposed here,
    organized by the Primer3Plus UI tabs:

    === Main ===
        PRIMER_TASK, PRIMER_PICK_LEFT_PRIMER, PRIMER_PICK_RIGHT_PRIMER,
        PRIMER_PICK_INTERNAL_OLIGO, PRIMER_NUM_RETURN, PRIMER_PRODUCT_SIZE_RANGE

    === General Settings ===
        PRIMER_MIN_SIZE, PRIMER_OPT_SIZE, PRIMER_MAX_SIZE,
        PRIMER_MIN_TM, PRIMER_OPT_TM, PRIMER_MAX_TM,
        PRIMER_MIN_GC, PRIMER_OPT_GC_PERCENT, PRIMER_MAX_GC,
        PRIMER_MAX_POLY_X

    === Advanced Settings ===
        PRIMER_DNA_CONC, PRIMER_SALT_MONOVALENT, PRIMER_SALT_DIVALENT,
        PRIMER_DNTP_CONC, PRIMER_SALT_CORRECTIONS, PRIMER_TM_FORMULA,
        PRIMER_ANNEALING_TEMP, PRIMER_MAX_END_STABILITY, PRIMER_MAX_END_GC,
        PRIMER_GC_CLAMP, PRIMER_MAX_NS_ACCEPTED, PRIMER_LIBERAL_BASE,
        PRIMER_LOWERCASE_MASKING, PRIMER_FIRST_BASE_INDEX,
        PRIMER_MAX_SELF_ANY, PRIMER_MAX_SELF_ANY_TH, PRIMER_MAX_SELF_END,
        PRIMER_MAX_SELF_END_TH, PRIMER_MAX_HAIRPIN_TH,
        PRIMER_MAX_TEMPLATE_MISPRIMING, PRIMER_MAX_TEMPLATE_MISPRIMING_TH,
        PRIMER_MAX_LIBRARY_MISPRIMING, PRIMER_MISPRIMING_LIBRARY,
        PRIMER_PAIR_MAX_COMPL_ANY, PRIMER_PAIR_MAX_COMPL_ANY_TH,
        PRIMER_PAIR_MAX_COMPL_END, PRIMER_PAIR_MAX_COMPL_END_TH,
        PRIMER_PAIR_MAX_DIFF_TM, PRIMER_PAIR_MAX_LIBRARY_MISPRIMING,
        PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING, PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH,
        PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT, PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT,
        PRIMER_SECONDARY_STRUCTURE_ALIGNMENT,
        PRIMER_DMSO_CONC, PRIMER_DMSO_FACTOR, PRIMER_FORMAMIDE_CONC,
        PRIMER_PRODUCT_OPT_SIZE, PRIMER_PRODUCT_MAX_TM, PRIMER_PRODUCT_MIN_TM,
        PRIMER_PRODUCT_OPT_TM, PRIMER_PICK_ANYWAY,
        PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE, PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE,
        PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION, PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION,
        PRIMER_QUALITY_RANGE_MIN, PRIMER_QUALITY_RANGE_MAX,
        PRIMER_SEQUENCING_ACCURACY, PRIMER_SEQUENCING_INTERVAL,
        PRIMER_SEQUENCING_LEAD, PRIMER_SEQUENCING_SPACING

    === Internal Oligo ===
        PRIMER_INTERNAL_MIN_SIZE, PRIMER_INTERNAL_OPT_SIZE, PRIMER_INTERNAL_MAX_SIZE,
        PRIMER_INTERNAL_MIN_TM, PRIMER_INTERNAL_OPT_TM, PRIMER_INTERNAL_MAX_TM,
        PRIMER_INTERNAL_MIN_GC, PRIMER_INTERNAL_OPT_GC_PERCENT, PRIMER_INTERNAL_MAX_GC,
        PRIMER_INTERNAL_MAX_POLY_X, PRIMER_INTERNAL_MAX_NS_ACCEPTED,
        PRIMER_INTERNAL_MAX_SELF_ANY, PRIMER_INTERNAL_MAX_SELF_ANY_TH,
        PRIMER_INTERNAL_MAX_SELF_END, PRIMER_INTERNAL_MAX_SELF_END_TH,
        PRIMER_INTERNAL_MAX_HAIRPIN_TH, PRIMER_INTERNAL_MAX_LIBRARY_MISHYB,
        PRIMER_INTERNAL_MISHYB_LIBRARY,
        PRIMER_INTERNAL_DNA_CONC, PRIMER_INTERNAL_SALT_MONOVALENT,
        PRIMER_INTERNAL_SALT_DIVALENT, PRIMER_INTERNAL_DNTP_CONC,
        PRIMER_INTERNAL_DMSO_CONC, PRIMER_INTERNAL_DMSO_FACTOR,
        PRIMER_INTERNAL_FORMAMIDE_CONC,
        PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION,
        PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION

    === Penalties ===
        PRIMER_INSIDE_PENALTY, PRIMER_OUTSIDE_PENALTY,
        PRIMER_WT_TM_GT, PRIMER_WT_TM_LT,
        PRIMER_WT_SIZE_GT, PRIMER_WT_SIZE_LT,
        PRIMER_WT_GC_PERCENT_GT, PRIMER_WT_GC_PERCENT_LT,
        PRIMER_WT_SELF_ANY, PRIMER_WT_SELF_ANY_TH,
        PRIMER_WT_SELF_END, PRIMER_WT_SELF_END_TH,
        PRIMER_WT_HAIRPIN_TH,
        PRIMER_WT_END_STABILITY, PRIMER_WT_END_QUAL,
        PRIMER_WT_NUM_NS, PRIMER_WT_SEQ_QUAL,
        PRIMER_WT_LIBRARY_MISPRIMING, PRIMER_WT_TEMPLATE_MISPRIMING,
        PRIMER_WT_TEMPLATE_MISPRIMING_TH,
        PRIMER_WT_POS_PENALTY, PRIMER_WT_BOUND_GT, PRIMER_WT_BOUND_LT,
        PRIMER_PAIR_WT_DIFF_TM, PRIMER_PAIR_WT_COMPL_ANY,
        PRIMER_PAIR_WT_COMPL_ANY_TH, PRIMER_PAIR_WT_COMPL_END,
        PRIMER_PAIR_WT_COMPL_END_TH,
        PRIMER_PAIR_WT_PR_PENALTY, PRIMER_PAIR_WT_IO_PENALTY,
        PRIMER_PAIR_WT_PRODUCT_SIZE_GT, PRIMER_PAIR_WT_PRODUCT_SIZE_LT,
        PRIMER_PAIR_WT_PRODUCT_TM_GT, PRIMER_PAIR_WT_PRODUCT_TM_LT,
        PRIMER_PAIR_WT_LIBRARY_MISPRIMING,
        PRIMER_PAIR_WT_TEMPLATE_MISPRIMING, PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH,
        PRIMER_INTERNAL_WT_TM_GT, PRIMER_INTERNAL_WT_TM_LT,
        PRIMER_INTERNAL_WT_SIZE_GT, PRIMER_INTERNAL_WT_SIZE_LT,
        PRIMER_INTERNAL_WT_GC_PERCENT_GT, PRIMER_INTERNAL_WT_GC_PERCENT_LT,
        PRIMER_INTERNAL_WT_SELF_ANY, PRIMER_INTERNAL_WT_SELF_ANY_TH,
        PRIMER_INTERNAL_WT_SELF_END, PRIMER_INTERNAL_WT_SELF_END_TH,
        PRIMER_INTERNAL_WT_HAIRPIN_TH,
        PRIMER_INTERNAL_WT_NUM_NS, PRIMER_INTERNAL_WT_SEQ_QUAL,
        PRIMER_INTERNAL_WT_END_QUAL,
        PRIMER_INTERNAL_WT_LIBRARY_MISHYB, PRIMER_INTERNAL_WT_BOUND_GT,
        PRIMER_INTERNAL_WT_BOUND_LT

    === Advanced Seq. ===
        PRIMER_MUST_MATCH_FIVE_PRIME, PRIMER_MUST_MATCH_THREE_PRIME,
        PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME,
        PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME,
        PRIMER_MIN_QUALITY, PRIMER_MIN_END_QUALITY,
        PRIMER_INTERNAL_MIN_QUALITY,
        PRIMER_MIN_BOUND, PRIMER_MAX_BOUND, PRIMER_OPT_BOUND,
        PRIMER_INTERNAL_MIN_BOUND, PRIMER_INTERNAL_MAX_BOUND,
        PRIMER_INTERNAL_OPT_BOUND
    """

    # ── Main ────────────────────────────────────────────────────────────────
    primer_task: str = "generic"
    primer_pick_left_primer: int = 1
    primer_pick_right_primer: int = 1
    primer_pick_internal_oligo: int = 0
    primer_num_return: int = 200  # Return all valid pairs (was 10)
    # Product size range as a list of (min, max) tuples
    primer_product_size_range: List[tuple] = field(
        default_factory=lambda: [(60, 150)]
    )

    # ── General Settings ────────────────────────────────────────────────────
    # Primer length
    primer_min_size: int = 18
    primer_opt_size: int = 20
    primer_max_size: int = 23

    # Tm — REQUIRED user input. No defaults: the pipeline refuses to run
    # until the user sets these three values. This forces an explicit decision
    # about the Tm window for each run.
    # Example: primer_min_tm=58.0, primer_opt_tm=60.0, primer_max_tm=62.0
    primer_min_tm: Optional[float] = None
    primer_opt_tm: Optional[float] = None
    primer_max_tm: Optional[float] = None

    # GC content
    primer_min_gc: float = 30.0
    primer_opt_gc_percent: float = 50.0
    primer_max_gc: float = 80.0

    # Max poly-X (consecutive identical bases)
    primer_max_poly_x: int = 5

    # ── Advanced Settings ────────────────────────────────────────────────────
    # Salt conditions — default to Primer3Plus settings (50/3/1.2/250nM)
    # Use salt_preset="roche_dlc" for Roche DLC conditions (100/4.5/1.2/250nM)
    salt_preset: str = "primer3plus"
    primer_dna_conc: float = 250.0       # primer concentration nM
    primer_salt_monovalent: float = 50.0 # Na+ mM
    primer_salt_divalent: float = 3.0    # Mg2+ mM
    primer_dntp_conc: float = 1.2        # dNTPs mM
    primer_salt_corrections: int = 1     # 1 = SantaLucia, 2 = Owczarzy 2008
    primer_tm_formula: int = 1           # 1 = SantaLucia 1998

    # Annealing temperature (for Tm alignment calculations)
    primer_annealing_temp: float = 52.0

    # End stability and GC clamp
    primer_max_end_stability: float = 9.0
    primer_max_end_gc: int = 2
    primer_gc_clamp: int = 1

    # Ns and ambiguity
    primer_max_ns_accepted: int = 0
    primer_liberal_base: int = 1
    primer_lowercase_masking: int = 0
    primer_first_base_index: int = 1

    # Self-complementarity and hairpin (TH = thermodynamic alignment)
    primer_max_self_any: float = 8.0
    primer_max_self_any_th: float = 47.0
    primer_max_self_end: float = 3.0
    primer_max_self_end_th: float = 47.0
    primer_max_hairpin_th: float = 47.0

    # Template mispriming
    primer_max_template_mispriming: float = 12.0
    primer_max_template_mispriming_th: float = 47.0

    # Library mispriming
    primer_max_library_mispriming: float = 12.0
    primer_mispriming_library: str = "humrep_and_simple.txt"

    # Pair complementarity
    primer_pair_max_compl_any: float = 8.0
    primer_pair_max_compl_any_th: float = 47.0
    primer_pair_max_compl_end: float = 3.0
    primer_pair_max_compl_end_th: float = 47.0
    primer_pair_max_diff_tm: float = 1.0
    primer_pair_max_library_mispriming: float = 24.0
    primer_pair_max_template_mispriming: float = 24.0
    primer_pair_max_template_mispriming_th: float = 47.0

    # Thermodynamic alignment
    primer_thermodynamic_oligo_alignment: int = 0
    primer_thermodynamic_template_alignment: int = 0
    primer_secondary_structure_alignment: int = 1

    # Additives
    primer_dmso_conc: float = 0.0
    primer_dmso_factor: float = 0.6
    primer_formamide_conc: float = 0.0

    # Product Tm
    primer_product_opt_size: int = 0
    primer_product_max_tm: float = 1000000.0
    primer_product_min_tm: float = -1000000.0
    primer_product_opt_tm: float = 0.0

    # Pick anyway (ignore quality flags)
    primer_pick_anyway: int = 0

    # Three-prime distance
    primer_min_left_three_prime_distance: int = 3
    primer_min_right_three_prime_distance: int = 3

    # Junction overlap
    primer_min_3_prime_overlap_of_junction: int = 4
    primer_min_5_prime_overlap_of_junction: int = 7

    # Quality range
    primer_quality_range_min: int = 0
    primer_quality_range_max: int = 100

    # Sequencing parameters
    primer_sequencing_accuracy: int = 20
    primer_sequencing_interval: int = 250
    primer_sequencing_lead: int = 50
    primer_sequencing_spacing: int = 500

    # ── Internal Oligo (Hydrolysis Probe) ────────────────────────────────────
    primer_internal_min_size: int = 14
    primer_internal_opt_size: int = 20
    primer_internal_max_size: int = 25

    primer_internal_min_tm: float = 63.0
    primer_internal_opt_tm: float = 65.0
    primer_internal_max_tm: float = 68.0

    primer_internal_min_gc: float = 20.0
    primer_internal_opt_gc_percent: float = 50.0
    primer_internal_max_gc: float = 80.0

    primer_internal_max_poly_x: int = 5
    primer_internal_max_ns_accepted: int = 0

    primer_internal_max_self_any: float = 12.0
    primer_internal_max_self_any_th: float = 47.0
    primer_internal_max_self_end: float = 12.0
    primer_internal_max_self_end_th: float = 47.0
    primer_internal_max_hairpin_th: float = 47.0
    primer_internal_max_library_mishyb: float = 12.0
    primer_internal_mishyb_library: str = "humrep_and_simple.txt"

    # Internal oligo salt conditions
    primer_internal_dna_conc: float = 100.0
    primer_internal_salt_monovalent: float = 50.0
    primer_internal_salt_divalent: float = 3.0
    primer_internal_dntp_conc: float = 1.2
    primer_internal_dmso_conc: float = 0.0
    primer_internal_dmso_factor: float = 0.6
    primer_internal_formamide_conc: float = 0.0

    # Internal oligo junction overlap
    primer_internal_min_3_prime_overlap_of_junction: int = 4
    primer_internal_min_5_prime_overlap_of_junction: int = 7

    # ── Penalties ────────────────────────────────────────────────────────────
    # Position penalties
    primer_inside_penalty: float = -1.0
    primer_outside_penalty: float = 0.0

    # Tm penalties
    primer_wt_tm_gt: float = 1.0
    primer_wt_tm_lt: float = 1.0

    # Size penalties
    primer_wt_size_gt: float = 1.0
    primer_wt_size_lt: float = 1.0

    # GC penalties
    primer_wt_gc_percent_gt: float = 0.5
    primer_wt_gc_percent_lt: float = 0.5

    # Self-complementarity penalties
    primer_wt_self_any: float = 0.0
    primer_wt_self_any_th: float = 0.0
    primer_wt_self_end: float = 0.0
    primer_wt_self_end_th: float = 0.0
    primer_wt_hairpin_th: float = 0.0

    # End stability / quality penalties
    primer_wt_end_stability: float = 0.0
    primer_wt_end_qual: float = 0.0
    primer_wt_num_ns: float = 0.0
    primer_wt_seq_qual: float = 0.0

    # Library/template mispriming penalties
    primer_wt_library_mispriming: float = 0.0
    primer_wt_template_mispriming: float = 0.0
    primer_wt_template_mispriming_th: float = 0.0

    # Position / bound penalties
    primer_wt_pos_penalty: float = 0.0
    primer_wt_bound_gt: float = 0.0
    primer_wt_bound_lt: float = 0.0

    # Pair penalties
    primer_pair_wt_diff_tm: float = 0.0
    primer_pair_wt_compl_any: float = 0.0
    primer_pair_wt_compl_any_th: float = 0.0
    primer_pair_wt_compl_end: float = 0.0
    primer_pair_wt_compl_end_th: float = 0.0
    primer_pair_wt_pr_penalty: float = 1.0
    primer_pair_wt_io_penalty: float = 0.0
    primer_pair_wt_product_size_gt: float = 0.0
    primer_pair_wt_product_size_lt: float = 0.0
    primer_pair_wt_product_tm_gt: float = 0.0
    primer_pair_wt_product_tm_lt: float = 0.0
    primer_pair_wt_library_mispriming: float = 0.0
    primer_pair_wt_template_mispriming: float = 0.0
    primer_pair_wt_template_mispriming_th: float = 0.0

    # Internal oligo penalties
    primer_internal_wt_tm_gt: float = 1.0
    primer_internal_wt_tm_lt: float = 1.0
    primer_internal_wt_size_gt: float = 1.0
    primer_internal_wt_size_lt: float = 1.0
    primer_internal_wt_gc_percent_gt: float = 0.0
    primer_internal_wt_gc_percent_lt: float = 0.0
    primer_internal_wt_self_any: float = 0.0
    primer_internal_wt_self_any_th: float = 0.0
    primer_internal_wt_self_end: float = 0.0
    primer_internal_wt_self_end_th: float = 0.0
    primer_internal_wt_hairpin_th: float = 0.0
    primer_internal_wt_num_ns: float = 0.0
    primer_internal_wt_seq_qual: float = 0.0
    primer_internal_wt_end_qual: float = 0.0
    primer_internal_wt_library_mishyb: float = 0.0
    primer_internal_wt_bound_gt: float = 0.0
    primer_internal_wt_bound_lt: float = 0.0

    # ── Advanced Seq. ────────────────────────────────────────────────────────
    # Must-match constraints (empty string = no constraint)
    primer_must_match_five_prime: str = ""
    primer_must_match_three_prime: str = ""
    primer_internal_must_match_five_prime: str = ""
    primer_internal_must_match_three_prime: str = ""

    # Quality thresholds
    primer_min_quality: int = 0
    primer_min_end_quality: int = 0
    primer_internal_min_quality: int = 0

    # Bound (Tm alignment) parameters
    primer_min_bound: float = -10.0
    primer_max_bound: float = 110.0
    primer_opt_bound: float = 97.0
    primer_internal_min_bound: float = -10.0
    primer_internal_max_bound: float = 110.0
    primer_internal_opt_bound: float = 97.0

    # ── Bisulfite-specific settings (not in Primer3Plus) ─────────────────────
    # These are pipeline-specific extensions for bisulfite MSP primer design
    min_cpg_per_primer: int = 2          # hard minimum CpGs per primer
    recommended_cpg_per_primer: int = 4  # below this = WARNING flag only
    min_cpg_pair_total: int = 4          # hard minimum total CpGs
    require_terminal_cpg: bool = True    # at least one primer must have dist=0
    primer_flank_bp: int = 200           # flanking bases around block for Primer3

    # CpG placement mode for methylation discrimination:
    #   "strict_terminal" — CpGs counted in last 5 nt of 3' end (default,
    #       bisulfite-PCR mode)
    #   "three_prime_half" — CpGs counted in last 50% of primer length
    #       (MSP mode, broader discrimination window)
    cpg_placement_mode: str = "strict_terminal"

    # ── New QC filter thresholds ─────────────────────────────────────────────
    # Secondary structure (RNAfold DNA params / primer3 calc_hairpin)
    structure_mfe_cutoff: float = -1.5   # kcal/mol; reject if MFE <= this

    # Primer-dimer (DimerDetective)
    dimer_end_min_dg_cutoff: float = -1.0  # kcal/mol; reject if end_min_dg <= this
    # DimerDetective thermodynamic parameters
    dimer_mv_conc: float = 50.0     # mM monovalent
    dimer_dv_conc: float = 1.5      # mM divalent
    dimer_dntp_conc: float = 0.6    # mM dNTP
    dimer_dna_conc: float = 50.0    # nM DNA
    dimer_temp_c: float = 37.0      # temperature C

    # Common SNP screening
    snp_maf_threshold: float = 0.01  # MAF >= 1% = common
    snp_reject_3prime: bool = True   # reject SNPs in 3'-terminal 5 nt
    snp_reject_any: bool = False     # also reject SNPs anywhere in primer

    def __post_init__(self):
        """Apply salt preset and validate required Tm settings."""
        if self.salt_preset in SALT_PRESETS:
            preset = SALT_PRESETS[self.salt_preset]
            self.primer_salt_monovalent = preset["PRIMER_SALT_MONOVALENT"]
            self.primer_salt_divalent = preset["PRIMER_SALT_DIVALENT"]
            self.primer_dntp_conc = preset["PRIMER_DNTP_CONC"]
            self.primer_dna_conc = preset["PRIMER_DNA_CONC"]
            self.primer_salt_corrections = preset["PRIMER_SALT_CORRECTIONS"]
            self.primer_tm_formula = preset["PRIMER_TM_FORMULA"]

    def validate_tm(self) -> None:
        """
        Validate that Tm window is set. Call this before any primer design.

        Raises ValueError if primer_min_tm, primer_opt_tm, or primer_max_tm
        is None. The Tm window is a required user input — there are no
        defaults because the optimal window depends on the platform and
        experiment design.

        Also validates that min <= opt <= max.
        """
        missing = []
        if self.primer_min_tm is None:
            missing.append("primer_min_tm")
        if self.primer_opt_tm is None:
            missing.append("primer_opt_tm")
        if self.primer_max_tm is None:
            missing.append("primer_max_tm")
        if missing:
            raise ValueError(
                f"Tm window not set. The following required parameters are None: "
                f"{', '.join(missing)}. "
                f"Set them before running the pipeline, e.g.:\n"
                f"  config.primer3.primer_min_tm = 58.0\n"
                f"  config.primer3.primer_opt_tm = 60.0\n"
                f"  config.primer3.primer_max_tm = 62.0\n"
                f"Or load from a Primer3Plus settings file:\n"
                f"  config.primer3 = Primer3PlusConfig.from_primer3plus_file('settings.txt')"
            )
        if not (self.primer_min_tm <= self.primer_opt_tm <= self.primer_max_tm):
            raise ValueError(
                f"Tm window invalid: min={self.primer_min_tm}, opt={self.primer_opt_tm}, "
                f"max={self.primer_max_tm}. Must satisfy min <= opt <= max."
            )

    def to_primer3_global_args(self) -> Dict[str, Any]:
        """
        Convert this config to the dict format expected by
        primer3.bindings.design_primers(global_args=...).

        Every parameter name matches the Primer3/Primer3Plus naming convention.
        Calls validate_tm() first — raises ValueError if Tm window is not set.
        """
        self.validate_tm()
        args = {
            # ── Main ─────────────────────────────────────────────────────────
            "PRIMER_TASK": self.primer_task,
            "PRIMER_PICK_LEFT_PRIMER": self.primer_pick_left_primer,
            "PRIMER_PICK_RIGHT_PRIMER": self.primer_pick_right_primer,
            "PRIMER_PICK_INTERNAL_OLIGO": self.primer_pick_internal_oligo,
            "PRIMER_NUM_RETURN": self.primer_num_return,
            "PRIMER_PRODUCT_SIZE_RANGE": [
                [lo, hi] for lo, hi in self.primer_product_size_range
            ],

            # ── General Settings ─────────────────────────────────────────────
            "PRIMER_MIN_SIZE": self.primer_min_size,
            "PRIMER_OPT_SIZE": self.primer_opt_size,
            "PRIMER_MAX_SIZE": self.primer_max_size,
            "PRIMER_MIN_TM": self.primer_min_tm,
            "PRIMER_OPT_TM": self.primer_opt_tm,
            "PRIMER_MAX_TM": self.primer_max_tm,
            "PRIMER_MIN_GC": self.primer_min_gc,
            "PRIMER_OPT_GC_PERCENT": self.primer_opt_gc_percent,
            "PRIMER_MAX_GC": self.primer_max_gc,
            "PRIMER_MAX_POLY_X": self.primer_max_poly_x,

            # ── Advanced Settings ────────────────────────────────────────────
            "PRIMER_DNA_CONC": self.primer_dna_conc,
            "PRIMER_SALT_MONOVALENT": self.primer_salt_monovalent,
            "PRIMER_SALT_DIVALENT": self.primer_salt_divalent,
            "PRIMER_DNTP_CONC": self.primer_dntp_conc,
            "PRIMER_SALT_CORRECTIONS": self.primer_salt_corrections,
            "PRIMER_TM_FORMULA": self.primer_tm_formula,
            "PRIMER_ANNEALING_TEMP": self.primer_annealing_temp,
            "PRIMER_MAX_END_STABILITY": self.primer_max_end_stability,
            "PRIMER_MAX_END_GC": self.primer_max_end_gc,
            "PRIMER_GC_CLAMP": self.primer_gc_clamp,
            "PRIMER_MAX_NS_ACCEPTED": self.primer_max_ns_accepted,
            "PRIMER_LIBERAL_BASE": self.primer_liberal_base,
            "PRIMER_LOWERCASE_MASKING": self.primer_lowercase_masking,
            "PRIMER_FIRST_BASE_INDEX": self.primer_first_base_index,
            "PRIMER_MAX_SELF_ANY": self.primer_max_self_any,
            "PRIMER_MAX_SELF_ANY_TH": self.primer_max_self_any_th,
            "PRIMER_MAX_SELF_END": self.primer_max_self_end,
            "PRIMER_MAX_SELF_END_TH": self.primer_max_self_end_th,
            "PRIMER_MAX_HAIRPIN_TH": self.primer_max_hairpin_th,
            "PRIMER_MAX_TEMPLATE_MISPRIMING": self.primer_max_template_mispriming,
            "PRIMER_MAX_TEMPLATE_MISPRIMING_TH": self.primer_max_template_mispriming_th,
            "PRIMER_MAX_LIBRARY_MISPRIMING": self.primer_max_library_mispriming,
            "PRIMER_PAIR_MAX_COMPL_ANY": self.primer_pair_max_compl_any,
            "PRIMER_PAIR_MAX_COMPL_ANY_TH": self.primer_pair_max_compl_any_th,
            "PRIMER_PAIR_MAX_COMPL_END": self.primer_pair_max_compl_end,
            "PRIMER_PAIR_MAX_COMPL_END_TH": self.primer_pair_max_compl_end_th,
            "PRIMER_PAIR_MAX_DIFF_TM": self.primer_pair_max_diff_tm,
            "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING": self.primer_pair_max_library_mispriming,
            "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING": self.primer_pair_max_template_mispriming,
            "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH": self.primer_pair_max_template_mispriming_th,
            "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": self.primer_thermodynamic_oligo_alignment,
            "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT": self.primer_thermodynamic_template_alignment,
            "PRIMER_SECONDARY_STRUCTURE_ALIGNMENT": self.primer_secondary_structure_alignment,
            "PRIMER_DMSO_CONC": self.primer_dmso_conc,
            "PRIMER_DMSO_FACTOR": self.primer_dmso_factor,
            "PRIMER_FORMAMIDE_CONC": self.primer_formamide_conc,
            "PRIMER_PRODUCT_OPT_SIZE": self.primer_product_opt_size,
            "PRIMER_PRODUCT_MAX_TM": self.primer_product_max_tm,
            "PRIMER_PRODUCT_MIN_TM": self.primer_product_min_tm,
            "PRIMER_PRODUCT_OPT_TM": self.primer_product_opt_tm,
            "PRIMER_PICK_ANYWAY": self.primer_pick_anyway,
            "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE": self.primer_min_left_three_prime_distance,
            "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE": self.primer_min_right_three_prime_distance,
            "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION": self.primer_min_3_prime_overlap_of_junction,
            "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION": self.primer_min_5_prime_overlap_of_junction,
            "PRIMER_QUALITY_RANGE_MIN": self.primer_quality_range_min,
            "PRIMER_QUALITY_RANGE_MAX": self.primer_quality_range_max,
            "PRIMER_SEQUENCING_ACCURACY": self.primer_sequencing_accuracy,
            "PRIMER_SEQUENCING_INTERVAL": self.primer_sequencing_interval,
            "PRIMER_SEQUENCING_LEAD": self.primer_sequencing_lead,
            "PRIMER_SEQUENCING_SPACING": self.primer_sequencing_spacing,

            # ── Internal Oligo ───────────────────────────────────────────────
            "PRIMER_INTERNAL_MIN_SIZE": self.primer_internal_min_size,
            "PRIMER_INTERNAL_OPT_SIZE": self.primer_internal_opt_size,
            "PRIMER_INTERNAL_MAX_SIZE": self.primer_internal_max_size,
            "PRIMER_INTERNAL_MIN_TM": self.primer_internal_min_tm,
            "PRIMER_INTERNAL_OPT_TM": self.primer_internal_opt_tm,
            "PRIMER_INTERNAL_MAX_TM": self.primer_internal_max_tm,
            "PRIMER_INTERNAL_MIN_GC": self.primer_internal_min_gc,
            "PRIMER_INTERNAL_OPT_GC_PERCENT": self.primer_internal_opt_gc_percent,
            "PRIMER_INTERNAL_MAX_GC": self.primer_internal_max_gc,
            "PRIMER_INTERNAL_MAX_POLY_X": self.primer_internal_max_poly_x,
            "PRIMER_INTERNAL_MAX_NS_ACCEPTED": self.primer_internal_max_ns_accepted,
            "PRIMER_INTERNAL_MAX_SELF_ANY": self.primer_internal_max_self_any,
            "PRIMER_INTERNAL_MAX_SELF_ANY_TH": self.primer_internal_max_self_any_th,
            "PRIMER_INTERNAL_MAX_SELF_END": self.primer_internal_max_self_end,
            "PRIMER_INTERNAL_MAX_SELF_END_TH": self.primer_internal_max_self_end_th,
            "PRIMER_INTERNAL_MAX_HAIRPIN_TH": self.primer_internal_max_hairpin_th,
            "PRIMER_INTERNAL_MAX_LIBRARY_MISHYB": self.primer_internal_max_library_mishyb,
            "PRIMER_INTERNAL_DNA_CONC": self.primer_internal_dna_conc,
            "PRIMER_INTERNAL_SALT_MONOVALENT": self.primer_internal_salt_monovalent,
            "PRIMER_INTERNAL_SALT_DIVALENT": self.primer_internal_salt_divalent,
            "PRIMER_INTERNAL_DNTP_CONC": self.primer_internal_dntp_conc,
            "PRIMER_INTERNAL_DMSO_CONC": self.primer_internal_dmso_conc,
            "PRIMER_INTERNAL_DMSO_FACTOR": self.primer_internal_dmso_factor,
            "PRIMER_INTERNAL_FORMAMIDE_CONC": self.primer_internal_formamide_conc,
            "PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION": self.primer_internal_min_3_prime_overlap_of_junction,
            "PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION": self.primer_internal_min_5_prime_overlap_of_junction,

            # ── Penalties ────────────────────────────────────────────────────
            "PRIMER_INSIDE_PENALTY": self.primer_inside_penalty,
            "PRIMER_OUTSIDE_PENALTY": self.primer_outside_penalty,
            "PRIMER_WT_TM_GT": self.primer_wt_tm_gt,
            "PRIMER_WT_TM_LT": self.primer_wt_tm_lt,
            "PRIMER_WT_SIZE_GT": self.primer_wt_size_gt,
            "PRIMER_WT_SIZE_LT": self.primer_wt_size_lt,
            "PRIMER_WT_GC_PERCENT_GT": self.primer_wt_gc_percent_gt,
            "PRIMER_WT_GC_PERCENT_LT": self.primer_wt_gc_percent_lt,
            "PRIMER_WT_SELF_ANY": self.primer_wt_self_any,
            "PRIMER_WT_SELF_ANY_TH": self.primer_wt_self_any_th,
            "PRIMER_WT_SELF_END": self.primer_wt_self_end,
            "PRIMER_WT_SELF_END_TH": self.primer_wt_self_end_th,
            "PRIMER_WT_HAIRPIN_TH": self.primer_wt_hairpin_th,
            "PRIMER_WT_END_STABILITY": self.primer_wt_end_stability,
            "PRIMER_WT_END_QUAL": self.primer_wt_end_qual,
            "PRIMER_WT_NUM_NS": self.primer_wt_num_ns,
            "PRIMER_WT_SEQ_QUAL": self.primer_wt_seq_qual,
            "PRIMER_WT_LIBRARY_MISPRIMING": self.primer_wt_library_mispriming,
            "PRIMER_WT_TEMPLATE_MISPRIMING": self.primer_wt_template_mispriming,
            "PRIMER_WT_TEMPLATE_MISPRIMING_TH": self.primer_wt_template_mispriming_th,
            "PRIMER_WT_POS_PENALTY": self.primer_wt_pos_penalty,
            "PRIMER_WT_BOUND_GT": self.primer_wt_bound_gt,
            "PRIMER_WT_BOUND_LT": self.primer_wt_bound_lt,
            "PRIMER_PAIR_WT_DIFF_TM": self.primer_pair_wt_diff_tm,
            "PRIMER_PAIR_WT_COMPL_ANY": self.primer_pair_wt_compl_any,
            "PRIMER_PAIR_WT_COMPL_ANY_TH": self.primer_pair_wt_compl_any_th,
            "PRIMER_PAIR_WT_COMPL_END": self.primer_pair_wt_compl_end,
            "PRIMER_PAIR_WT_COMPL_END_TH": self.primer_pair_wt_compl_end_th,
            "PRIMER_PAIR_WT_PR_PENALTY": self.primer_pair_wt_pr_penalty,
            "PRIMER_PAIR_WT_IO_PENALTY": self.primer_pair_wt_io_penalty,
            "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": self.primer_pair_wt_product_size_gt,
            "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": self.primer_pair_wt_product_size_lt,
            "PRIMER_PAIR_WT_PRODUCT_TM_GT": self.primer_pair_wt_product_tm_gt,
            "PRIMER_PAIR_WT_PRODUCT_TM_LT": self.primer_pair_wt_product_tm_lt,
            "PRIMER_PAIR_WT_LIBRARY_MISPRIMING": self.primer_pair_wt_library_mispriming,
            "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING": self.primer_pair_wt_template_mispriming,
            "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH": self.primer_pair_wt_template_mispriming_th,
            "PRIMER_INTERNAL_WT_TM_GT": self.primer_internal_wt_tm_gt,
            "PRIMER_INTERNAL_WT_TM_LT": self.primer_internal_wt_tm_lt,
            "PRIMER_INTERNAL_WT_SIZE_GT": self.primer_internal_wt_size_gt,
            "PRIMER_INTERNAL_WT_SIZE_LT": self.primer_internal_wt_size_lt,
            "PRIMER_INTERNAL_WT_GC_PERCENT_GT": self.primer_internal_wt_gc_percent_gt,
            "PRIMER_INTERNAL_WT_GC_PERCENT_LT": self.primer_internal_wt_gc_percent_lt,
            "PRIMER_INTERNAL_WT_SELF_ANY": self.primer_internal_wt_self_any,
            "PRIMER_INTERNAL_WT_SELF_ANY_TH": self.primer_internal_wt_self_any_th,
            "PRIMER_INTERNAL_WT_SELF_END": self.primer_internal_wt_self_end,
            "PRIMER_INTERNAL_WT_SELF_END_TH": self.primer_internal_wt_self_end_th,
            "PRIMER_INTERNAL_WT_HAIRPIN_TH": self.primer_internal_wt_hairpin_th,
            "PRIMER_INTERNAL_WT_NUM_NS": self.primer_internal_wt_num_ns,
            "PRIMER_INTERNAL_WT_SEQ_QUAL": self.primer_internal_wt_seq_qual,
            "PRIMER_INTERNAL_WT_END_QUAL": self.primer_internal_wt_end_qual,
            "PRIMER_INTERNAL_WT_LIBRARY_MISHYB": self.primer_internal_wt_library_mishyb,
            "PRIMER_INTERNAL_WT_BOUND_GT": self.primer_internal_wt_bound_gt,
            "PRIMER_INTERNAL_WT_BOUND_LT": self.primer_internal_wt_bound_lt,

            # ── Advanced Seq. ────────────────────────────────────────────────
            "PRIMER_MUST_MATCH_FIVE_PRIME": self.primer_must_match_five_prime,
            "PRIMER_MUST_MATCH_THREE_PRIME": self.primer_must_match_three_prime,
            "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME": self.primer_internal_must_match_five_prime,
            "PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME": self.primer_internal_must_match_three_prime,
            "PRIMER_MIN_QUALITY": self.primer_min_quality,
            "PRIMER_MIN_END_QUALITY": self.primer_min_end_quality,
            "PRIMER_INTERNAL_MIN_QUALITY": self.primer_internal_min_quality,
            "PRIMER_MIN_BOUND": self.primer_min_bound,
            "PRIMER_MAX_BOUND": self.primer_max_bound,
            "PRIMER_OPT_BOUND": self.primer_opt_bound,
            "PRIMER_INTERNAL_MIN_BOUND": self.primer_internal_min_bound,
            "PRIMER_INTERNAL_MAX_BOUND": self.primer_internal_max_bound,
            "PRIMER_INTERNAL_OPT_BOUND": self.primer_internal_opt_bound,

            # ── Explain flag (always return explanation) ─────────────────────
            "PRIMER_EXPLAIN_FLAG": 1,
        }

        # Remove empty must-match strings (primer3 requires exactly 5 chars
        # or the key must be absent entirely)
        for k in ["PRIMER_MUST_MATCH_FIVE_PRIME", "PRIMER_MUST_MATCH_THREE_PRIME",
                  "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME",
                  "PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME"]:
            if k in args and not args[k]:
                del args[k]

        # Remove mispriming library references if the library file is not
        # available (primer3 will error if the file doesn't exist)
        for k in ["PRIMER_MISPRIMING_LIBRARY", "PRIMER_INTERNAL_MISHYB_LIBRARY"]:
            if k in args:
                del args[k]

        return args

    @classmethod
    def from_primer3plus_file(cls, filepath: str) -> "Primer3PlusConfig":
        """
        Load settings from a Primer3Plus settings file (P3_FILE_TYPE=settings).

        This allows loading the exact same settings you use in the
        Primer3Plus web interface.
        """
        config = cls()
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("P3_"):
                    continue
                if "=" not in line:
                    continue
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip()

                # Map Primer3Plus keys to dataclass fields
                field_map = {
                    "PRIMER_TASK": "primer_task",
                    "PRIMER_PICK_LEFT_PRIMER": "primer_pick_left_primer",
                    "PRIMER_PICK_RIGHT_PRIMER": "primer_pick_right_primer",
                    "PRIMER_PICK_INTERNAL_OLIGO": "primer_pick_internal_oligo",
                    "PRIMER_NUM_RETURN": "primer_num_return",
                    "PRIMER_MIN_SIZE": "primer_min_size",
                    "PRIMER_OPT_SIZE": "primer_opt_size",
                    "PRIMER_MAX_SIZE": "primer_max_size",
                    "PRIMER_MIN_TM": "primer_min_tm",
                    "PRIMER_OPT_TM": "primer_opt_tm",
                    "PRIMER_MAX_TM": "primer_max_tm",
                    "PRIMER_MIN_GC": "primer_min_gc",
                    "PRIMER_OPT_GC_PERCENT": "primer_opt_gc_percent",
                    "PRIMER_MAX_GC": "primer_max_gc",
                    "PRIMER_MAX_POLY_X": "primer_max_poly_x",
                    "PRIMER_DNA_CONC": "primer_dna_conc",
                    "PRIMER_SALT_MONOVALENT": "primer_salt_monovalent",
                    "PRIMER_SALT_DIVALENT": "primer_salt_divalent",
                    "PRIMER_DNTP_CONC": "primer_dntp_conc",
                    "PRIMER_SALT_CORRECTIONS": "primer_salt_corrections",
                    "PRIMER_TM_FORMULA": "primer_tm_formula",
                    "PRIMER_ANNEALING_TEMP": "primer_annealing_temp",
                    "PRIMER_MAX_END_STABILITY": "primer_max_end_stability",
                    "PRIMER_MAX_END_GC": "primer_max_end_gc",
                    "PRIMER_GC_CLAMP": "primer_gc_clamp",
                    "PRIMER_MAX_NS_ACCEPTED": "primer_max_ns_accepted",
                    "PRIMER_LIBERAL_BASE": "primer_liberal_base",
                    "PRIMER_LOWERCASE_MASKING": "primer_lowercase_masking",
                    "PRIMER_FIRST_BASE_INDEX": "primer_first_base_index",
                    "PRIMER_MAX_SELF_ANY": "primer_max_self_any",
                    "PRIMER_MAX_SELF_ANY_TH": "primer_max_self_any_th",
                    "PRIMER_MAX_SELF_END": "primer_max_self_end",
                    "PRIMER_MAX_SELF_END_TH": "primer_max_self_end_th",
                    "PRIMER_MAX_HAIRPIN_TH": "primer_max_hairpin_th",
                    "PRIMER_MAX_TEMPLATE_MISPRIMING": "primer_max_template_mispriming",
                    "PRIMER_MAX_TEMPLATE_MISPRIMING_TH": "primer_max_template_mispriming_th",
                    "PRIMER_MAX_LIBRARY_MISPRIMING": "primer_max_library_mispriming",
                    "PRIMER_PAIR_MAX_COMPL_ANY": "primer_pair_max_compl_any",
                    "PRIMER_PAIR_MAX_COMPL_ANY_TH": "primer_pair_max_compl_any_th",
                    "PRIMER_PAIR_MAX_COMPL_END": "primer_pair_max_compl_end",
                    "PRIMER_PAIR_MAX_COMPL_END_TH": "primer_pair_max_compl_end_th",
                    "PRIMER_PAIR_MAX_DIFF_TM": "primer_pair_max_diff_tm",
                    "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING": "primer_pair_max_library_mispriming",
                    "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING": "primer_pair_max_template_mispriming",
                    "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH": "primer_pair_max_template_mispriming_th",
                    "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": "primer_thermodynamic_oligo_alignment",
                    "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT": "primer_thermodynamic_template_alignment",
                    "PRIMER_SECONDARY_STRUCTURE_ALIGNMENT": "primer_secondary_structure_alignment",
                    "PRIMER_DMSO_CONC": "primer_dmso_conc",
                    "PRIMER_DMSO_FACTOR": "primer_dmso_factor",
                    "PRIMER_FORMAMIDE_CONC": "primer_formamide_conc",
                    "PRIMER_PRODUCT_OPT_SIZE": "primer_product_opt_size",
                    "PRIMER_PRODUCT_MAX_TM": "primer_product_max_tm",
                    "PRIMER_PRODUCT_MIN_TM": "primer_product_min_tm",
                    "PRIMER_PRODUCT_OPT_TM": "primer_product_opt_tm",
                    "PRIMER_PICK_ANYWAY": "primer_pick_anyway",
                    "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE": "primer_min_left_three_prime_distance",
                    "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE": "primer_min_right_three_prime_distance",
                    "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION": "primer_min_3_prime_overlap_of_junction",
                    "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION": "primer_min_5_prime_overlap_of_junction",
                    "PRIMER_QUALITY_RANGE_MIN": "primer_quality_range_min",
                    "PRIMER_QUALITY_RANGE_MAX": "primer_quality_range_max",
                    "PRIMER_SEQUENCING_ACCURACY": "primer_sequencing_accuracy",
                    "PRIMER_SEQUENCING_INTERVAL": "primer_sequencing_interval",
                    "PRIMER_SEQUENCING_LEAD": "primer_sequencing_lead",
                    "PRIMER_SEQUENCING_SPACING": "primer_sequencing_spacing",
                    "PRIMER_INSIDE_PENALTY": "primer_inside_penalty",
                    "PRIMER_OUTSIDE_PENALTY": "primer_outside_penalty",
                    "PRIMER_WT_TM_GT": "primer_wt_tm_gt",
                    "PRIMER_WT_TM_LT": "primer_wt_tm_lt",
                    "PRIMER_WT_SIZE_GT": "primer_wt_size_gt",
                    "PRIMER_WT_SIZE_LT": "primer_wt_size_lt",
                    "PRIMER_WT_GC_PERCENT_GT": "primer_wt_gc_percent_gt",
                    "PRIMER_WT_GC_PERCENT_LT": "primer_wt_gc_percent_lt",
                    "PRIMER_WT_SELF_ANY": "primer_wt_self_any",
                    "PRIMER_WT_SELF_ANY_TH": "primer_wt_self_any_th",
                    "PRIMER_WT_SELF_END": "primer_wt_self_end",
                    "PRIMER_WT_SELF_END_TH": "primer_wt_self_end_th",
                    "PRIMER_WT_HAIRPIN_TH": "primer_wt_hairpin_th",
                    "PRIMER_WT_END_STABILITY": "primer_wt_end_stability",
                    "PRIMER_WT_END_QUAL": "primer_wt_end_qual",
                    "PRIMER_WT_NUM_NS": "primer_wt_num_ns",
                    "PRIMER_WT_SEQ_QUAL": "primer_wt_seq_qual",
                    "PRIMER_WT_LIBRARY_MISPRIMING": "primer_wt_library_mispriming",
                    "PRIMER_WT_TEMPLATE_MISPRIMING": "primer_wt_template_mispriming",
                    "PRIMER_WT_TEMPLATE_MISPRIMING_TH": "primer_wt_template_mispriming_th",
                    "PRIMER_WT_POS_PENALTY": "primer_wt_pos_penalty",
                    "PRIMER_WT_BOUND_GT": "primer_wt_bound_gt",
                    "PRIMER_WT_BOUND_LT": "primer_wt_bound_lt",
                    "PRIMER_PAIR_WT_DIFF_TM": "primer_pair_wt_diff_tm",
                    "PRIMER_PAIR_WT_COMPL_ANY": "primer_pair_wt_compl_any",
                    "PRIMER_PAIR_WT_COMPL_ANY_TH": "primer_pair_wt_compl_any_th",
                    "PRIMER_PAIR_WT_COMPL_END": "primer_pair_wt_compl_end",
                    "PRIMER_PAIR_WT_COMPL_END_TH": "primer_pair_wt_compl_end_th",
                    "PRIMER_PAIR_WT_PR_PENALTY": "primer_pair_wt_pr_penalty",
                    "PRIMER_PAIR_WT_IO_PENALTY": "primer_pair_wt_io_penalty",
                    "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": "primer_pair_wt_product_size_gt",
                    "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": "primer_pair_wt_product_size_lt",
                    "PRIMER_PAIR_WT_PRODUCT_TM_GT": "primer_pair_wt_product_tm_gt",
                    "PRIMER_PAIR_WT_PRODUCT_TM_LT": "primer_pair_wt_product_tm_lt",
                    "PRIMER_PAIR_WT_LIBRARY_MISPRIMING": "primer_pair_wt_library_mispriming",
                    "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING": "primer_pair_wt_template_mispriming",
                    "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH": "primer_pair_wt_template_mispriming_th",
                    "PRIMER_INTERNAL_MIN_SIZE": "primer_internal_min_size",
                    "PRIMER_INTERNAL_OPT_SIZE": "primer_internal_opt_size",
                    "PRIMER_INTERNAL_MAX_SIZE": "primer_internal_max_size",
                    "PRIMER_INTERNAL_MIN_TM": "primer_internal_min_tm",
                    "PRIMER_INTERNAL_OPT_TM": "primer_internal_opt_tm",
                    "PRIMER_INTERNAL_MAX_TM": "primer_internal_max_tm",
                    "PRIMER_INTERNAL_MIN_GC": "primer_internal_min_gc",
                    "PRIMER_INTERNAL_OPT_GC_PERCENT": "primer_internal_opt_gc_percent",
                    "PRIMER_INTERNAL_MAX_GC": "primer_internal_max_gc",
                    "PRIMER_INTERNAL_MAX_POLY_X": "primer_internal_max_poly_x",
                    "PRIMER_INTERNAL_MAX_NS_ACCEPTED": "primer_internal_max_ns_accepted",
                    "PRIMER_INTERNAL_MAX_SELF_ANY": "primer_internal_max_self_any",
                    "PRIMER_INTERNAL_MAX_SELF_ANY_TH": "primer_internal_max_self_any_th",
                    "PRIMER_INTERNAL_MAX_SELF_END": "primer_internal_max_self_end",
                    "PRIMER_INTERNAL_MAX_SELF_END_TH": "primer_internal_max_self_end_th",
                    "PRIMER_INTERNAL_MAX_HAIRPIN_TH": "primer_internal_max_hairpin_th",
                    "PRIMER_INTERNAL_MAX_LIBRARY_MISHYB": "primer_internal_max_library_mishyb",
                    "PRIMER_INTERNAL_DNA_CONC": "primer_internal_dna_conc",
                    "PRIMER_INTERNAL_SALT_MONOVALENT": "primer_internal_salt_monovalent",
                    "PRIMER_INTERNAL_SALT_DIVALENT": "primer_internal_salt_divalent",
                    "PRIMER_INTERNAL_DNTP_CONC": "primer_internal_dntp_conc",
                    "PRIMER_INTERNAL_DMSO_CONC": "primer_internal_dmso_conc",
                    "PRIMER_INTERNAL_DMSO_FACTOR": "primer_internal_dmso_factor",
                    "PRIMER_INTERNAL_FORMAMIDE_CONC": "primer_internal_formamide_conc",
                    "PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION": "primer_internal_min_3_prime_overlap_of_junction",
                    "PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION": "primer_internal_min_5_prime_overlap_of_junction",
                    "PRIMER_INTERNAL_WT_TM_GT": "primer_internal_wt_tm_gt",
                    "PRIMER_INTERNAL_WT_TM_LT": "primer_internal_wt_tm_lt",
                    "PRIMER_INTERNAL_WT_SIZE_GT": "primer_internal_wt_size_gt",
                    "PRIMER_INTERNAL_WT_SIZE_LT": "primer_internal_wt_size_lt",
                    "PRIMER_INTERNAL_WT_GC_PERCENT_GT": "primer_internal_wt_gc_percent_gt",
                    "PRIMER_INTERNAL_WT_GC_PERCENT_LT": "primer_internal_wt_gc_percent_lt",
                    "PRIMER_INTERNAL_WT_SELF_ANY": "primer_internal_wt_self_any",
                    "PRIMER_INTERNAL_WT_SELF_ANY_TH": "primer_internal_wt_self_any_th",
                    "PRIMER_INTERNAL_WT_SELF_END": "primer_internal_wt_self_end",
                    "PRIMER_INTERNAL_WT_SELF_END_TH": "primer_internal_wt_self_end_th",
                    "PRIMER_INTERNAL_WT_HAIRPIN_TH": "primer_internal_wt_hairpin_th",
                    "PRIMER_INTERNAL_WT_NUM_NS": "primer_internal_wt_num_ns",
                    "PRIMER_INTERNAL_WT_SEQ_QUAL": "primer_internal_wt_seq_qual",
                    "PRIMER_INTERNAL_WT_END_QUAL": "primer_internal_wt_end_qual",
                    "PRIMER_INTERNAL_WT_LIBRARY_MISHYB": "primer_internal_wt_library_mishyb",
                    "PRIMER_INTERNAL_WT_BOUND_GT": "primer_internal_wt_bound_gt",
                    "PRIMER_INTERNAL_WT_BOUND_LT": "primer_internal_wt_bound_lt",
                    "PRIMER_MUST_MATCH_FIVE_PRIME": "primer_must_match_five_prime",
                    "PRIMER_MUST_MATCH_THREE_PRIME": "primer_must_match_three_prime",
                    "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME": "primer_internal_must_match_five_prime",
                    "PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME": "primer_internal_must_match_three_prime",
                    "PRIMER_MIN_QUALITY": "primer_min_quality",
                    "PRIMER_MIN_END_QUALITY": "primer_min_end_quality",
                    "PRIMER_INTERNAL_MIN_QUALITY": "primer_internal_min_quality",
                    "PRIMER_MIN_BOUND": "primer_min_bound",
                    "PRIMER_MAX_BOUND": "primer_max_bound",
                    "PRIMER_OPT_BOUND": "primer_opt_bound",
                    "PRIMER_INTERNAL_MIN_BOUND": "primer_internal_min_bound",
                    "PRIMER_INTERNAL_MAX_BOUND": "primer_internal_max_bound",
                    "PRIMER_INTERNAL_OPT_BOUND": "primer_internal_opt_bound",
                }

                if key in field_map:
                    attr = field_map[key]
                    # Type conversion using dataclass field type annotations
                    # This handles None defaults (e.g. required Tm fields) correctly
                    from dataclasses import fields as _dc_fields
                    _field_types = {f.name: f.type for f in _dc_fields(Primer3PlusConfig)}
                    _ftype = _field_types.get(attr)
                    # Resolve string type annotations (Optional[float] etc.)
                    _type_str = str(_ftype) if _ftype else ""
                    try:
                        if "int" in _type_str and "float" not in _type_str:
                            setattr(config, attr, int(float(value)))
                        elif "float" in _type_str or "Optional[float]" in _type_str:
                            setattr(config, attr, float(value))
                        elif "int" in _type_str:
                            setattr(config, attr, int(float(value)))
                        else:
                            # Fallback: try float, then int, then str
                            try:
                                if "." in value or "e" in value.lower():
                                    setattr(config, attr, float(value))
                                else:
                                    setattr(config, attr, int(float(value)))
                            except (ValueError, TypeError):
                                setattr(config, attr, value)
                    except (ValueError, TypeError):
                        pass

                # Special handling for product size range
                if key == "PRIMER_PRODUCT_SIZE_RANGE":
                    # Parse "60-150" or "60-150 200-300"
                    ranges = []
                    for part in value.split():
                        if "-" in part:
                            lo, hi = part.split("-")
                            ranges.append((int(lo), int(hi)))
                    if ranges:
                        config.primer_product_size_range = ranges

        # Mark salt preset as custom since we loaded from file
        config.salt_preset = "custom"
        return config


# ─────────────────────────────────────────────────────────────────────────────
# Pipeline-level configuration
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class PipelineConfig:
    """Top-level pipeline configuration."""

    # Target
    target_cell_type: str = "Blood-Monocytes"
    background_cell_types: Optional[List[str]] = None

    # Input
    fasta_sequence: Optional[str] = None
    region: Optional[str] = None
    genome: str = "hg19"

    # Atlas paths
    atlas_dir: Optional[str] = None
    blocks_file: Optional[str] = None
    groups_csv: Optional[str] = None
    markers_output_dir: str = "markers_output"

    # DMR selection thresholds
    delta_means: float = 0.3
    delta_quants: float = 0.3
    pval: float = 0.05
    min_cpg: int = 5              # hard minimum CpGs per DMR block
    preferred_min_cpg: int = 7    # preferred minimum; blocks with >=7 are flagged as preferred
    only_hypo: bool = True
    top_n_markers: int = 1000
    sort_by: str = "delta_means"
    threads: int = 8
    max_primers_per_block: int = 200  # cap primer pairs returned per block per template (was 10)
    min_cpg_per_primer: int = 2       # min CpGs in at least one primer of the pair
    min_cpg_pair_total: int = 4       # min total CpGs across both primers

    # CpG placement mode: "strict_terminal" (last 5 nt) or "three_prime_half"
    # (last 50% of primer length — used for MSP-style discrimination)
    cpg_placement_mode: str = "strict_terminal"

    # Sub-configs
    primer3: Primer3PlusConfig = field(default_factory=Primer3PlusConfig)

    # Output
    output_dir: str = "pipeline_output"
    report_title: str = "Cell-Specific Epigenetic dPCR Assay Design Report"

    # dbSNP
    dbsnp_vcf: Optional[str] = None

    # Segmental duplication & repeat overlap (UCSC tracks)
    segdup_bed: Optional[str] = None
    repeat_bed: Optional[str] = None

    def __post_init__(self):
        if self.background_cell_types is None:
            self.background_cell_types = []
