#!/usr/bin/env python3
"""data_dictionary.py - Build a data dictionary DataFrame for the master atlas xlsx.

Explains every column in every sheet with maximum detail.
"""
import pandas as pd


def build_data_dictionary():
    """Build and return a DataFrame with detailed explanations of every column."""
    DD = []
    
    def add(sheet, no, col, dtype, desc, examples=''):
        DD.append({
            'Sheet': sheet, 'No': no, 'Column': col, 'Type': dtype,
            'Description': desc, 'Example_Values': examples
        })
    
    # === Sheet descriptions (No = -1 marks sheet-level rows) ===
    add('ALL SHEETS', -1, '---', '---',
        'This file contains 6 sheets:\n'
        '1. All_Primers: Every primer pair designed (391,570 rows x 65 cols). One row per primer-template combination.\n'
        '2. Best_Per_DMR: The best-scoring primer per DMR (46,977 rows x 60 cols). Subset selected by the design pipeline.\n'
        '3. Selected_Panel: The 976 primers that passed ePCR screening and were selected for the final panel (65 cols). Subset of All_Primers where selected_for_panel=1.\n'
        '4. DMR_Blocks: DMR-level summary, one row per differentially methylated region (46,977 rows x 18 cols).\n'
        '5. Summary: Per-cell-type aggregate statistics (161 rows x 10 cols).\n'
        '6. Data_Dictionary: This sheet.\n'
        'Data source: Loyfer et al. Nature 2023 (GSE186458), hg19/GRCh37 assembly, 82 cell types (14 blood + 68 non-blood).',
        '')
    
    # === All_Primers / Best_Per_DMR / Selected_Panel shared columns ===
    S = 'All_Primers, Best_Per_DMR, Selected_Panel'
    
    add(S, 0, 'cell_type', 'string',
        'The immune or tissue cell type that this DMR discriminates. 82 cell types total from Loyfer et al. Nature 2023 (GSE186458). '
        'Includes 14 blood cell types (e.g., Blood-T-CD4, Blood-B, Blood-NK) and 68 non-blood tissue types '
        '(e.g., Adipocytes, Neuron, Liver-Endothelium, Pancreas-Beta). Each DMR is specific to one cell type.',
        'Adipocytes, Blood-T-CD4, Neuron, Pancreas-Beta')
    
    add(S, 1, 'source', 'string',
        'Pipeline version that generated this DMR and primer design entry.',
        'v2.2.8')
    
    add(S, 2, 'rank', 'integer',
        'DMR rank within the cell type, ordered by delta_means (methylation difference between target and background) descending. '
        'Rank 1 = the DMR with the highest discrimination power for that cell type. '
        'Used for prioritizing which DMRs to screen first (top 200 per cell type are screened).',
        '1, 2, 3, ...')
    
    add(S, 3, 'chrom', 'string',
        'Chromosome name in hg19/GRCh37 coordinates. All genomic coordinates in this atlas use the hg19 assembly.',
        'chr1, chr17, chrX')
    
    add(S, 4, 'dmr_start', 'integer',
        'Start position (1-based) of the DMR (differentially methylated region) on the chromosome. '
        'The DMR is the genomic region where methylation differs significantly between the target cell type and all other cell types.',
        '4315368')
    
    add(S, 5, 'dmr_end', 'integer',
        'End position (1-based, inclusive) of the DMR on the chromosome. DMR length = dmr_end - dmr_start + 1.',
        '4315523')
    
    add(S, 6, 'gene', 'string',
        'Name of the gene nearest to or overlapping the DMR. Used for biological interpretation of the methylation marker. '
        'May be blank for intergenic DMRs with no nearby gene.',
        'SPNS3, LHPP, SNAI1')
    
    add(S, 7, 'annotation', 'string',
        'Genomic feature annotation describing the DMR relationship to the nearest gene. '
        'Intergenic = not within any gene body. intron = within an intron. exon = within an exon. '
        'promoter = near transcription start site.',
        'Intergenic, intron, exon')
    
    add(S, 8, 'delta_means', 'float',
        'Difference in mean methylation level (0-1) between the target cell type and the mean of all background cell types. '
        'Higher values indicate stronger discrimination power. This is the PRIMARY metric for ranking DMRs within each cell type. '
        'Calculated from WGBS data (Loyfer et al. Nature 2023, GSE186458). Range in this atlas: 0.50-0.95.',
        '0.779, 0.742, 0.728')
    
    add(S, 9, 'dms_score', 'float',
        'Differential Methylation Suitability score (0-1). Composite metric evaluating how '
        'suitable a DMR block is for bisulfite PCR assay design. Computed from per-CpG WGBS '
        'methylation data across target and background cell types.\n\n'
        'Formula: DMS_score = mean(A, B, C, D) * E\n\n'
        'A (score_A_tg_near0) = 1 - max(per-CpG mean methylation in target samples).\n'
        '  Measures how completely unmethylated the target is at every CpG. Uses max() because '
        'a single CpG with high target methylation would cause false positives.\n'
        'B (score_B_bg_near1) = min(per-CpG mean methylation in background samples).\n'
        '  Measures how completely methylated the background is at every CpG. Uses min() because '
        'a single CpG with low background methylation would cause false positives.\n'
        'C (score_C_tg_consist) = max(0, 1 - 2 * mean(per-CpG std across target samples)).\n'
        '  Measures donor-to-donor consistency in the target. Factor of 2 scales so std=0.5 gives 0.\n'
        'D (score_D_bg_consist) = max(0, 1 - 2 * mean(per-CpG std across background samples)).\n'
        '  Measures consistency across background cell types and donors.\n'
        'E (score_E_coverage) = mean(fraction CpGs with data in >=50% target samples, '
        'fraction CpGs with data in >=50% background samples).\n'
        '  Coverage multiplier. Scores based on sparse data are unreliable even if available data looks perfect.\n\n'
        'E is a multiplier (not additive) because data reliability gates all other components. '
        'Threshold: DMS_score >= 0.6 = good candidate for primer design. '
        'Data source: Loyfer et al. Nature 2023 (GSE186458), hg19.',
        '0.8147, 0.7952, 0.6034')
    
    add(S, 10, 'score_A_tg_near0', 'float',
        'DMS sub-score A: target near-zero methylation. Formula: A = 1 - max(per-CpG mean methylation in target samples). '
        'For each CpG in the block, mean methylation is computed across all target samples (nanmean). '
        'The maximum of these per-CpG means is subtracted from 1. '
        'Range: 1.0 = all CpGs fully unmethylated in target. 0.0 = at least one CpG 100% methylated in target. '
        'Uses max() not mean() because a single CpG with high target methylation would cause false positives '
        'if a primer spans it.',
        '1.0, 0.95, 0.8')
    
    add(S, 11, 'score_B_bg_near1', 'float',
        'DMS sub-score B: background near-one methylation. Formula: B = min(per-CpG mean methylation in background samples). '
        'For each CpG, mean methylation is computed across all background samples (nanmean). '
        'The minimum of these per-CpG means is taken. '
        'Range: 1.0 = all CpGs fully methylated in all background cell types. 0.0 = at least one CpG 0% methylated in background. '
        'Uses min() because a single CpG with low background methylation would cause false positives.',
        '1.0, 0.9, 0.7')
    
    add(S, 12, 'score_C_tg_consist', 'float',
        'DMS sub-score C: target consistency. Formula: C = max(0, 1 - 2 * mean(per-CpG std across target samples)). '
        'For each CpG, std of methylation is computed across target samples (nanstd). '
        'The mean of these per-CpG stds is scaled by factor 2 and subtracted from 1. '
        'The factor of 2 scales so that max possible std (0.5, when half samples are 0 and half are 1) gives score 0. '
        'Range: 1.0 = all donors identical. 0.0 = maximum inconsistency. '
        'Inconsistent markers are unreliable as clinical assays.',
        '1.0, 0.967, 0.8')
    
    add(S, 13, 'score_D_bg_consist', 'float',
        'DMS sub-score D: background consistency. Formula: D = max(0, 1 - 2 * mean(per-CpG std across background samples)). '
        'Identical to C but using the background matrix. Measures consistency across all background cell types and donors.',
        '1.0, 0.9, 0.7')
    
    add(S, 14, 'score_E_coverage', 'float',
        'DMS sub-score E: data coverage. Formula: E = mean(fraction CpGs with data in >=50% target samples, '
        'fraction CpGs with data in >=50% background samples). '
        'For each CpG, the fraction of samples with non-NaN data is computed. '
        'A CpG is considered covered if >=50% of samples have data at that position. '
        'E is the average of the target and background coverage fractions. '
        'E is a multiplier in the DMS formula (not additive) because scores based on sparse data '
        'are unreliable even if the available data looks perfect. '
        'Range: 1.0 = all CpGs have data in all samples. 0.0 = no data.',
        '1.0, 0.85, 0.5')
    
    add(S, 15, 'target_mean_meth', 'float',
        'Mean methylation level (0-1) of all CpG sites within the DMR, measured in the target cell type from WGBS data. '
        'For U-direction DMRs: low (hypomethylated, typically <0.15). For M-direction DMRs: high (hypermethylated, typically >0.7).',
        '0.102 (U-DMR), 0.752 (M-DMR)')
    
    add(S, 16, 'background_mean_meth', 'float',
        'Mean methylation level (0-1) of all CpG sites within the DMR, averaged across all non-target (background) cell types from WGBS data. '
        'For U-direction DMRs: high (typically >0.8). For M-direction DMRs: low (typically <0.2).',
        '0.881 (U-DMR), 0.141 (M-DMR)')
    
    add(S, 17, 'num_cpgs', 'integer',
        'Number of CpG dinucleotide sites within the DMR block. More CpGs provide more methylation measurement points '
        'and potentially better discrimination between cell types.',
        '8, 5, 3, 10')
    
    add(S, 18, 'block_len', 'integer',
        'Length of the DMR block in base pairs. Longer blocks may contain more CpGs but require longer amplicons '
        'or multiple primer pairs to cover.',
        '155, 46, 127, 171')
    
    add(S, 19, 'assay_id', 'string',
        'Unique identifier for each primer pair (assay). Format: {cell_type}_{rank:04d}_{M or U}_{seq:03d}. '
        'M = DMR hypermethylated in target cell type. U = DMR hypomethylated in target cell type. '
        'seq (001, 002, ...) numbers multiple primer designs targeting the same DMR. '
        'Parse with rsplit("_", 2) to handle multi-word cell types (e.g., Bone_marrow-Erythrocyte_progenitors).',
        'Adipocytes_0001_U_001, Neuron_0027_M_001')
    
    add(S, 20, 'seq_id', 'string',
        'DMR-level identifier. Format: {cell_type}_{rank:04d}_{M or U}. '
        'Multiple primers (assay_ids) share the same seq_id when they target the same DMR with different primer designs. '
        'The M/U suffix indicates DMR methylation state in the target cell type, NOT which genome state the product should be on.',
        'Adipocytes_0001_U, Neuron_0027_M')
    
    add(S, 21, 'template_used', 'string',
        'Four-character code specifying the primer strand orientation and methylation specificity. '
        'First letter: S = Sense (Watson/forward strand), A = Antisense (Crick/reverse strand). '
        'Second letter: M = Methylated (primer designed for CT_M converted genome, CpG cytosines preserved), '
        'U = Unmethylated (primer designed for CT converted genome, all cytosines converted to thymine). '
        'SM = Sense Methylated, AM = Antisense Methylated, SU = Sense Unmethylated, AU = Antisense Unmethylated. '
        'The template determines which genome state the PCR product should be on: SM/AM -> methylated (CT_M), SU/AU -> unmethylated (CT). '
        'This is confirmed by the ePCR_verdict column. All four templates are generated as candidates for each assay_id.',
        'SM, AM, SU, AU')
    
    add(S, 22, 'assay', 'string',
        'Assay configuration descriptor. May combine template and direction information or describe the assay design approach.',
        '')
    
    add(S, 23, 'left_primer', 'string',
        'Forward primer sequence (5 to 3), designed to match the bisulfite-converted genome (NOT the reference genome). '
        'At CpG positions: C for methylated-specific templates (SM/AM), T for unmethylated-specific templates (SU/AU). '
        'Non-CpG cytosines are converted to T in all bisulfite-converted primers. Length typically 18-25 bp.',
        'ATAGGCGTGAGTTATTGC')
    
    add(S, 24, 'right_primer', 'string',
        'Reverse primer sequence (5 to 3), designed to match the bisulfite-converted antisense strand. '
        'The reverse complement of this sequence is used for ePCR screening. Length typically 18-25 bp.',
        'CCATCATCTACTACGCATAC')
    
    add(S, 25, 'left_primer_cpg_display', 'string',
        'Annotation of the left primer showing CpG positions in the genome at the primer binding region. '
        'Uses "yg" markers where y marks the cytosine position and g marks the guanine position of each CpG dinucleotide. '
        'Non-CpG positions show the actual primer base. This is a GENOME property (same for SM and SU primers covering '
        'the same genomic positions), not a primer sequence property.',
        'GAAAATAGTTTAGGTTGTAGATG (no CpG), CygygCATTATTTygCTCATTC (3 CpG sites)')
    
    add(S, 26, 'right_primer_cpg_display', 'string',
        'Same as left_primer_cpg_display but for the right primer binding region.',
        'TCTygCygAygAACCACAAC (3 CpG sites)')
    
    add(S, 27, 'left_tm_C', 'float',
        'Melting temperature (Tm) of the left primer in degrees Celsius, calculated using nearest-neighbor thermodynamic parameters. '
        'Optimal range for bisulfite PCR: 55-65 C. Both primers in a pair should have similar Tm (within 2-3 C).',
        '58.3, 62.1')
    
    add(S, 28, 'right_tm_C', 'float',
        'Melting temperature of the right primer in degrees Celsius.',
        '57.9, 61.5')
    
    add(S, 29, 'product_size_bp', 'integer',
        'Expected PCR product size in base pairs when the primer pair amplifies its intended target. '
        'Must be <= 1000 bp for ePCR screening (MAX_PCR_LEN parameter). Typical range 40-500 bp for bisulfite PCR.',
        '155, 80, 97')
    
    add(S, 30, 'c_total', 'integer',
        'Total number of CpG sites in both primer binding regions combined (left_c_total + right_c_total). '
        'More CpG sites in the primers = better methylation discrimination, because each CpG acts as a '
        'methylation-specific discrimination site. Range in this atlas: 3-13. '
        'This is the secondary ranking criterion for primer selection (after delta_means).',
        '4, 3, 5, 7')
    
    add(S, 31, 'left_c_total', 'integer',
        'Number of CpG sites in the left primer binding region. This is a genome property (counts CpG sites in the '
        'genome at the primer binding location), not a primer sequence property. Same for SM and SU primers covering '
        'the same genomic positions.',
        '0, 1, 2, 3')
    
    add(S, 32, 'right_c_total', 'integer',
        'Number of CpG sites in the right primer binding region. Same genome property as left_c_total.',
        '1, 2, 3, 4')
    
    add(S, 33, 'c_total_tail', 'integer',
        'Number of CpG sites in the 3-prime tail regions of both primers combined. CpG sites near the 3 end are most '
        'critical for methylation-specific priming: a mismatch at the 3-terminal CpG prevents polymerase extension, '
        'providing methylation specificity.',
        '0, 1, 2')
    
    add(S, 34, 'left_c_tail', 'integer',
        'Number of CpG sites in the 3-prime tail region of the left primer.',
        '0, 1')
    
    add(S, 35, 'right_c_tail', 'integer',
        'Number of CpG sites in the 3-prime tail region of the right primer.',
        '0, 1, 2')
    
    add(S, 36, 'penalty', 'float',
        'Overall primer design penalty score. Combines Tm mismatch between primers, GC content deviation, '
        'self-complementarity, 3-prime end stability, and CpG placement. Lower = better primer design. '
        'Used for selecting the best primer per DMR during design.',
        '0.5, 1.2, 2.8')
    
    add(S, 37, 'left_gc_percent', 'float',
        'GC content of the left primer as a percentage (0-100). Bisulfite conversion reduces GC content because '
        'cytosines are converted to thymines. Typical range 20-60% for bisulfite primers.',
        '42.1, 50.0, 33.3')
    
    add(S, 38, 'right_gc_percent', 'float',
        'GC content of the right primer as a percentage (0-100).',
        '38.9, 47.6')
    
    add(S, 39, 'left_structure_mfe_kcal_mol', 'float',
        'Minimum free energy (MFE) of the left primer most stable secondary structure (hairpin, self-fold) in kcal/mol. '
        'More negative = more stable structure = potentially worse primer (may self-fold instead of binding template). '
        'Good primers typically have MFE > -5 kcal/mol.',
        '-2.1, -5.3, -7.8')
    
    add(S, 40, 'right_structure_mfe_kcal_mol', 'float',
        'MFE of the right primer secondary structure in kcal/mol. Same interpretation as left_structure_mfe_kcal_mol.',
        '-1.8, -4.2, -6.5')
    
    add(S, 41, 'left_self_dimer_dg', 'float',
        'Free energy (delta-G) of the left primer self-dimer (two copies of the left primer binding to each other) in kcal/mol. '
        'More negative = more stable self-dimer = worse primer. Good primers typically have dG > -5 kcal/mol.',
        '-1.2, -3.5, -8.0')
    
    add(S, 42, 'right_self_dimer_dg', 'float',
        'Free energy of the right primer self-dimer in kcal/mol.',
        '-0.8, -2.9, -6.1')
    
    add(S, 43, 'heterodimer_dg', 'float',
        'Free energy of the heterodimer (left primer + right primer binding to each other) in kcal/mol. '
        'More negative = more stable heterodimer = worse primer pair (primers may bind each other instead of template). '
        'Critical metric for primer pair quality.',
        '-1.5, -4.0, -9.2')
    
    add(S, 44, 'amplicon_mfe', 'float',
        'Minimum free energy of the PCR amplicon secondary structure in kcal/mol. '
        'Stable amplicon secondary structures can reduce PCR efficiency.',
        '-5.0, -12.3, -25.1')
    
    add(S, 45, 'primer_dimer_prediction', 'string',
        'Categorical prediction of primer dimer formation risk based on thermodynamic calculations. '
        'low = minimal dimer risk. medium = moderate risk, may need optimization. high = significant risk, '
        'likely to produce artifact products.',
        'low, medium, high')
    
    add(S, 46, 'primer_dimer_end_min_dg', 'float',
        'Minimum free energy at the 3-prime ends of primer dimer configurations in kcal/mol. '
        'More negative = stronger 3-prime dimer = worse for PCR, because 3-prime dimers are extendable by '
        'DNA polymerase and produce artifact products.',
        '-0.5, -2.3, -5.1')
    
    add(S, 47, 'heterodimer_full_dg', 'float',
        'Full-length free energy of heterodimer formation between left and right primers in kcal/mol. '
        'More comprehensive than heterodimer_dg, considering the full primer-primer alignment.',
        '-2.0, -5.5, -10.8')
    
    add(S, 48, 'left_homodimer_dg', 'float',
        'Free energy of the left primer homodimer (two identical left primers binding in parallel) in kcal/mol.',
        '-1.0, -3.2, -7.5')
    
    add(S, 49, 'right_homodimer_dg', 'float',
        'Free energy of the right primer homodimer in kcal/mol.',
        '-0.9, -2.8, -6.3')
    
    add(S, 50, 'common_variant_score', 'integer',
        'Count of common SNPs (minor allele frequency >= 1%) overlapping either primer binding site '
        '(forward or reverse). Computed from UCSC hg19 dbSNP snp151Common track. '
        'Primer binding sites are determined by searching the bisulfite-converted genome within the DMR block '
        '(plus 500bp flanking) for each primer sequence. '
        'Common SNPs in primer binding regions can cause allele-specific PCR failure: the primer may not bind '
        'to the variant allele, leading to false-negative methylation calls in carriers. '
        '0 = no common SNPs in primer binding regions (ideal). '
        'Only computed for primers that pass ePCR screening. NaN for unscreened or failing primers.',
        '0, 1, 3, 5')

    add(S, 51, 'snp_count', 'integer',
        'Total count of known SNPs (any MAF, including rare variants) overlapping either primer binding site. '
        'Computed from UCSC hg19 dbSNP snp151 track (all known SNPs, not just common ones). '
        'Includes rare variants not in snp151Common. '
        '0 = no known SNPs in primer binding regions (ideal). '
        'Only computed for primers that pass ePCR screening. NaN for unscreened or failing primers.',
        '0, 1, 2, 4')

    add(S, 52, 'segdup_overlap', 'boolean',
        'True if either primer binding site overlaps a segmental duplication region. '
        'Computed from UCSC hg19 genomicSuperDups track. '
        'Segmental duplications are large DNA blocks (>1kb) with high sequence similarity (>90%) '
        'found at multiple genomic locations. Primers binding within segmental duplications have '
        'increased risk of non-specific amplification at paralogous sites. '
        'False = no segmental duplication overlap (ideal). '
        'Only computed for primers that pass ePCR screening. NaN for unscreened or failing primers.',
        'True, False')

    add(S, 53, 'repeat_overlap', 'boolean',
        'True if either primer binding site overlaps a repetitive element. '
        'Computed from UCSC hg19 RepeatMasker (rmsk) track. '
        'Repetitive elements (SINE/Alu, LINE/L1, LTR/ERV, simple repeats, low complexity) are found at '
        'many genomic locations and can cause non-specific priming. '
        'False = no repeat overlap (ideal). '
        'Only computed for primers that pass ePCR screening. NaN for unscreened or failing primers.',
        'True, False')

    add(S, 54, 'repeat_family', 'string',
        'Name of the repetitive element family overlapping the primer binding region with the largest overlap. '
        'Computed from UCSC hg19 RepeatMasker (rmsk) track repFamily field. '
        'If both forward and reverse primers overlap repeats, the forward primer family is reported. '
        'Empty string if no repeat overlap. '
        'Common families: Alu (SINE), L1/L2 (LINE), MIR (SINE), ERVL/ERVK/ERVL-MaLR (LTR), '
        'simple_repeat, Low_complexity. '
        'Only computed for primers that pass ePCR screening. NaN for unscreened or failing primers.',
        'Alu, L1, MIR, ERVL, (empty)')

    # === ePCR columns (All_Primers and Selected_Panel only) ===
    S2 = 'All_Primers, Selected_Panel'

    add(S2, 0, 'ePCR_verdict', 'string',
        'PASS or FAIL verdict from the ePCR genome-wide specificity screen (phase4_bisulfite_ePCR.py). '
        'PASS = exactly one PCR product found genome-wide AND the product is on the genome state matching the primer template '
        '(SM/AM -> methylated, SU/AU -> unmethylated). '
        'FAIL = zero products (no amplification), multiple products (off-target amplification), '
        'one product on the wrong genome state (wrong_genome), or product on both genome states (not_methylation_specific). '
        'Blank for primers that were not screened (not in top 200 per cell type).',
        'PASS, FAIL, (blank)')

    add(S2, 1, 'ePCR_product_count', 'integer',
        'Number of PCR products found genome-wide by the ePCR algorithm. '
        '0 = no amplification. 1 = specific amplification (PASS if on correct genome state). '
        '>1 = off-target amplification (FAIL). '
        'Products are counted across both genome states (methylated CT_M and unmethylated CT) and both strands '
        '(Watson and Crick), with deduplication so a locus found in both states counts once. '
        'Blank for unscreened primers.',
        '0, 1, 2, 3, (blank)')

    add(S2, 2, 'ePCR_products', 'string',
        'Semicolon-separated list of PCR product coordinates found by ePCR. '
        'Format per product: chrom:start-end(size_bp,genome_state,strand,total_mismatches). '
        'genome_state: "methylated" (CT_M, CpG cytosines preserved) or "unmethylated" (CT, all cytosines converted to T) '
        'or "methylated+unmethylated" (product found on both states). '
        'strand: "W" (Watson/forward) or "C" (Crick/reverse). '
        'total_mismatches (mm=N): total mismatches across both primers at this product site. '
        'Blank for unscreened primers or primers with 0 products.',
        'chr17:4315368-4315523(155bp,methylated,W,mm=0)')

    add(S2, 3, 'ePCR_error', 'string',
        'Error message if the primer failed ePCR screening. Types: '
        '"wrong_genome: product on unmethylated, template expects methylated" (or vice versa) = product on wrong genome state. '
        '"not_methylation_specific: product on both methylated and unmethylated" = primer amplifies on both states, '
        'cannot distinguish methylation status. '
        '"Invalid base X in primer" = non-ACGT character in primer sequence. '
        '"Primer too short (min 16 bp)" = primer shorter than the 16-mer seed length. '
        'Blank for PASS and for product-count FAILs (0 or >1 products).',
        'wrong_genome: product on unmethylated, template expects methylated, (blank)')

    add(S2, 4, 'ePCR_elapsed_sec', 'float',
        'Wall-clock time in seconds for the ePCR screening of this primer pair. '
        'Includes genome index lookup, seed matching, full-primer verification, and product pairing '
        'across all chromosomes, both strands, and both methylation states. '
        'Blank for unscreened primers.',
        '0.05, 12.3, 45.6')

    add(S2, 5, 'selected_for_panel', 'integer (0/1)',
        'Binary flag indicating whether this primer was selected for the final panel. '
        '1 = selected. 0 = not selected. '
        'Selection criteria: ePCR_verdict = PASS (exactly 1 product on the intended genome state only). '
        'ALL passing primers are kept (no DMR deduplication, no cell-type cap). '
        'Primers are ranked by delta_means descending, then cpg_total descending within each cell type.',
        '0, 1')

    # === DMR_Blocks sheet columns ===
    S3 = 'DMR_Blocks'

    add(S3, 0, 'cell_type', 'string', 'Same as All_Primers.cell_type. The cell type this DMR discriminates.', 'Adipocytes, Neuron')
    add(S3, 1, 'source', 'string', 'Pipeline version.', 'v2.2.8')
    add(S3, 2, 'seq_id', 'string', 'DMR-level identifier. Format: {cell_type}_{rank:04d}_{M or U}.', 'Adipocytes_0001_U')
    add(S3, 3, 'rank', 'integer', 'DMR rank within cell type by delta_means descending. Rank 1 = best discrimination.', '1, 2, 3')
    add(S3, 4, 'chrom', 'string', 'Chromosome (hg19).', 'chr1, chr17')
    add(S3, 5, 'dmr_start', 'integer', 'DMR start position (1-based, hg19).', '4315368')
    add(S3, 6, 'dmr_end', 'integer', 'DMR end position (1-based, inclusive, hg19).', '4315523')
    add(S3, 7, 'gene', 'string', 'Nearest gene to the DMR.', 'SPNS3, LHPP')
    add(S3, 8, 'annotation', 'string', 'Genomic feature annotation relative to nearest gene.', 'Intergenic, intron')
    add(S3, 9, 'delta_means', 'float', 'Methylation difference (target vs background). Primary DMR ranking metric.', '0.779, 0.742')
    add(S3, 10, 'dms_score', 'float', 'DMS score = mean(A,B,C,D) * E. See All_Primers.dms_score for full formula.', '0.8147, 0.7952')
    add(S3, 11, 'target_mean_meth', 'float', 'Mean methylation in target cell type (0-1).', '0.102, 0.752')
    add(S3, 12, 'background_mean_meth', 'float', 'Mean methylation across background cell types (0-1).', '0.881, 0.141')
    add(S3, 13, 'num_cpgs', 'integer', 'Number of CpG sites in the DMR block.', '8, 5, 3')
    add(S3, 14, 'block_len', 'integer', 'DMR block length in base pairs.', '155, 46, 127')
    add(S3, 15, 'num_primers', 'integer',
        'Number of primer pairs designed for this specific DMR block. '
        'More primers = more candidates for ePCR screening and selection.', '13, 7, 1, 10')

    # === Summary sheet columns ===
    S4 = 'Summary'

    add(S4, 0, 'cell_type', 'string', 'Cell type name.', 'Adipocytes, Blood-T-CD4')
    add(S4, 1, 'source', 'string', 'Pipeline version.', 'v2.2.8')
    add(S4, 2, 'num_blocks', 'integer', 'Number of DMR blocks identified for this cell type.', '575, 1048')
    add(S4, 3, 'num_primers', 'integer', 'Total primer pairs designed for this cell type across all DMRs.', '2775, 18738')
    add(S4, 4, 'num_low_dimer', 'integer', 'Primers with low primer dimer formation risk.', '13494, 5000')
    add(S4, 5, 'num_medium_dimer', 'integer', 'Primers with medium primer dimer formation risk.', '251548, 30000')
    add(S4, 6, 'num_high_dimer', 'integer', 'Primers with high primer dimer formation risk.', '126528, 15000')

    return pd.DataFrame(DD, columns=['Sheet', 'No', 'Column', 'Type', 'Description', 'Example_Values'])


if __name__ == '__main__':
    df = build_data_dictionary()
    print('Data dictionary: ' + str(len(df)) + ' rows')
    print(df[['Sheet', 'Column', 'Type']].to_string())
