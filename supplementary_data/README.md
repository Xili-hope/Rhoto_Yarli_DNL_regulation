# Supplementary data

These files accompany Li et al., *Lipid accumulation in nitrogen and phosphorus-limited yeast is caused by less growth-related dilution*, Metabolic Engineering 93 (2026), 60–72. DOI: [10.1016/j.ymben.2025.08.010](https://doi.org/10.1016/j.ymben.2025.08.010).

| File | Contents |
| --- | --- |
| [Full supplementary data](Li_et_al_Supplementary_Data_2026_push.xlsx) | The author's supplied 34-sheet upload copy, preserved byte for byte. Includes figure data, protein abundance, reporter signals and baseline calibration. |
| [Protein abundance](Li_et_al_Protein_Abundance_2026.xlsx) | Six plain worksheets: a contents/units guide, absolute abundance and proteome fractions for both organisms, and sample information. The five data worksheets preserve the supplied workbook's values. |

Release date: 29 October 2026. File integrity hashes are in [CHECKSUMS.sha256](CHECKSUMS.sha256).

## Protein abundance

Rt denotes *Rhodotorula toruloides* NBRC 0880; Yl denotes *Yarrowia lipolytica* W29. B, N and P denote nutrient-replete, nitrogen-limited and phosphorus-limited conditions. Each condition has three TMT samples.

Absolute abundance is reported as **g individual protein / g dry cell weight**. Proteome fractions are **g individual protein / g quantified proteome**. Sample-specific S7/S8 biomass weights are provided in `Table 1. sample_information`. Means use three replicates and SD uses ddof = 1. Replicate pairing follows the original plotting code and is not independently documented.

The data contain 2,487 quantified Rt proteins and 2,319 quantified Yl proteins. All 254 unmatched Rt and 70 unmatched Yl identifications remain present with blank abundance; blanks are not zero. Each quantified proteome-fraction column sums to one. Raw TMT reporter signals and UPS2-calibrated baseline concentrations are available in Tables 6 and 7 of the full workbook. [PRIDE accession PXD063019](https://www.ebi.ac.uk/pride/archive/projects/PXD063019) provides the public proteomics archive.

## Data notes

The supplied upload copy omits the working sheets `source_notes`, `source_files`, `reconstruction_checks` and `lipid_fit_check`. Some workbook text still refers to those sheets or their issue identifiers. The source index and relevant notes are recorded here. The supplied workbook itself has not been rewritten.

- **P01:** Rt protein A0A0K3C6X4 was explicitly excluded by the original analysis notebook. Its raw reporter signals remain in Table 6a.
- **P02:** Yl Figure 4 uses condition-average biomass weights, while Tables 4a/4b use the sample-specific S7/S8 weights. Both weight sets are listed in Table 1; they should not be silently interchanged.
- **P03:** The 324 unmatched TMT identifications lack a matched absolute baseline. Their abundance values remain blank.
- **F01, recovered:** Yl P1 lipid was recovered from `Si_data_revision2025.xlsx`, `biomass_pool!R6` (110.0825925984112 µg/OD600 mL), and confirmed against the sum of 16 named fatty acids in `yl_P_pool!B2:B17`. All three Yl P biomass replicates are complete.
- **F02, recovered:** `growth_rate_revision!A14:E16` explicitly maps lipid replicates Yl N4/N5/N6 to growth variables Yl_N1/Yl_N2/Yl_N3 at 0.10/0.11/0.11 h⁻¹. X-condition growth rates are also recovered from this sheet. Table 13 includes 18 complete B/N/P observations; historical Yl B4 remains a row with missing pool/growth values and is excluded from the model calculations.
- **Model comparison:** Table 13 gives fixed-growth R² = 0.16811, agreeing with the paper's reported 0.17 at its displayed precision. Fixed-flux R² = 0.49841 differs from the reported 0.51 and remains to be reconciled. These recovered model values should not be described as an exact reproduction of every published fit.
- **F03:** Yl control growth measurements come from the glucose-control source sheet referenced by the plotting notebook. The notebook's final execution order is not recorded.
- **F04/F05:** Table 16c reproduces an earlier protein-flux calculation with condition-average pool weights and four control replicates. Table 16d preserves the later three-replicate literal model arrays; its Yl B1/P1 pool values differ from the S8 weights. These historical versions are kept distinct.
- **F06:** Table 15a retains an unlabeled numerical total row in the original Yl P composition input. The original normalization includes this row, doubling that denominator. This historical normalization has not been silently corrected.
- **F07:** Table 14a distinguishes observed MID from archived fitted MID. Some archived fit labels/parameters differ from the later plotting inputs; archived fitted rows are not certified as the final published curves.

The supplied copy removes the source-sheet column from Tables 16b/16c. The source blocks below identify those retained rows; row numbers include the header row.

| Source sheet in `Si_data_revision2025.xlsx` | Table 16b rows | Table 16c rows |
| --- | --- | --- |
| `rt_b_aa` | 2–13 | 2–5 |
| `rt_n_aa_24may` | 14–25 | 6–9 |
| `rt_p_aa_dec23` | 26–34 | 10–12 |
| `yl_b_aa` | 35–46 | 13–16 |
| `yl_n_aa` | 47–55 | 17–19 |
| `yl_p_aa` | 56–64 | 20–22 |

## Verification

All retained measurement values were compared with the verified reconstruction; differences consist of removed review fields/sheets, text whitespace cleanup and removal of the unquantified historical Yl B4 row from Table 15c. The protein extract was checked against the supplied workbook, including identifier uniqueness, all nine replicate columns, proteome-fraction sums, biomass scaling, means and sample SD. The Yl P1 and Yl nitrogen recoveries remain present. Remaining historical/version issues are listed above.

## Source identifiers

Source IDs in the full workbook's table of contents identify the original analysis files listed below. Paths are relative to the archived project and do not imply that every source file is included in this repository.

| Source ID | Original analysis file |
| --- | --- |
| S001 | `DNLpp_figs_all/Yl_proteomics_alignment_paper/2025_new_data/output/2025_rhoto_proteomics_fraction_mapped.xlsx` |
| S002 | `DNLpp_figs_all/Yl_proteomics_alignment_paper/2025_new_data/protein_quant/RhoToruprotein_quant_2048.tsv` |
| S003 | `revision1/scripts/heatmap.ipynb` |
| S004 | `DNLpp_figs_all/Yl_proteomics_alignment_paper/rt/112018_RT_spikein_lf_proteinquant.xlsx` |
| S005 | `DNLpp_figs_all/proteomics/2025_final_scrips_revisit/2025revision_rt_lipid.csv` |
| S006 | `DNLpp_figs_all/Yl_proteomics_alignment_paper/2025_new_data/output/20250408_Yl_proteomics_fraction_mapped.xlsx` |
| S007 | `DNLpp_figs_all/Yl_proteomics_alignment_paper/2025_new_data/protein_quant/CanLipo_protein_quant_2049.tsv` |
| S008 | `DNLpp_figs_all/Yl_proteomics_alignment_paper/112018_YL_spikein_lf_proteinquant.xlsx` |
| S009 | `DNLpp_figs_all/proteomics/2025_final_scrips_revisit/2025yl_revision.csv` |
| S010 | `DNLpp_figs_all/proteomics/2025_final_scrips_revisit/2025_final_rt_prot_allocation.ipynb` |
| S011 | `DNLpp_figs_all/proteomics/2025_final_scrips_revisit/2025_final_yl_prot_allocation.ipynb` |
| S012 | `revision_data/Si_data_revision2025.xlsx` |
| S013 | `revision_data/Revision_figures/Fig2/dnl_flux_df_with_growth.csv` |
| S014 | `revision1/AcCoA/20231225_yl_P_MID.xlsx` |
| S015 | `revision1/AcCoA/202507_batch8_2_rtBN_dnl_pblksub_cor_MID.xlsx` |
| S016 | `revision1/AcCoA/202507_batch8_RtP_DNL_pblksub_cor_MID.xlsx` |
| S017 | `revision1/AcCoA/2025yL_b_REVISIT_pblksub_cor_MID.xlsx` |
| S018 | `revision1/AcCoA/2025ylX_revision_MID.xlsx` |
| S019 | `revision1/AcCoA/RTB_even_chain_fa_pblksub_cor_MID.xlsx` |
| S020 | `revision1/AcCoA/batch6_lipid_kinetics_pblksub_cor_MID.xlsx` |
| S021 | `revision_data/2025july_revision_evenodd_plot_20250707.ipynb` |
| S022 | `DNLpp_figs_all/dnl_all_figs.ipynb` |
| S023 | `DNLpp_figs_all/SI_data_rt_yl_result_summary_XI_to_Daniel.xlsx` |
| S024 | `DNLpp_figs_all/growth_rate_estimation_plot.xlsx` |
| S025 | `DNLpp_figs_all/chemo_from_publication.xlsx` |
| S026 | `revision1/mu_carbon_accoa_MID.csv` |
