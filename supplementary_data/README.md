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

- **P01:** Rt protein A0A0K3C6X4 was explicitly excluded by the original analysis notebook. Its raw reporter signals remain in Table 6a.
- **P02:** Yl Figure 4 uses condition-average biomass weights, while Tables 4a/4b use the sample-specific S7/S8 weights. Both weight sets are listed in Table 1; they should not be silently interchanged.
- **P03:** The 324 unmatched TMT identifications lack a matched absolute baseline. Their abundance values remain blank.
- **F01, recovered:** Yl P1 lipid is 110.0825925984112 µg/OD600 mL, confirmed against the sum of 16 named fatty acids. All three Yl P biomass replicates are complete.
- **F02, recovered:** Lipid replicates Yl N4/N5/N6 correspond to growth variables Yl_N1/Yl_N2/Yl_N3 at 0.10/0.11/0.11 h⁻¹. X-condition growth rates are also included. Table 13 includes 18 complete B/N/P observations; historical Yl B4 remains a row with missing pool/growth values and is excluded from the model calculations.
- **Model comparison:** Table 13 gives fixed-growth R² = 0.16811, agreeing with the paper's reported 0.17 at its displayed precision. Fixed-flux R² = 0.49841 differs from the reported 0.51 and remains to be reconciled. These recovered model values should not be described as an exact reproduction of every published fit.
- **F03:** Yl control growth measurements come from the glucose-control source sheet referenced by the plotting notebook. The notebook's final execution order is not recorded.
- **F04/F05:** Table 16c reproduces an earlier protein-flux calculation with condition-average pool weights and four control replicates. Table 16d preserves the later three-replicate literal model arrays; its Yl B1/P1 pool values differ from the S8 weights. These historical versions are kept distinct.
- **F06:** Table 15a retains an unlabeled numerical total row in the original Yl P composition input. The original normalization includes this row, doubling that denominator. This historical normalization has not been silently corrected.
- **F07:** Table 14a distinguishes observed MID from archived fitted MID. Some archived fit labels/parameters differ from the later plotting inputs; archived fitted rows are not certified as the final published curves.

## Verification

All retained measurement values were compared with the verified reconstruction; differences consist of removed review fields/sheets, text whitespace cleanup and removal of the unquantified historical Yl B4 row from Table 15c. The protein extract was checked against the supplied workbook, including identifier uniqueness, all nine replicate columns, proteome-fraction sums, biomass scaling, means and sample SD. The Yl P1 and Yl nitrogen recoveries remain present. Remaining historical/version issues are listed above.
