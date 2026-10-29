# Rhoto_Yarli_DNL_regulation
Python and MATLAB scripts and supplementary data for Li et al., **Lipid accumulation in nitrogen and phosphorus-limited yeast is caused by less growth-related dilution**, *Metabolic Engineering* 93 (2026), 60–72. [Published paper](https://doi.org/10.1016/j.ymben.2025.08.010).

## Supplementary data

- [Full supplementary data workbook](supplementary_data/Li_et_al_Supplementary_Data_2026_push.xlsx): figure data, protein abundance, raw TMT signals and calibration tables.
- [Protein abundance workbook](supplementary_data/Li_et_al_Protein_Abundance_2026.xlsx): absolute abundance and proteome fractions for *Rhodotorula toruloides* and *Yarrowia lipolytica*, with sample definitions and units.
- [Data notes](supplementary_data/README.md): units, processing, and remaining figure/version differences.

The protein workbook includes 2,487 quantified Rt proteins and 2,319 quantified Yl proteins, each across nine samples. Unmatched identifications retain blank abundance values. Raw proteomics are deposited in [PRIDE, PXD063019](https://www.ebi.ac.uk/pride/archive/projects/PXD063019).

## Analysis code

- `MID_fitting/`: mass-isotopologue distribution fitting.
- `protein_function_assignment_Rhoto/`: protein function assignment for *R. toruloides*.
- `protein_function_assignment_Yarli/`: protein function assignment for *Y. lipolytica*.
