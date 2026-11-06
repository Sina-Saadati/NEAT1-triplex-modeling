# NEAT1 ENE+A-rich triplex: structure prediction, MD and pocket analysis

Computational pipeline and **preliminary** results for the 3′-end triple-helix element of the lncRNA NEAT1 (101 nt, GQ859162.1: 22643–22743), for which no experimental or computational 3D structure exists.

Presented as a poster: *Structural modeling and dynamics of the NEAT1 lncRNA triplex reveal a candidate druggable region*, 19th National and 4th International Conference of Iran Biophysical Chemistry.

## Summary of results
- **Conservation (19 mammals):** U-rich loop (nt 14–23) 98.9%, A-rich tract (nt 82–101) 97.6%; U15 and U19 invariant; nt 85–101 identical in all sequences. Conservation is high across the whole domain (93.2%).
- **Thermodynamics (RNAfold):** MFE −35.00 kcal/mol (free) vs −24.80 kcal/mol (both regions single-stranded); pre-organization estimate ΔΔG = +10.2 kcal/mol.
- **Four template-free 3D predictors** (RhoFold+ with and without MSA, trRosettaRNA, RNAComposer) diverge strongly: pairwise C1′ RMSD 11.13–35.96 Å. Mean pLDDT: 61.2 (RhoFold+ no MSA), 68.7 (RhoFold+ with MSA), 41.1 (trRosettaRNA).
- **MD (GROMACS, AMBER χOL3, TIP3P):** 100 ps production plus a 500 ps extension, started from the RhoFold+ (no MSA) model, the model with the most intact backbone. No RMSD plateau; the two regions remain linked by an edge-localized hydrogen-bond network (mean 7.06 bonds, range 3–11, never zero in the extended run).
- **Pockets (fpocketR, fpocket 4.0.3):** best druggability 0.141 (static RhoFold+ model) vs 0.521 (MD frame at 100 ps; A-rich tract nt 86, 88–90, plus nt 7 and 10). The static score is unchanged (0.141) when the B-factor column is zeroed, so the increase is not an artifact of B-factor handling. The RNAComposer model has a pocket scoring 0.555 outside both functional regions; trRosettaRNA has none.
- **Docking (AutoDock Vina) and MM-GBSA:** Vina scores −8.86 (kanamycin A), −8.59 (quercetin) and −6.77 (caffeine) kcal/mol; single-snapshot MM-GBSA for quercetin gives ΔG = −1.72 kcal/mol (no entropy term).

## Limitations
- MD is short (100 + 500 ps) and single-replica, so structural conclusions are exploratory.
- fpocket's druggability score was calibrated on proteins and is only a heuristic for RNA; the 0.5 cutoff is a conventional guide, not a validated threshold for RNA.
- Vina scores do not demonstrate specificity (larger ligands score better; see ligand efficiency in the figure). The quercetin score is compared only with a reported ΔG for a different RNA (MALAT1) and is not a validation.
- MM-GBSA is a single snapshot and supports binding only weakly.
- The 3D predictors disagree strongly, and none reproduces the triple-helix tertiary contacts.

## Repository structure
- `notebooks/` – Google Colab notebook for sequence retrieval, folding, BLAST/MUSCLE, RhoFold+, MD, fpocketR and docking
- `analysis/` – R scripts for the figures
- `data/` – sequences, alignment, model PDB files, small MD outputs (.xvg), fpocketR tables and docking summaries
- `figures/` – final figures
- `environment.yml` – software versions

## Software
GROMACS 2021.3 · AMBER parmbsc0+χOL3 for GROMACS ([srnas/ff](https://github.com/srnas/ff)) · RhoFold+ · trRosettaRNA · RNAComposer · ViennaRNA · BLAST+ / MUSCLE · fpocketR 1.3.6 (fpocket 4.0.3) · AutoDock Vina 1.2.7 · R (ggplot2, patchwork)

## How to reproduce
1. Run the notebook in `notebooks/` in order.
2. Run the scripts in `analysis/`, setting `DATA_DIR` at the top of each script to your local data folder.

## Citation
Saadati, S., Rashtbari, S., et al. Structural modeling and dynamics of the NEAT1 lncRNA triplex reveal a candidate druggable region. Conference poster, 2026.

## License
Code: MIT. Data and figures: CC BY 4.0.

## Contact
Sina Saadati – sina.saadati.ac@gmail.com – GitHub: [Sina-Saadati](https://github.com/Sina-Saadati)
