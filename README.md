# RuZr-BCC-DFT

Data and analysis code for

> A. Mahata, *Lattice Instability and Alloying Design of the Zr-Rich BCC Matrix in Zr–Ru BCC+B2 Eutectic Alloys from First Principles* (submitted to Acta Materialia).

The repository contains the DFT structures and raw stress/energy/DOS data for 14 special quasirandom structures (54-atom BCC cells) of the Zr–7.4 at.% Ru matrix and 12 alloyed variants, together with the Python scripts that produce every number, table and figure in the paper.

## Contents

```
dft/
  inputs/           INCAR files (relaxation, elastic, DOS) and KPOINTS files used for all calculations
  inputs/POTCAR_info.md   PAW-PBE datasets used (POTCAR files are not distributed; VASP license)
  make_potcar.sh    builds a POTCAR in POSCAR species order from a local VASP PAW library
  structures/<id>/
    POSCAR_ideal            starting SQS on ideal BCC sites (a0 = 3.57 Å)
    CONTCAR_relaxed         fully relaxed structure (positions, shape, volume)
    strained_cells_stress.csv   energy, stress (kbar, VASP order xx yy zz xy yz zx) and lattice of the 12 strained cells
data/
  dataset.json      all calculations in one file: species, relaxed cell, positions, relaxation summary,
                    relaxation trajectory (energy and positions at every ionic step), strained-cell stresses
  dos/<id>.npz      total and site-projected (s, p, t2g, eg) DOS, summed per element, E relative to E_F
code/
  ruzr_analysis.py  elastic tensors, moduli, static displacement field, omega-mode projection,
                    random baseline, local Ru environment, DOS descriptors, classification and fits
  make_figures.py   Figures 1–6 and Supplementary Figures S1–S2
  make_tables.py    LaTeX tables of the paper and the Supplementary Material
results/            output of the scripts (CSV, fits.txt, latex/)
figures/            output of make_figures.py
```

Configuration ids are `<composition>_c<n>`; `c1` and `c2` are two independent SQS configurations of the same composition (computed for Zr50Ru4). Only relaxations that met the force criterion of 0.02 eV/Å are included.

| Family | Composition | X |
|---|---|---|
| base | Zr50Ru4 | – |
| Z (X for Zr) | Zr48Ru4X2 | Al, Hf, Re |
| R (X for Ru) | Zr50Ru2X2 | Al, Sc, Ti, V, Hf, Y |
| R* (Ru-lean, X-rich) | Zr48Ru2X4 | Ti, Ta, Y |

## Reproducing the results

```bash
pip install -r requirements.txt
python code/ruzr_analysis.py     # -> results/*.csv, results/fits.txt
python code/make_figures.py      # -> figures/
python code/make_tables.py       # -> results/latex/
```

Run the scripts from the repository root. The full pipeline takes under a minute on a laptop.

## Calculation settings (summary)

VASP 6, PAW-PBE, non-spin-polarized, 54-atom 3×3×3 BCC SQS.

| | Relaxation | Elastic (strained cells) | DOS |
|---|---|---|---|
| ENCUT | 320 eV | 320 eV | 450 eV |
| EDIFF | 1e-4 eV | 1e-6 eV | 1e-5 eV |
| Ionic | ISIF 3, EDIFFG −0.02 eV/Å, NSW 150 | static, ISYM 0 | static |
| Smearing | MP-1, 0.2 eV | MP-1, 0.1 eV | tetrahedron (Blöchl) |
| k-points | Γ 4×4×4 | Γ 4×4×4 | MP 8×8×8 |

Elastic constants are obtained from ±1% normal and engineering shear strains by a least-squares fit of the stress–strain relation including a residual stress, with σ = −σ_VASP (tension positive). They are clamped-ion values.

## License

Code: MIT (see `LICENSE`). Data: CC BY 4.0.

## Citation

If you use these data or scripts, please cite the paper above.
