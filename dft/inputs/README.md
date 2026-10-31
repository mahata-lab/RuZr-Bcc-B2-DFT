# VASP inputs

| File | Used for |
|---|---|
| `INCAR.relax` + `KPOINTS` | full relaxation of each SQS from `POSCAR_ideal` |
| `INCAR.elastic` + `KPOINTS` | 12 static strained cells (±1% e1–e6) built from `CONTCAR_relaxed` |
| `INCAR.dos` + `KPOINTS.dos` | DOS of the relaxed structure |

The strained lattices are given in `structures/<id>/strained_cells_stress.csv` (columns a00–a22, row vectors in Å);
atomic positions in direct coordinates are those of `CONTCAR_relaxed`. See `POTCAR_info.md` for the PAW-PBE datasets.
