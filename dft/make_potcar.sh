#!/usr/bin/env bash
# Build a POTCAR in the species order of a POSCAR.
# Usage: PAW_PBE=/path/to/potpaw_PBE  dft/make_potcar.sh  POSCAR  [POTCAR]
set -euo pipefail
poscar=${1:?POSCAR}; out=${2:-POTCAR}; lib=${PAW_PBE:?set PAW_PBE to the potpaw_PBE directory}
declare -A DS=([Zr]=Zr [Ru]=Ru_pv [Al]=Al [B]=B [Sc]=Sc [Ti]=Ti_pv [V]=V_pv [Nb]=Nb_pv [Mo]=Mo [Hf]=Hf_pv [Ta]=Ta_pv [Y]=Y_sv [Re]=Re_pv)
: > "$out"
for el in $(sed -n 6p "$poscar"); do
  cat "$lib/${DS[$el]}/POTCAR" >> "$out"
done
grep TITEL "$out"
