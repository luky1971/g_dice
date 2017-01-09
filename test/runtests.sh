#!/bin/bash

# TODO: residue eta tests

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# Clear previous test outputs
rm -f eta_atom_all.dat
rm -f eta_atom_grp1.dat

# All atoms
g_ensemble_comp -f1 ../tutor/PDZ2_frag_apo.pdb -f2 ../tutor/PDZ2_frag_bound.pdb -eta_atom eta_atom_all.dat

# Index group 1
g_ensemble_comp -f1 ../tutor/PDZ2_frag_apo.pdb -n1 ../tutor/PDZ2_frag_apo.ndx -f2 ../tutor/PDZ2_frag_bound.pdb -n2 ../tutor/PDZ2_frag_bound.ndx -eta_atom eta_atom_grp1.dat <<< "11"

if [[ $(diff eta_atom_all.dat exp_eta_atom_all.dat) ]]; then
    echo -e "${RED}Test case all atoms failed"
else
    echo -e "${GREEN}Test case all atoms passed"
fi

if [[ $(diff eta_atom_grp1.dat exp_eta_atom_grp1.dat) ]]; then
    echo -e "${RED}Test case index group 1 failed"
else
    echo -e "${GREEN}Test case index group 1 passed"
fi

# Turn off color
printf "${NC}"
