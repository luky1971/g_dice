#!/bin/bash

{ time g_ensemble_comp -f1 ../test/PDZ2_frag_apo.pdb -f2 ../test/PDZ2_frag_bound.pdb ; } 2>&1 | tail -3
