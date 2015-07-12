grompp -p water.top -f 0_1ps.mdp -o junk.tpr -c 0_1ns_boxedited.gro 
trjconv -f junk.xtc -o junk.gro -s junk.tpr

grompp and trjconv are tools of gromacs stored in gromacs-4.5.3/src/tools/

tools that read xtc and trr files:
trjconv, g_density, g_angle, g_dih

gmxdump will dump contents of binary file
ex. gmxdump -f <file>.xtc | less
