#!/bin/bash
#grompp EM0
gmx grompp -f martini_em0.mdp -c $1.gro -p $1"2.top" -o em0_$1.tpr -maxwarn 1 > em0_grompp_$1.log 2>&1
#EM0
gmx mdrun -ntomp 1 -deffnm em0_$1 -nt 1

#grompp EM1
gmx grompp -f martini_em1.mdp -c em0_$1.gro -p $1"2.top" -o em1_$1.tpr -maxwarn 1 > em1_grompp_$1.log 2>&1
#EM1
gmx mdrun -ntomp 1 -deffnm em1_$1 -nt 1

#grompp EM2
gmx grompp -f martini_em2.mdp -c em1_$1.gro -p $1"2.top" -o em2_$1.tpr -maxwarn 1 > em2_grompp_$1.log 2>&1
#EM2
gmx mdrun -ntomp 1 -deffnm em2_$1 -nt 1

#Edit topology to include charges
python UpdateTopolCharge.py -p $1"2.top"

#grompp EM3
gmx grompp -f martini_em3.mdp -c em2_$1.gro -p $1"2.top" -o em3_$1.tpr > em3_grompp_$1.log 2>&1
#EM3
gmx mdrun -ntomp 1 -deffnm em3_$1 -nt 1

#grompp NVT eq 
gmx grompp -f martini_nvt.mdp -c em3_$1.gro -p $1"2.top" -o nvteq_$1.tpr > nvteq_grompp_$1.log 2>&1

#NVT Eq
gmx mdrun -ntomp 1 -deffnm nvteq_$1 

#grompp NPT md
gmx grompp -f martini_npt.mdp -c nvteq_$1.gro -p $1"2.top" -o md_$1.tpr  > md_grompp_$1.log 2>&1
cp md.mdp md_$1.mdp
#gmx grompp -f ST_$1.mdp -c nvteq_$1.gro -p $1.top -o ST_$1.tpr  > ST_grompp_$1.log 2>&1
#cp ST.mdp ST_$1.mdp

