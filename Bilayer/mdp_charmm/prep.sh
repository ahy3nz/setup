#!/bin/bash

gmx grompp -f em.mdp -c $1.gro -p $1.top -o em -maxwarn 2 &> em_grompp.log
gmx mdrun -deffnm em

gmx grompp -f nvt.mdp -c em.gro -p $1.top -o nvt -maxwarn 2 &> nvt_grompp.log
gmx mdrun -deffnm nvt

gmx grompp -f npt_500ps.mdp -c nvt.gro -p $1.top -o npt -t nvt.cpt -maxwarn 2 &> npt_500ps_grompp.log
gmx mdrun -deffnm npt_500ps

gmx grompp -f npt.mdp -c npt_500ps.gro -p $1.top -o npt -t npt_500ps.cpt -maxwarn 2 &> npt_grompp.log
