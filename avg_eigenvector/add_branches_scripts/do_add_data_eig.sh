#!/bin/bash

analysis="seg"
dtype="_err" # OR nothing
bpms="" # OR "_allbpms"

if [ $# -ge 1 ] ; then
  analysis=$1;
fi
if [ $# -ge 2 ] ; then
  dtype=$2;
fi
if [ $# -ge 3 ] ; then
  bpms=$3;
fi

hadd -f /lustre19/expphy/volatile/halla/parity/crex-respin2/LagrangeOutput/rootfiles${dtype}_eigen_${analysis}_avg${bpms}/respin2-${analysis}_avg_sens-eigen_reg_lagr-basic.root /lustre19/expphy/volatile/halla/parity/crex-respin2/LagrangeOutput/rootfiles${dtype}_eigen_${analysis}_avg${bpms}/crexRespin2_eigen_*.000.root
root -l -b -q Add_RCDB_from_Agg_All.C'('\"/lustre19/expphy/volatile/halla/parity/crex-respin2/aggRootfiles/ErrorFlag/slugRootfiles/CREX-All-miniruns_units_maindet.root\"', '\"/lustre19/expphy/volatile/halla/parity/crex-respin2/LagrangeOutput/rootfiles${dtype}_eigen_${analysis}_avg${bpms}/respin2-${analysis}_avg_sens-eigen_reg_lagr-basic.root\"', '\"/lustre19/expphy/volatile/halla/parity/crex-respin2/LagrangeOutput/rootfiles${dtype}_eigen_${analysis}_avg${bpms}/respin2-rcdb-${analysis}_avg_sens-eigen_reg_lagr.root\"', '\"mini\"')'
root -l -b -q Add_crex_segment_eig.C'('\"/lustre19/expphy/volatile/halla/parity/crex-respin2/LagrangeOutput/rootfiles${dtype}_eigen_${analysis}_avg${bpms}/respin2-rcdb-${analysis}_avg_sens-eigen_reg_lagr.root\"', '\"/lustre19/expphy/volatile/halla/parity/crex-respin2/LagrangeOutput/rootfiles${dtype}_eigen_${analysis}_avg${bpms}/respin2-rcdb-parts-${analysis}_avg_sens-eigen_reg_lagr.root\"', '\"mini\"')'
root -l -b -q Add_slowControls_eig.C'('\"/lustre19/expphy/volatile/halla/parity/crex-respin2/LagrangeOutput/rootfiles${dtype}_eigen_${analysis}_avg${bpms}/respin2-rcdb-parts-${analysis}_avg_sens-eigen_reg_lagr.root\"', '\"/lustre19/expphy/volatile/halla/parity/crex-respin2/LagrangeOutput/rootfiles${dtype}_eigen_${analysis}_avg${bpms}/respin2-rcdb-parts-slowc-${analysis}_avg_sens-eigen_reg_lagr.root\"', '\"mini\"')'
root -l -b -q Add_dit_segment_eig.C'('\"/lustre19/expphy/volatile/halla/parity/crex-respin2/bmodOutput/slopes_run_avg_1X/respin2_All_dithering_slopes_run_avg_rcdb_segmented_pruned_usavg.root\"', '\"/lustre19/expphy/volatile/halla/parity/crex-respin2/LagrangeOutput/rootfiles${dtype}_eigen_${analysis}_avg${bpms}/respin2-rcdb-parts-slowc-${analysis}_avg_sens-eigen_reg_lagr.root\"', '\"/lustre19/expphy/volatile/halla/parity/crex-respin2/LagrangeOutput/rootfiles${dtype}_eigen_${analysis}_avg${bpms}/respin2-rcdb-parts-slowc-segments-${analysis}_avg_sens-eigen_reg_lagr.root\"', '\"mini\"')'
