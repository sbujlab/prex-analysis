#! /bin/bash
root -b -q 'EigenDitSlope6bpm.C('$1')';
root -b -q 'DrawEigVecFromDance6bpm.C('$1')';
root -b -q 'CompareSlopes6bpm.C('$1')';
root -b -q 'EigenDitSlope.C('$1')';
root -b -q 'DrawEigVecFromDance.C('$1')';
root -b -q 'CompareSlopes.C('$1')';
