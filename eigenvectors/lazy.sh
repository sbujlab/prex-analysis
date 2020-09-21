#! /bin/bash
# root -b -q 'EigenDitSlope.C(0,'$1')';
# root -b -q 'DrawEigVecFromDance.C(0,'$1')';
root -b -q 'CompareSlopes.C(0,'$1')';

# root -b -q 'EigenDitSlope.C(1,'$1')';
# root -b -q 'DrawEigVecFromDance.C(1,'$1')';
root -b -q 'CompareSlopes.C(1,'$1')';

# root -b -q 'EigenDitSlope.C(2,'$1')';
# root -b -q 'DrawEigVecFromDance.C(2,'$1')';
root -b -q 'CompareSlopes.C(2,'$1')';
