#!/bin/sh
for i in {1..94};
do 
    root -b -q 'MergeSumTree.C('$i')';
    root -b -q 'MergePostpan.C('$i')';
done
