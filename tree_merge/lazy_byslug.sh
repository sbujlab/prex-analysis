#!/bin/sh
root -b -q 'MergeSumTree.C('$1')';
root -b -q 'MergePostpan.C('$1')';

