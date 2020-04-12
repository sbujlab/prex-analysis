#! /bin/bash
for i in {101..167}; do
    root -b -q 'ped_noise.cc('$i')';
done
