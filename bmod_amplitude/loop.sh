#! /bin/sh
while IFS= read -r line; do
    root -b -q 'CheckBmodAmplitude.C('$line')';
done < all.list
