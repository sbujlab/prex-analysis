#!/bin/bash
filelist=""
shopt -s extglob
while IFS= read -r line; do
    all_files=$(ls -1 ./output/prexRespin2_bmw_extra_$line.root);
    for file in $all_files
    do
	echo $file;
	filelist+=" "$file;
    done

done < ./prex-runlist/slug$1.list
echo $filelist;

hadd -f ./treeMergeOutput/plus_bmw/slug$1.root $filelist;

shopt -u extglob
