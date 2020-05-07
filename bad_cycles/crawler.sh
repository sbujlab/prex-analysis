#! /bin/sh
for i in {1..94}; do
    while IFS= read -r line; do
	if [[ $line -gt 0 ]]; then
            runnum=$line
            root -b -q -l './CycleCrawler.C('$runnum')';
	fi
    done < ../prex-runlist/simple_list/slug$i.list;
done










