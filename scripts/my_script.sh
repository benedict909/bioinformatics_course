#!/usr/bin/env bash

mkdir -p data analysis docs

for my_dir in scripts analysis docs data;do
	touch ${my_dir}/README.md 
	echo "# ${my_dir}" >> ${my_dir}/README.md 
done


