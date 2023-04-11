#!/bin/bash

cd analysis-scripts
for script in 0*.R;
do
	echo "Running Rscript:"$script
	Rscript $script
done 

cd ..
