#!/bin/bash

treemix -seed 123 -i maerl_SNPs_treemix.gz -o results -root Fal

for m in {1..7}
	do
	treemix -seed 123 -i maerl_SNPs_treemix.gz -o resultsM${m} -m ${m} -root Fal
done
