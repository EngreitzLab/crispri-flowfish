## Estimate effect sizes of each guide by calling the MLE on bin data

import os
import pandas as pd
import numpy as np

def get_sortparams_file(wildcards):
	currSamples = samplesheet.loc[(samplesheet['ExperimentIDReplicates'] == wildcards.ExperimentIDReplicates) & (samplesheet['Bin'].isin(binList))]
	Batch = currSamples['Batch'].unique()
	# SampleNumber = currSamples['SampleNumber'].unique()
	SampleNumber = currSamples['Sample'].unique()
	if (len(Batch) != 1) or (len(SampleNumber) != 1):
		print(currSamples['SampleID'])
		raise ValueError("Found more than one possible sort params file path. Correct the samplesheet and rerun.")
	return os.path.join(config['sortparamsdir'], Batch[0] + "_" + SampleNumber[0] + ".txt")


# run mle to calculate the effect sizes
rule calculate_allelic_effect_sizes:
	input:
		counts='results/byExperimentRep/{ExperimentIDReplicates}.bin_counts.txt',
		design=config["design_file"],
		sortparams=get_sortparams_file 
	output:
		'results/byExperimentRep/{ExperimentIDReplicates}.raw_effects.txt'
	log:
		'results/byExperimentRep/{ExperimentIDReplicates}.mle_log.txt'
	shell:
		"""
		Rscript crispri-flowfish/workflow/scripts/estimate_effect_sizes.R \
			 --designDocLocation {input.design} \
			 --countsLocation {input.counts} \
			 --sortParamsloc {input.sortparams} \
			 --outputmle {output} --log {log}
		"""
