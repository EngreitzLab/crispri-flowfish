from snakemake.utils import validate
import pandas as pd
import os
import glob

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# singularity: "docker://continuumio/miniconda3"



###########################################################################################
##### load config and sample sheets #####

# configfile: "config/config.yaml"   ## Read from command line instead
# validate(config, schema="../schemas/config.schema.yaml")


def find_fastq_files(samplesheet, fastqdir):
	## Adds columns 'fastqR1' to the sample sheet, only if they do not already exist

	for read in ["1"]:
		colName = 'fastqR' + read
		if not colName in samplesheet.columns:
			samplesheet[colName] = ""
			for i in samplesheet.index:
				currSample = samplesheet.at[i,'SampleID']
				file = glob.glob("{}_*_R{}_*.fastq.gz".format(os.path.join(fastqdir, currSample), read))
				if len(file) < 1:
					file = glob.glob("{}*.fastq.gz".format(os.path.join(fastqdir, currSample), read))
				if len(file) > 1:
					raise ValueError("Found more than one FASTQ file for sample :" + currSample)
				if len(file) == 1:
					samplesheet.at[i,colName] = file[0]
	
	return samplesheet


def add_experiment_names(samplesheet):
	if ('ExperimentIDReplicates' in samplesheet.columns) or ('ExperimentID' in samplesheet.columns) or ('ExperimentIDPCRRep' in samplesheet.columns):
		print("Warning: ExperimentID columns found and will be overwritten in the sample sheet")

	## Experiments at the level of PCR replicates 
	s = samplesheet[keyCols + repCols + ['PCRRep']].drop_duplicates()
	s['ExperimentIDPCRRep'] = ['-'.join([str(v) for v in list(row[keyCols])]) + "-Rep" + '-'.join([str(v) for v in list(row[repCols])]) + "-PCR" + row['PCRRep'] for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	## Experiments at the level of specified replicate columns
	s = samplesheet[keyCols + repCols].drop_duplicates()
	s['ExperimentIDReplicates'] = ['-'.join([str(v) for v in list(row[keyCols])]) + "-Rep" + '-'.join([str(v) for v in list(row[repCols])]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	## Experiments at the level of experiments (combined across replicates)
	s = samplesheet[keyCols].drop_duplicates()
	s['ExperimentID'] = ['-'.join([str(v) for v in list(row.values)]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	return(samplesheet)	


def add_outputs(samplesheet):
	samplesheet['ExperimentIDPCRRep_BinCounts'] = ['results/byPCRRep/{}.bin_counts.txt'.format(e) for e in samplesheet['ExperimentIDPCRRep']]
	samplesheet['ExperimentIDReplicates_BinCounts'] = ['results/byExperimentRep/{}.bin_counts.txt'.format(e) for e in samplesheet['ExperimentIDReplicates']]
	return samplesheet


def validate_sample_sheet(samplesheet):
	print("Validating the Sample Sheet ...\n")

	for col in requiredCols:
		if not col in samplesheet.columns:
			raise ValueError("Missing required column in sample sheet: " + col)

	if not samplesheet['SampleID'].is_unique:
		raise ValueError("SampleID column in samplesheet must not contain duplicates.")

	for col in keyCols:
		if not col in samplesheet.columns:
			raise ValueError("Missing column in sample sheet that is provided in experiment_keycols in the config file: " + col)

	for col in repCols:
		if not col in samplesheet.columns:
			raise ValueError("Missing column in sample sheet that is provided in replicate_keycols in the config file: " + col)

	if any(samplesheet['Bin'] == "Water"):
		raise ValueError("Found 'Water' value in Bin column. Please instead enter a blank in this slot.")

	print("Found all Experiment Key columns. Generating comparisons for the following experiments:")
	print('\t'.join(keyCols))
	for index, row in samplesheet[keyCols].drop_duplicates().iterrows():
		print('\t'.join([str(v) for v in row.values]))


def load_sample_sheet(samplesheetFile):
	samplesheet = pd.read_table(samplesheetFile, dtype=str)
	validate_sample_sheet(samplesheet)
	samplesheet.index = samplesheet['SampleID']  ## Requires that SampleID is unique
	samplesheet = find_fastq_files(samplesheet, fastqdir)
	samplesheet = add_experiment_names(samplesheet)
	samplesheet = add_outputs(samplesheet)
	samplesheet.index = samplesheet['SampleID']  ## Repeat, because this seems to be lost in merge steps above
	return samplesheet


def get_bin_list():
	binList = samplesheet['Bin'].drop_duplicates()
	binList = binList[(binList != "All") & (binList != "Neg") & (binList.notnull())]
	binList = [str(b) for b in list(binList)]
	print("Processing unique bins: " + ' '.join(binList))
	return(binList)

def get_gene(wildcards):
	currSamples = samplesheet.loc[(samplesheet['ExperimentIDReplicates'] == wildcards.ExperimentIDReplicates) & (samplesheet['Bin'].isin(binList))]
	Gene = currSamples['FlowFISHGene'].unique()
	if (len(Gene) != 1):
		print(currSamples['SampleID'])
		raise ValueError("Found more than one possible Gene. Correct the samplesheet and rerun.")
	return Gene[0]



# global variables
requiredCols = ['SampleID','Bin','PCRRep']
keyCols = config['experiment_keycols'].split(',')
repCols = config['replicate_keycols'].split(',')
codedir = config['codedir']
sortparamsdir = config['sortparamsdir']
fastqdir = config['fastqdir']

# defining some other globals here
window = 10 # had to change this because not enough guides at prom
mineffectsize = 0
minGuides = 3
fdr = 0.05

samplesheet = load_sample_sheet(config['sample_sheet'])
samplesheet.to_csv("SampleList.snakemake.tsv", index=False, header=True, sep='\t')
binList = get_bin_list()



#######################################################################################
####### helpers ###########

def all_input(wildcards):

	wanted_input = []

	## Guide counts per sample:
	wanted_input.extend(['results/counts/{SampleID}.count.txt'.format(SampleID=s) for s in samplesheet['SampleID']])
	
	## Output files for PCR replicates
	wanted_input.extend(list(samplesheet['ExperimentIDPCRRep_BinCounts'].unique()))
	wanted_input.append("results/summary/GuideCounts.flat.tsv.gz")

	## Output files for replicate experiments
	if len(repCols) > 0:
		wanted_input.extend(list(samplesheet['ExperimentIDReplicates_BinCounts'].unique()))
		wanted_input.extend([
			'results/byExperimentRep/{}.raw_effects.txt'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])
		wanted_input.extend([
			'results/byExperimentRep/{}.scaled.txt'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])
		wanted_input.extend([
			'results/byExperimentRep/{}.KnownEnhancers.FlowFISH.txt'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])

	## Output files for experiments (with replicates merged)
	wanted_input.extend(["results/summary/PeakCallingSummary.tsv"])

	return wanted_input


