# Snakemake workflow: CRISPRi-FlowFISH

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}}.svg?branch=master)](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}})

This snakemake workflow is for analysis of CRISPRi-FlowFISH data.


## Authors

* Ben Doughty (@bdoughty)
* Jesse Engreitz (@engreitz)

## Description

To do

## Usage

### Step 1: Clone this github repository

[Clone](https://help.github.com/en/articles/cloning-a-repository) this to your local system, into the place where you want to perform the data analysis.

### Step 2: Install conda environment

Install Snakemake and conda environment using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda env create --file envs/EngreitzLab.yml [TODO: replace with pipeline-specific conda environment]

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 3.1: Set up the Sample Sheet 

(Updated 5/7/21)

The Sample Sheet lists all of the sequencing libraries that will be included in the analysis, and describes their relationships and groupings.

Required columns:
    
    SampleID          Unique name for each amplicon library.
    Bin               Name of a FACS-sorted bin (e.g.: A B C D E F). 'All' for input samples. 'Neg' or blank if not applicable
    PCRRep            PCR replicate number or name

Optional columns:

    [Experiment Keys] Provide any number of additional columns (e.g., CellLine) that distinguish different samples.
                        Key columns are defined as such by the 'experiment_keycols' parameter in the config file.
                        These columns will be combined to form a unique experiment key.
                        Replicates for a given unique experiment key will be combined.
    [Replicate Keys]  Provide any number of additional columns (e.g., FlowFISHRep) that distinguish different experimental replicates (not including PCR replicates)
                        Replicate columns are defined as such by the 'replicate_keycols' parameter in the config file.
                        These columns will be combined to form a unique replicate id.
                        PCR replicate counts for each unique replicate key will be summed at the level of this replicate ID.
                        MLE estimates will also be performed at the level of this replicate ID, 
                        then compared according to grouping of the experiment key.
    fastqR1           If provided in the Sample Sheet, overwrites the default value (config['fastqdir']/{SampleID}_*_R1_*fastq.gz)

### Step 3.2: Set up the guide design file

This file lists information about the gRNAs included in the experiment, and should in theory be output by the Engreitz Lab CRISPRi design pipeline (https://github.com/broadinstitute/CRISPRiTilingDesign)

Columns (most are legacy from the gRNA designer and are not used by this pipeline):
    
    chr                 chromosome for gRNA spacer (can be blank, e.g. for non-targeting control)
    start               start coordinate for gRNA spacer (can be blank, e.g. for non-targeting control)
    end                 end coordinate for gRNA spacer (can be blank, e.g. for non-targeting control)
    name                unique ID for the gRNA 
    score               A score for the gRNA (not used by this pipeline)
    strand              strand of the gRNA spacer
    GuideSequence       Spacer sequence (variable length)
    GuideSequenceMinusG Spacer sequence, minus a leading G (variable length)
    MappingSequence     Fixed length spacer sequence used for mapping
    OffTargetScore      gRNA off-target score
    target              Required column used for grouping gRNAs into targets (e.g. promoter of a given gene) or 'negative_control' for negative control gRNAs
    subpool             Column specifying a group of gRNAs 
    OligoID             unique ID for the gRNA [required]


### Step 3.3: Set up the TSS qPCR file (optional)

This file specifies an optional scaling parameter, representing the effect on gene expression when perturbing the TSS of a given gene, used to linearly adjust the effects estimated by the FlowFISH assay. For example, some probesets have some detectable background fluorescence that means that the knockdown from FlowFISH appears to be only 50%, when in reality it is 80% knockdown plus some background. To remove this background signal, we pass in an externally measured TSS knockdown value (e.g. 80%) into the pipeline, which is then used to shift/adjust the FlowFISH estimated effects to match (based on the best 20-gRNA window around the TSS). 

If this file is not provided, then this TSS qPCR shifting calculation is not done (see workflow/scripts/FlowFISHtssKD.py)

This file includes the following columns:
    
    [KeyCols]           Any number of columns matching the names of columns in the Sample Sheet, for table join
    name                Gene symbol to merge into the config.genelist file
    TSS_Override        If provided, overrides the TSS coordinate provided in 'tss' column of config.genelist
    

### Step 4: Provide sort params files

To do.

### Step 5: Configure workflow

[TODO]:  For now, edit `workflow/config.json` to point to the right files


### Step 6: Execute workflow

Activate the conda environment:

    conda activate EngreitzLab 
    ## TODO: Create specific environment for thie pipeline and check in yml file to workflows/envs/

Test your configuration by performing a dry-run via

    snakemake -s crispri-flowfish/workflow/Snakefile --configfile config/config.json -n

Execute the workflow locally via

    snakemake -s crispri-flowfish/workflow/Snakefile --configfile config/config.json

using `$N` cores or run it in a cluster environment (Stanford Sherlock SLURM) via

`
snakemake \
  -s crispri-flowfish/workflow/Snakefile \
  --configfile config/config.json \
  --cores 1 \
  --jobs 200 \
  --cluster "sbatch -n 1 -c 1 --mem 8G -t 4:00:00 -p owners -J VFF_{rule} -o log/{rule}_{wildcards} -e log/{rule}_{wildcards}"
`

For more about cluster configuration using snakemake, see [here](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/)

### Step 7: Interpreting outputs

Description of output files, in sequential order (updated JME 11/13/21):

    {s}.bin_counts.txt              Guide counts per FlowFISH bin
    {s}.bin_freq.txt                Guide count frequencies per FlowFISH bin (counts for gRNA G / sum of counts for all gRNAs)
    {s}.raw_effects.txt             Two estimates of gRNA expression levels, output by estimate_effect_sizes.R:  WeightedAvg represents weighted average expression across bins, and logMean and logSD give the MLE estimate (in log10 space)
    {s}.real_space.{txt,bedgraph}   gRNA effects converted to "fraction expression remaining" in real space, normalized to negative control gRNAs (e.g. 1 = gRNA does not change expression vs control gRNAs; 0.4 = gRNA reduces expression by 60% vs control gRNAs)
    {s}.windows.*                   Average effects of 20-gRNA windows (max span: 750bp). Most relevant for comprehensive tiling screens, not as relevant for peak-focused screens except it is used for the qPCR adjustment.
    {s}.scaled.*                    gRNA effects ("fraction expression remaining") linearly shifted so that knockdown at the TSS of the target gene matches the provided qPCR data, and clamped to limit the maximum gRNA effect size (currently, to 5) [To do: Add flag to skip this adjustment for making downstream files]
    {s}.collapse.bed                BED file containing peak coordinates and the scores of overlapping RNAs (comma-separated list)
    {s}.FullEnhancerScore.txt       Table containing full statistics for comparing gRNAs in each peak to negative controls (effect size for each peak, p-value, n guides, etc.)
    {s}.PeakCallingSummary.txt      Summary of peak calling statistics for a sample, to facilitate comparisons across samples
    {s}.tfdr0.05.bed                BED file containing significant peaks (where t-test adjusted p < 0.05 for comparing gRNAs in the peak vs negative control gRNAs)
    {s}.ScreenData.txt              Formatted metadata for creating the KnownEnhancers.txt formatted file
    {s}.KnownEnhancers.FlowFISH.txt Format screen data for integrating across samples for comparison to models, including annotating peaks with gene TSS coordinate, distance to TSS, study name, etc.
    
Recommended analysis:
* Load in the .real_space.bedgraph files into IGV to look at MLE effect size estimates for individual gRNAs, before qPCR scaling and effect size clamping
* Use .FullEnhancerScore.txt (or .KnownEnhancers.FlowFISH.txt file, if created) for analysis of results aggregated by peaks [note that this includes the qPCR scaling]

