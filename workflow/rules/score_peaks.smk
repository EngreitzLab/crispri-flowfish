## Take guide effects and overlap with DHS peaks, score enhancers

import os
import pandas as pd
import numpy as np


# turn neighborhoods/enhancerlist bed into enhancer file for the flowfish
rule make_enhancers:
    input:
        # "{}/EnhancerList.bed".format(path.join(neighborhoods, cellLine))
        config['enhancer_list']
    output:
        "config/enhancers_from_neighborhoods.bed"
    shell:
        "cat {input} | cut -f1-3 | bedtools sort > {output}"

# collapse scores to peaks
rule collapse:
    input:
        scaled='results/byExperimentRep/{ExperimentIDReplicates}.scaled.txt',
        enhancers="config/enhancers_from_neighborhoods.bed"
    output:
        'results/byExperimentRep/{ExperimentIDReplicates}.collapse.bed'
    shell:
        'tail -n+2 {input.scaled} | grep -vw "negative_control" | awk \'$2!="NA"\' | awk \'$2!="NaN"\' | \
            bedtools sort | bedtools map -o collapse -c 12 -a {input.enhancers} -b stdin | awk \'$4!="."\' > {output}'

# score individual peaks
rule score_dhs:
    input:
        scaled='results/byExperimentRep/{ExperimentIDReplicates}.scaled.txt',
        collapse='results/byExperimentRep/{ExperimentIDReplicates}.collapse.bed'
        # info='{screen}/{screen}.ScreenInfo.txt'
    output:
        score='results/byExperimentRep/{ExperimentIDReplicates}.FullEnhancerScore.txt',
        summary='results/byExperimentRep/{ExperimentIDReplicates}.PeakCallingSummary.txt'        
    params:
        peaks=lambda wildcards: 'results/byExperimentRep/{}.peaks.tfdr{}.bed'.format(wildcards.ExperimentIDReplicates, fdr),
        mineffectsize=mineffectsize,
        minguides=minGuides,
        fdr=fdr
    shell:
        "python crispri-flowfish/workflow/scripts/ScoreEnhancers.py -c {input.collapse} -s {input.scaled} \
            -u {output.summary} -o {output.score} -p {params.peaks} -m {params.minguides} \
            -e {params.mineffectsize} -f {params.fdr} --exptName {wildcards.ExperimentIDReplicates}" 

        # $SCREEN $PROJECTDIR $POWER $MINEFFECTSIZE"
        # "{params.score} {params.screen} {params.mineffectsize} {params.enhancers} . {params.power} {params.codedir}"


# Combine peak calling summary statistics across samples
rule summarize_peak_calling:
    input:
        lambda wildcards:
            ['results/byExperimentRep/{ExperimentIDReplicates}.PeakCallingSummary.txt'.format(ExperimentIDReplicates=e) for e in samplesheet['ExperimentIDReplicates'].drop_duplicates()]
    output:
        summary='results/summary/PeakCallingSummary.tsv'
    shell:
        "csvtk concat -t {input} > {output.summary}"


# format screen for prediction
rule format_data:
    input:
        score='results/byExperimentRep/{ExperimentIDReplicates}.FullEnhancerScore.txt'
    output:
        screendata='results/byExperimentRep/{ExperimentIDReplicates}.ScreenData.txt'
    params:
        gene=get_target_gene_symbol,
        score='results/byExperimentRep/{ExperimentIDReplicates}.scaled.txt' # need to add projectdir for full path here?
    run:
        shell('echo -e "Screen\tScreenData\tRNAReadoutMethod\tReference" > {output.screendata}')
        #shell('echo -e "{params.gene}\t{params.score}\tFlowFISH Screen\tThisStudy" >> {output.screendata}')
        shell('echo -e "{params.gene}\t{input.score}\tFlowFISH Screen\tThisStudy" >> {output.screendata}')
        shell('echo "{params.gene}"')

# format screen for prediction
rule known_enhancers:
    input:
        screendata='results/byExperimentRep/{ExperimentIDReplicates}.ScreenData.txt'
    output:
        screen="results/byExperimentRep/{ExperimentIDReplicates}.KnownEnhancers.FlowFISH.txt"
    params:
        format=os.path.join(config['ep_code_dir'], 'src/CRISPRScreen/FormatCRISPRiScreensForPredictions.R'), ## TODO: Move this script into this repo, or point to github version
        cellLine = config['cell_line'],
        genes=config['genelist'],
        enhancers=config['enhancer_list'],
        epcode=config['ep_code_dir'],
        fdr=fdr
    shell:
        "Rscript {params.format} \
            --output {output.screen} \
            --genesFile {params.genes} \
            --enhancersFile {params.enhancers} \
            --screenData {input.screendata} \
            --cellLine {params.cellLine} \
            --codeDir {params.epcode} \
            --addCRISPRiExtension TRUE \
            --collapseToModelElements TRUE \
            --buffer3p 2000 \
            --fdrCutoff {params.fdr}"

