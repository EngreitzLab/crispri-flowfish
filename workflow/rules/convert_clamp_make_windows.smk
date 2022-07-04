## Convert guide effects to real space, 

import os
import pandas as pd
import numpy as np


# convert the log-effect sizes into real space and normalize by negative controls
rule convert_to_real_space:
    input:
        mleouts='results/byExperimentRep/{ExperimentIDReplicates}.raw_effects.txt',
        design=config['design_file']
    params:
        minsum=.050, # config['min_guide_count_sum']
        minbins=4
    output:
        #bed='results/byExperimentRep/{ExperimentIDReplicates}.real_space.bed',
        bedgraph='results/byExperimentRep/{ExperimentIDReplicates}.real_space.bedgraph',
        real='results/byExperimentRep/{ExperimentIDReplicates}.real_space.txt',
        log='results/byExperimentRep/{ExperimentIDReplicates}.convert.log'
    # params:
    #     background=lambda wildcards: float(background.loc[wildcards.screen].values[0])
    shell:
        "python crispri-flowfish/workflow/scripts/convert_to_real_space.py \
            -m {input.mleouts} \
            -d {input.design} \
            -b {output.bedgraph} \
            -o {output.real} \
            -l {output.log} \
            --minsum {params.minsum} \
            --minbins {params.minbins}"


# calculate effect per 20-guide window
rule make_windows:
    input:
        real='results/byExperimentRep/{ExperimentIDReplicates}.real_space.txt'
    output:
        window='results/byExperimentRep/{ExperimentIDReplicates}.windows.txt'
    params:
        window=window,
        epcode=config['ep_code_dir']
    shell:
        "Rscript crispri-flowfish/workflow/scripts/CalculateTilingStatistic.R --input {input.real} \
        	--output {output.window} --grnaWindowSize {params.window} --scoreColumn mleAvg --minOffTargetScore 50 \
        	--maxSpan 750 --minT0count 0 --T0columns sum1 --codeDir {params.epcode}"

# calcuate the knockdown at the TSS
rule tss:
    input:
        window='results/byExperimentRep/{ExperimentIDReplicates}.windows.txt'
    output:
        info='results/byExperimentRep/{ExperimentIDReplicates}.ScreenInfo.txt'
    params:
        qpcr=config['qpcr'],
        genelist=config['genelist'],
        gene=get_gene_for_qPCR_adjustment
    shell:
        "python crispri-flowfish/workflow/scripts/FlowFISHtssKD.py -g {params.gene} -w {input.window} -q {params.qpcr} -l {params.genelist} -o {output.info}"

# re-scale effects based on the TSS levels
rule rescale:
    input:
        info='results/byExperimentRep/{ExperimentIDReplicates}.ScreenInfo.txt',
        clamp='results/byExperimentRep/{ExperimentIDReplicates}.real_space.txt'
    output:
        scaled='results/byExperimentRep/{ExperimentIDReplicates}.scaled.txt',
        bedgraph='results/byExperimentRep/{ExperimentIDReplicates}.scaled.bedgraph',
    shell:
        "python crispri-flowfish/workflow/scripts/normalize_flowfish_to_qpcr.py -i {input.info} -w {input.clamp} -o {output.scaled} -b {output.bedgraph}"



