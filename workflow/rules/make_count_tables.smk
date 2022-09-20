## Create a count table from crispresso output that looks like:
## MappingSequence  A   B   C   D  ... other bin names

import os
import pandas as pd
import numpy as np


def make_count_table(samplesheet, group_col, group_id, bins, outfile_raw, outfile_frequencies):
    ## Function to make a count table at various layers of resolution (e.g., by experiment, or by replicate, or by PCR replicate)
    ## To do: Move the python code for these rules into separate python scripts so they can be run independently of the snakemake pipeline (at least, this makes it easier to test and debug the code)

    currSamples = samplesheet.loc[samplesheet[group_col]==group_id]

    count_tbls = []

    for idx, row in currSamples.iterrows():
        s=row['SampleID']
        file = "results/counts/{SampleID}.count.txt".format(SampleID=s)

        if (os.path.exists(file)):
            curr = pd.read_table(file, names=[s,'OligoID'])
            curr[s] = curr[s].astype(np.int32)
            count_tbls.append(curr)

    if len(count_tbls) > 0:
        count_tbl = count_tbls.pop()

        for tbl in count_tbls:
            count_tbl = count_tbl.merge(tbl, on='OligoID', how='outer')

        count_tbl = count_tbl.set_index('OligoID')
        count_tbl = count_tbl.fillna(0)
        count_tbl = count_tbl.astype(int)

        ## Now, sum counts per bin
        bin_list = bins + list(set(currSamples['Bin'].unique())-set(bins))
        for uniqBin in bin_list:
            samples = currSamples.loc[currSamples['Bin'] == uniqBin]
            if len(samples) > 0:
                count_tbl[uniqBin] = count_tbl[samples['SampleID']].sum(axis=1).values
            else:
                count_tbl[uniqBin] = 0
        count_tbl = count_tbl[bin_list]

    else:
        count_tbl = pd.DataFrame({'OligoID':[]})
        for uniqBin in bins:
            count_tbl[uniqBin] = []
        count_tbl = count_tbl.set_index('OligoID')

    ## Make compatible with multiple sorters
    ## drop bins from bin list that are not used for particular sorters

    count_tbl = count_tbl.loc[:, (count_tbl.sum(axis=0) != 0)]

    count_tbl.index.name = "OligoID"
    count_tbl.to_csv(outfile_raw, sep='\t')

    freq_tbl = count_tbl.div(count_tbl.sum(axis=0), axis=1)
    freq_tbl.to_csv(outfile_frequencies, sep='\t', float_format='%.6f')



def make_flat_table(samplesheet, outfile):
    count_tbls = []
    for idx, row in samplesheet.iterrows():
        file = "results/counts/{SampleID}.count.txt".format(SampleID=row['SampleID'])

        if (os.path.exists(file)):
            curr = pd.read_table(file, names=['count','OligoID'])
            curr['SampleID'] = row['SampleID']
            curr['count'] = curr['count'].astype(np.int32)
            count_tbls.append(curr)

    flat = pd.concat(count_tbls, axis='index', ignore_index=True)
    flat.to_csv(outfile, sep='\t', index=False, compression='gzip')



rule make_count_table_per_PCRrep:
    input: 
        lambda wildcards: 
            ['results/counts/{SampleID}.count.txt'.format(SampleID=s) 
                for s in samplesheet.loc[samplesheet['ExperimentIDPCRRep']==wildcards.ExperimentIDPCRRep]['SampleID']]
    output:
        counts='results/byPCRRep/{ExperimentIDPCRRep}.bin_counts.txt',
        freq='results/byPCRRep/{ExperimentIDPCRRep}.bin_freq.txt'
    run:
        make_count_table(samplesheet, 'ExperimentIDPCRRep', wildcards.ExperimentIDPCRRep, get_bin_list(), output.counts, output.freq)



rule make_count_table_per_experimentalRep:
    input: 
        lambda wildcards: 
            ['results/counts/{SampleID}.count.txt'.format(SampleID=s) 
                for s in samplesheet.loc[samplesheet['ExperimentIDReplicates']==wildcards.ExperimentIDReplicates]['SampleID']]
    output:
        counts='results/byExperimentRep/{ExperimentIDReplicates}.bin_counts.txt',
        freq='results/byExperimentRep/{ExperimentIDReplicates}.bin_freq.txt'
    run:
        make_count_table(samplesheet, 'ExperimentIDReplicates', wildcards.ExperimentIDReplicates, get_bin_list(), output.counts, output.freq)



rule make_flat_count_table_PCRrep:
    input: 
        lambda wildcards: ["results/counts/{SampleID}.count.txt".format(SampleID=s) for s in samplesheet['SampleID']]
    output:
        'results/summary/GuideCounts.flat.tsv.gz'
    run:
        make_flat_table(samplesheet, output[0])


rule plot_pcr_replicate_correlations:
    input:
        lambda wildcards:
            ['results/byPCRRep/{ExperimentIDPCRRep}.bin_counts.txt'.format(ExperimentIDPCRRep=e) for e in samplesheet['ExperimentIDPCRRep']]
    output:
        pdf='results/summary/ReplicateCorrelation.pdf'
    shell:
        "Rscript crispri-flowfish/workflow/scripts/plotReplicateCorrelations.R \
          --samplesheet SampleList.snakemake.tsv \
          --groupcol ExperimentIDReplicates \
          --filecol ExperimentIDPCRRep_BinCounts \
          --outsummary {output.pdf}"


