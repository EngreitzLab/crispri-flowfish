# take in mle outputs in log space
# convert back to real space
# normalize to negative controls
# clamp guides
# filter guides

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Converts MLE log outputs to real space, normalizes to NCs, and filters/clamps')
parser.add_argument("-m", "--mleouts", dest="mle", type=str,  help="MLE outputs (in log space)")
parser.add_argument("-d", "--design", dest="design", type=str,  help="Design file")
#parser.add_argument("-a", "--bed", dest="bed", type=str,  help="Bed file to write to")
parser.add_argument("-b", "--bedgraph", dest="bedgraph", type=str, help="Bedgraph file to write to")
parser.add_argument("-o", "--clampout", dest="outfile", type=str, help="Clamped bed file to write to")
parser.add_argument("-l", "--log", dest="log", type=str, help="Log file")
parser.add_argument("-c", "--clamp", dest="clamp", type=float, default=5.0, help="Max effect in real space")
parser.add_argument("-n", "--background", dest="background", type=float, default=0, help="Background subtraction number")

parser.add_argument("-g", "--maxg", dest="maxg", type=int, default=10, help="Max number of Gs in guide sequence")
parser.add_argument("-y", "--maxlen", dest="maxlen", type=int, default=21, help="Max guide length")
parser.add_argument("-z", "--minlen", dest="minlen", type=int, default=20, help="Min guide length")
parser.add_argument("-t", "--minofftarget", dest="minofftarget", type=int, default=50, help="Min off target score")
parser.add_argument("-s", "--minsum", dest="minsum", type=int, default=0, help="Min sum of cells across bins")

args = parser.parse_args()

with open(args.log, 'w') as log:
    # read in mle outs
    mle = pd.read_table(args.mle)

    # join with design file
    design = pd.read_table(args.design)

    mle = mle.merge(design, on='OligoID')

    log.write('Input guides: {}\n'.format(len(mle)))

    # convert to real space
    mle['mleAvg'] = np.power(10, mle['logMean'] + np.square(mle['logSD'])/2)

    ### WHAT TO DO ABOUT NORMALIZATION OF VARIANCE
    mle['mleSD'] = np.sqrt((np.power(10, np.square(mle['logSD']))-1)*np.square(mle['mleAvg']))

    # subtract background
    mle['mleAvg'] = mle['mleAvg'] - args.background      

    # normalize to NCs
    neg_ctrls = mle.loc[mle['target'] == 'negative_control']
    log.write('Negative control guides: {}\n'.format(len(neg_ctrls)))
    print('Negative control guides: {}\n'.format(len(neg_ctrls)))
    neg_ctrl_mean = np.mean(neg_ctrls['mleAvg'])
    print(neg_ctrl_mean)
    print(np.mean(mle['mleAvg']))

    # normalize to NCs (and clip at 0)
    # mle['mleAvg'] = (mle['mleAvg'] / neg_ctrl_mean).clip(lower=0)
    mle['mleAvg'] = mle['mleAvg'] / neg_ctrl_mean
    mle['WeightedAvg'] = mle['WeightedAvg'] / np.mean(neg_ctrls['WeightedAvg'])

    # write bed and bedgraph
    bed_cols = ["chr", "start", "end", "name", "score", "strand", "GuideSequence", "target", "OffTargetScore", "OligoID", "WeightedAvg", "mleAvg", "mleSD", "sum1"]
    bedgraph_cols = ["chr", "start", "end", "mleAvg"]

    #mle[bed_cols].to_csv(args.bed, sep='\t', index=False)
    for_bedgraph = mle.loc[~mle['chr'].isna(),bedgraph_cols]
    for_bedgraph['start'] = for_bedgraph['start'].astype(np.int32)
    for_bedgraph['end'] = for_bedgraph['end'].astype(np.int32)
    for_bedgraph.to_csv(args.bedgraph, sep='\t', index=False, header=False)

    for_bedgraph2 = mle.loc[~mle['chr'].isna(),["chr", "start", "end", "WeightedAvg"]]
    for_bedgraph2['start'] = for_bedgraph2['start'].astype(np.int32)
    for_bedgraph2['end'] = for_bedgraph2['end'].astype(np.int32)
    for_bedgraph2.to_csv(args.bedgraph + ".WeightedAvg.bedgraph", sep='\t', index=False, header=False)

    # clamp
    #mle.loc[mle['mleAvg'] > args.clamp, 'mleAvg'] = args.clamp

    print(sum(mle['target'] == 'negative_control'))

    # filter guides
    def is_valid_guide(row):
        return len(row['GuideSequence']) >= args.minlen and len(row['GuideSequence']) <= args.maxlen and row['OffTargetScore'] >= args.minofftarget and row['GuideSequence'].count('G') <= args.maxg and row['sum1'] >= args.minsum
        # return len(row['GuideSequence']) >= args.minlen and len(row['GuideSequence']) <= args.maxlen and row['GuideSequence'].count('G') <= args.maxg and row['sum1'] >= args.minsum


    filtered = mle.loc[mle.apply(lambda row: is_valid_guide(row), axis=1)]
    print(sum(filtered['target'] == 'negative_control'))

    log.write('Guides passing filter: {}\n'.format(len(filtered)))

    # write outputs
    filtered[bed_cols].to_csv(args.outfile, sep='\t', index=False)
