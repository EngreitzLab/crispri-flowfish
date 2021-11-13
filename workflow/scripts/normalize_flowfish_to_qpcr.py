import pandas as pd
import argparse
import numpy as np


def flowfish_to_qpcr(info_file, in_file, out_file, bedgraph, clamp, columns='mleAvg mleSD'.split()):
    assert in_file != out_file, "Can't save results of normalization in same file"

    # load screen info
    screen_info = pd.read_table(info_file).iloc[0]

    # compute transform from ff to qpcr
    lo_ff = screen_info.FlowFISH_at_TSS
    lo_qpcr = screen_info.TSS_qPCR
    if (np.isnan(lo_qpcr) or lo_qpcr == "NA"): ## If this value is not provided, then set this value so effectively no adjustment is performed
        print("qPCR TSS value missing, so we will not transform the data.\n")
        lo_qpcr = lo_ff
    assert (lo_qpcr >= 0 and lo_qpcr < 1)
    slope = (1 - lo_qpcr) / (1 - lo_ff)
    print(slope, lo_ff, lo_qpcr)

    # load data to be transformed
    in_data = pd.read_table(in_file)

    # transform
    out_data = in_data.copy()
    for c in columns:
        vals = out_data[c]
        if c.startswith('mleAvg'):
            vals = slope * (vals - lo_ff) + lo_qpcr
            vals = vals.clip(lower=0, upper=clamp)
        elif c.startswith('mleSD'):
            vals = slope * vals
        else:
            assert False, "Don't know how to normalize column  {}".format(c)
        out_data[c] = vals

    # fix coordinates for printing
    out_data['start'] = [str(int(x)) if not np.isnan(x) else '' for x in out_data['start']]
    out_data['end'] = [str(int(x)) if not np.isnan(x) else '' for x in out_data['end']]

    # save output
    out_data.to_csv(out_file, sep="\t", index=None)

    # also write bedgraph
    bedgraph_cols = ["chr", "start", "end", "mleAvg"]
    out_data.loc[~out_data['chr'].isna(),bedgraph_cols].to_csv(bedgraph, sep='\t', index=False, header=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Normalizes FlowFISH scores to qPRC scores')
    parser.add_argument("-i", "--infoFile", dest="INFO", type=str,  help="ScreenInfo.txt file")
    parser.add_argument("-w", "--windowsFile", dest="WINDOWS", type=str,  help="Windows file with mean & sem to be normalized")
    parser.add_argument("-o", "--outFile", dest="OUTFILE", type=str,  help="File for output.  Will have same columns as input")
    parser.add_argument("-b", "--bedgraph", dest="bedgraph", type=str,  help="File for output bedgraph (minimal columns)")
    parser.add_argument("-c", "--clamp", dest="clamp", type=float, help="max value in real space %change for guide", default=5.)

    args = parser.parse_args()

    flowfish_to_qpcr(args.INFO, args.WINDOWS, args.OUTFILE, args.bedgraph, args.clamp)
