#! /usr/bin/env python

import pandas as pd
import argparse
from pandas.io.parsers import read_table

#######################################################################################################################################################

def FlowFISHtssKD(GENE, WINDOWS, QPCR, GENELIST, OUTFILE, SLOP):
    
    #----------------------------------------- #
    # read input files
    #----------------------------------------- #
    windows=read_table(WINDOWS)
    SCREEN=GENE

    
    qpcr=read_table(QPCR)

    if '-' in SCREEN:
        # treat as normal screen rep
        qpcr=qpcr[qpcr["Screen"]==SCREEN]# match the qpcr to the screen from windows
    else:
        # treat as a gene
        qpcr=qpcr.loc[qpcr["qPCRGene"]==SCREEN].iloc[[0]]
        # qpcr=qpcr.loc[qpcr["name"]==SCREEN].iloc[[0]]
    #print qpcr.head(2)
    assert len(qpcr)==1, "Screen (e.g. GATA1-1) must appear exactly once in qPCR file: " + SCREEN
    genes=read_table(GENELIST)

    #----------------------------------------- #
    # Get TSS region (+/- 500 kb from TSS)
    #----------------------------------------- #
    data=pd.merge(qpcr, genes[['name','tss','chr','start','end','strand']], on='name')
    
    tss=pd.DataFrame()
    tss["chr"]=data["chr"]
    tss["start"]=data["tss"]-SLOP#.apply(lambda x: x-SLOP)
    tss["end"]=data["tss"]+SLOP#.apply(lambda x: x+SLOP)

    if "TSS_Override" in qpcr:
        if qpcr["TSS_Override"].isnull().any()==False: # If the TSS in the genes file is wrong
            print("TSS override: {}".format(SCREEN))
            tss["start"]=int(qpcr["TSS_Override"]-SLOP)
            tss["end"]=int(qpcr["TSS_Override"]+SLOP)
    print(tss[['start', 'end']])

    print(qpcr["TSS_qPCR"])
        
    #----------------------------------------- #
        # Get 20-guide window
        #----------------------------------------- #
    promoter=windows.copy()
    promoter=promoter[promoter["chr"]==tss["chr"].iloc[0]]
    # promoter=promoter[promoter["start"]>=tss["start"].iloc[0]]
    # promoter=promoter[promoter["end"]<=tss["end"].iloc[0]]
    promoter=promoter[promoter["start"]<tss["end"].iloc[0]]
    promoter=promoter[promoter["end"]>tss["start"].iloc[0]]
    print("20 guide windows in 1kb TSS region: ", len(promoter))
    best=promoter.sort_values(by="mean", ascending=True).iloc[0]

    #----------------------------------------- #
        # Calculate FlowFISH Background
        #----------------------------------------- #
    FF=best["mean"]
    qpcr["FlowFISH_at_TSS"]=1+best["mean"] # w20 gives KD as negative, so 20% remaining is -0.8
    try:
        qpcr["Background"]=(qpcr["FlowFISH_at_TSS"]-qpcr["TSS_qPCR"])/(1-qpcr["FlowFISH_at_TSS"])
    except:
        qpcr["Background"]="NA"

        #----------------------------------------- #
        # Write output
        #----------------------------------------- #
    for col in ["chr", "start", "end"]: qpcr[col]=tss[col].iloc[0]
    qpcr.to_csv(OUTFILE, sep="\t", header=True, index=False)#, quote=False)
    print(qpcr)

# -------------------------------------------------------------------------------------------------------------------------------------------------- #  

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Finds the best 20 guide window within 500 bp of the TSS and calculates FlowFISH Background')
    parser.add_argument("-g", "--gene", dest="GENE", type=str,  help="Gene name to search for ")
    parser.add_argument("-w", "--windowFile", dest="WINDOWS", type=str,  help="File of 20 guide windows")
    parser.add_argument("-q", "--qpcrFile", dest="QPCR", type=str,  help="Tab-delimited table of qPCR data. Screen|name|qPCR|[optional TSSOverride]")
    parser.add_argument("-l", "--genelist", dest="GENELIST", type=str,  help="Neighborhood GeneList.txt")
    parser.add_argument("-o", "--outFile", dest="OUTFILE", type=str,  help="File for output. Line of the qPCR file for this screen with added columns")
    parser.add_argument("-s", "--slop", dest="SLOP", type=int, default=500, help="Slop around annotated TSS to find best 20 guide window. default 500.")

    args = parser.parse_args()

    FlowFISHtssKD(args.GENE, args.WINDOWS, args.QPCR, args.GENELIST, args.OUTFILE, args.SLOP)





