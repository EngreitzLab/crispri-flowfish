from __future__ import division


# Usage: python ScoreEnhancers.py SCREENPREFIX PROJECTDIR POWER(True or False) MINEFFECTSIZE
#                                (i.e. HDAC6-170303) REQUIRES a HDAC6-170303.collapse.bed and HDAC6-170303.mean.bed files

# import stuff
import pandas as pd
import numpy as np
from pandas.io.parsers import read_table
import statsmodels.stats.multitest
import random
from scipy import stats
import sys
import argparse

def writeBed(df, OUTFILE, score=None,negative=False):
    a=df.copy()
    try:
        a["start"]=np.vectorize(str)(np.vectorize(int)(a["start"]))
        a["end"]=np.vectorize(str)(np.vectorize(int)(a["end"]))
        
        cols=["chr", "start", "end"]
        if score:
            cols=["chr", "start", "end", score]
        if negative:
            a[score]=-a[score]
        a[cols].to_csv(OUTFILE, sep='\t', header=False, index=False)
    except KeyError:
        a.to_csv(OUTFILE, sep='\t', header=False, index=False)


##############################################################

parser = argparse.ArgumentParser(description='Score Enhancers')
parser.add_argument("-c", "--collapsed", dest="collapsed", type=str,  help="Collapsed file")
parser.add_argument("-s", "--scaled", dest="scaled", type=str,  help="Scaled file")
parser.add_argument("-i", "--info", dest="info", type=str,  help="Screen info", default=None)
parser.add_argument("-u", "--summary", dest="summary", type=str,  help="Peak calling summary")
parser.add_argument("-o", "--fullscore", dest="fullscore", type=str,  help="Full enhancer score")
parser.add_argument("-p", "--peaks", dest="peaks", type=str,  help="Peaks")
parser.add_argument("-n", "--exptName", dest="exptName", type=str, help="Name of experiment to include in the peak calling summary")

parser.add_argument("-m", "--minguides", dest="minguides", type=int, default=20, help="Min guides")
parser.add_argument("-f", "--fdr", dest="fdr", type=float, default=0.05, help="FDR threshold")
parser.add_argument("-e", "--mineffect", dest="mineffect", type=float, default=0., help="Min effect size")

args = parser.parse_args()


# SCREEN=sys.argv[1]
# PROJECTDIR=sys.argv[2]
# POWER=sys.argv[3]

# FDRCUTTOFF=0.05 # MYC used 0.032 in original paper
# MINGUIDES=1 #20
# MINEFFECTSIZE=0
# SCORECOL="mleAvg"

# if len(sys.argv)>4:
#     MINEFFECTSIZE=float(sys.argv[4])
#print MINEFFECTSIZE, "Min Effect Size"

COLLAPSE = args.collapsed
ULTIMATEOUT = args.scaled
SCREENINFO = args.info
SUMOUT = args.summary
OUTFILE = args.fullscore
PEAKBED = args.peaks

MINGUIDES = args.minguides
FDRCUTTOFF = args.fdr
MINEFFECTSIZE = args.mineffect
SCORECOL = "mleAvg"

SCREEN = args.exptName


# COLLAPSE=PROJECTDIR+"/"+SCREEN+"/"+SCREEN+".collapse.bed"
# # ULTIMATEOUT=PROJECTDIR+"/"+SCREEN+"/"+SCREEN+".mean.bed"
# ULTIMATEOUT=PROJECTDIR+"/"+SCREEN+"/"+SCREEN+".mean.scaled.txt"
# SCREENINFO=PROJECTDIR+"/"+SCREEN+"/"+SCREEN+".ScreenInfo.txt"
# SUMOUT=PROJECTDIR+"/"+SCREEN+"/"+SCREEN+".PeakCallingSummary.txt"

# OUTFILE=PROJECTDIR+"/"+SCREEN+"/"+SCREEN+".FullEnhancerScore.txt"
# PEAKBED=PROJECTDIR+"/"+SCREEN+"/"+SCREEN+".peaks.fdr"+str(FDRCUTTOFF)+".bed"
# UPEAKBED=PROJECTDIR+"/"+SCREEN+"/"+SCREEN+".peaks.ufdr"+str(FDRCUTTOFF)+".bed"
# TPEAKBED=PROJECTDIR+"/"+SCREEN+"/"+SCREEN+".peaks.tfdr"+str(FDRCUTTOFF)+".bed"

header=["chr", "start", "end", "scores"]
data=read_table(COLLAPSE, names=header)


#print SCREEN, "Screen"
#print len(data), "Enhancers"

# convert scores to a list of floats from a string
data["scores"]=data["scores"].apply(lambda y: [float(x) for x in y.split(",")])

ultimateOut=read_table(ULTIMATEOUT)
negative_controls=ultimateOut[ultimateOut["target"]=="negative_control"]

#print len(ultimateOut), "Guides"
#print len(negative_controls), "negative_controls"

data["n"]=data["scores"].apply(lambda x: len(x))
data=data[data["n"]>=MINGUIDES]

data["mean"]=data["scores"].apply(lambda x: np.mean(x))
data["p.utest"]=data["scores"].apply(lambda x: stats.mannwhitneyu(x, negative_controls["mleAvg"]) [1])
data["fdr.utest"]=statsmodels.stats.multitest.fdrcorrection(data["p.utest"])[1]
data["p.ttest"]=data["scores"].apply(lambda x: stats.ttest_ind(x, negative_controls["mleAvg"], equal_var=True)[1])
data["fdr.ttest"]=statsmodels.stats.multitest.fdrcorrection(data["p.ttest"])[1]
data["sem"]=data["scores"].apply(lambda x: stats.sem(x))

# some screens have rare very large effectguides. Keep these in when significance testing, but remove them for showing the mean
#data["mean"]=data["scores"].apply(lambda x: stats.trim_mean(x, proportiontocut=0.05))

# Trimmed mean helps, but not enough. (PIM2-1 is better than PIM2-2, but still the worst by far...) ? try median

NC_mean=np.mean(negative_controls["mleAvg"])
#print NC_mean, "NC mean"

# Regulated is changing and putting gene down when CRISPRi'ed
reg=[False]*len(data)
sigu=[False]*len(data)
sigt=[False]*len(data)
for row in range(len(data)):
    tRow=data.iloc[row]
    ES=abs(tRow["mean"]-NC_mean)
    if (tRow["fdr.utest"]<FDRCUTTOFF) and (ES>=MINEFFECTSIZE):
        sigu[row]=True
        if tRow["mean"]<NC_mean:
            reg[row]=True
    if (tRow["fdr.ttest"]<FDRCUTTOFF) and (ES>=MINEFFECTSIZE):
        sigt[row]=True

#print sum(sigu), "Significant Elements (U-test)"
#print sum(sigt), "Significant Elements (T-test)"

data["Regulated"]=reg
#data["Significant"]=sigu
data["Significant"]=sigt # fixed 180507 1800

# add in control info 
data["mean.ctrl"]=NC_mean
data["n.ctrl"]=len(negative_controls)
data["sem.ctrl"]=stats.sem(negative_controls["mleAvg"])

## JE 5/5/18:  Writing these columns is unnecessary and will be overwritten by later scripts in the pipeline
data["CustomTargetGeneTSS"]=""

# write full data ##TASK: Add in all the columns to be used in model
# chr   start   end mean    sem n   mean.ctrl   sem.ctrl    n.ctrl  p.utest p.ttest fdr.utest   fdr.ttest
# the FormatCRISPR script will require 

# just write the stuff with enough guides to be informative (>14?)
## JE 5/5/18:  Writing the Regulated column is unnecessary and will be overwritten in later scripts in the pipeline
cols=["chr", "start", "end", "mean", "sem","n","mean.ctrl", "sem.ctrl","n.ctrl", "p.utest", "fdr.utest", "p.ttest", "fdr.ttest", "CustomTargetGeneTSS", "Significant", "Regulated"]
tw=data[data["n"]>=MINGUIDES]

tw[cols].to_csv(OUTFILE, sep="\t", index=False, header=True)

# write bed of peaks
peaks=tw[tw["Significant"]]
print(peaks)
writeBed(peaks, PEAKBED)
# writeBed(peaks, UPEAKBED)
# writeBed(tw[sigt], TPEAKBED)

## 180507: write some summary numbers 
summary=pd.DataFrame()
summary["Screen"]=[SCREEN]
summary["GuidePerDHS.min"]=[min(tw['n'])]
summary["GuidePerDHS.max"]=[max(tw['n'])]
summary["GuidePerDHS.median"]=[np.median(tw['n'])]
summary["GuidePerDHS.MINGUIDES"]=[MINGUIDES]
summary["mean.ctrl"]=[NC_mean]
summary["n.ctrl"]=[len(negative_controls)]
summary["sem.ctrl"]=[stats.sem(negative_controls["mleAvg"])]
summary["sig.elements.t"]=[sum(sigt)]
summary["all.elements"]=[len(tw)]

assert len(summary)==1
summary.to_csv(SUMOUT, sep="\t", index=False, header=True)

# ------------------------------------------------------------------- #
# POWER
if False: # POWER:
    import random
    POWEROUT=PROJECTDIR+"/"+SCREEN+"/"+SCREEN+".BH-Power-fdr"+str(FDRCUTTOFF)+".txt"
    def ff2reality(ES,BG):
            '''Rescales a given effect size based on assumption of constant background'''
            #return (BG+ES)/(BG+1)
            return 1-(BG+1-ES)/(BG+1)   
    
    info=read_table(SCREENINFO)
    assert len(info)==1, "Each Screen must be present exactly once"
    BG=info["Background"].iloc[0]
    
    FakeEnhancerN=int(data.median()["n"])
    eslist=[0.05, 0.10, 0.15, 0.2, 0.25, 0.3]
    power=[]
    for ES in eslist:
    
        # adjust ES based on linear background assumption
            ES=ff2reality(ES,BG)
    
            real=negative_controls[SCORECOL]
            shift=real-ES # note that this can create negative NC numbers...
            # is it better to cap at 0?

            results=[]
            reps=400
            results=[0]*reps
            for i in range(0, reps):
                sample=random.sample(shift.tolist(), FakeEnhancerN)

                # make the sample a dataframe row
                exp=pd.DataFrame(data={'chr': ["chrNC"], 'start': [0], "end": [0], "scores": [sample], "n":FakeEnhancerN})

                # cat with the real data
                exp=pd.concat([exp, data])

                exp["mean"]=exp["scores"].apply(lambda x: np.mean(x))
                exp["p.ttest"]=exp["scores"].apply(lambda x: stats.ttest_ind(x, negative_controls["mleAvg"], equal_var=True)[1])
                exp["fdr.ttest"]=statsmodels.stats.multitest.fdrcorrection(exp["p.ttest"])[1]
                
            # prior to 180511 used utest!
            exp["p.utest"]=exp["scores"].apply(lambda x: stats.mannwhitneyu(x, negative_controls["mleAvg"])[1])
            exp["fdr.utest"]=statsmodels.stats.multitest.fdrcorrection(exp["p.utest"])[1]

            tested=exp.iloc[0] # the fake enhancer made from NCs is in index 0
            assert tested["chr"]=="chrNC"

            result=tested["fdr.ttest"]<FDRCUTTOFF # prior to 180511 used utest!
            results[i]=(result)
            power.append(sum(results)/reps)

    pd.DataFrame(data={"EffectSize":eslist, "Power":power}).to_csv(POWEROUT, sep="\t", index=False, header=True)
