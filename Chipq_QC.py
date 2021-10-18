from ruffus import *
import sys
import os
import sqlite3
import shutil
import cgatcore.experiment as E
from cgatcore import pipeline as P
import re
import glob
import collections


# Pipeline configuration
P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"],
    defaults={
        'paired_end': False})

PARAMS = P.PARAMS

####merged fil####

@transform (PARAMS["raw_file.dir"], formatter(), r"{basename[0]}_r1.bed")
def bedile(infile,outfile):
    inputs = P.snip (infile,".ped")
    output = P.snip (outfile,".bed")
    faminfile = output+".b"
    ids = output+".f"
    removefam = output +".bim"
    removeid=output +".fam"
    module = PARAMS["Plink_module"]
    statement = '''module load %(module)s && plink --make-bed --file %(inputs)s --out  %(output)s && 
                   awk '$2 in a {$2=$2 "_" ++a[$2]}{a[$2];print}' %(removeid)s > %(ids)s &&
                   rm %(removeid)s && mv %(ids)s %(removeid)s 
                   '''
    P.run(statement)

@transform (bedile, regex(r"(.+)_r1.bed"), r"\1_r2.bed")
def strandcorrect (infile,outfile):
    inputs = P.snip (infile,".bed")
    output = P.snip (outfile,".bed")
    wryner = PARAMS["wryner_script"]
    strand = PARAMS["strand_correct"]
    module = PARAMS["Plink_module"]
    statement = '''module load %(module)s &&  bash update_build.sh %(inputs)s %(strand)s %(output)s ''' 
    P.run(statement)

@transform (strandcorrect, regex(r"(.+)_r2.bed"), r"\1_r3.bed")
def missing1 (infile,outfile):
    inputs = P.snip (infile,".bed")
    output = P.snip (outfile,".bed")
    module = PARAMS["Plink_module"]
    statement = '''module load %(module)s && plink \
                --bfile %(inputs)s \
                --geno 0.2 \
                --make-bed \
                --out %(output)s'''
    P.run(statement)



@transform (missing1, regex(r"(.+)_r3.bed"), r"\1_r4.bed")
def missing2 (infile,outfile):
    inputs = P.snip (infile,".bed")
    output = P.snip (outfile,".bed")
    module = PARAMS["Plink_module"]
    statement = '''module load %(module)s && plink \
                --bfile %(inputs)s \
                --mind 0.2 \
                --make-bed \
                --out %(output)s'''
    P.run(statement)


@transform (missing2, regex(r"(.+)_r4.bed"), r"\1_r5.bed")
def missing3 (infile,outfile):
    inputs = P.snip (infile,".bed")
    output = P.snip (outfile,".bed")
    module = PARAMS["Plink_module"]
    statement = '''module load %(module)s && plink \
                --bfile %(inputs)s \
                --geno 0.02 \
                --make-bed \
                --out %(output)s'''
    P.run(statement)

 
@transform (missing3, regex(r"(.+)_r5.bed"), r"\1_r6.bed")
def maf (infile,outfile):
    inputs = P.snip (infile,".bed")
    output = P.snip (outfile,".bed")
    module = PARAMS["Plink_module"]
    statement = '''module load %(module)s && plink \
                --bfile %(inputs)s \
                --maf 0.05 \
                --make-bed \
                --out %(output)s'''
    P.run(statement) 


@transform (maf, regex(r"(.+)_r6.bed"), r"\1_r7.bed")
def hwe (infile,outfile):
    inputs = P.snip (infile,".bed")
    output = P.snip (outfile,".bed")
    module = PARAMS["Plink_module"]
    statement = '''module load %(module)s && plink \
                --bfile %(inputs)s \
                --hwe 1e-6 \
                --make-bed \
                --out %(output)s'''
    P.run(statement)


@transform (hwe, regex(r"(.+)_r7.bed"),r"\1_r8.bed")

def prunned (infile,outfile):
    module = PARAMS["Plink_module"]
    output = P.snip(outfile,".bed")
    inputs = P.snip(infile,".bed")
    slinding = PARAMS["sliding_windows_size"]
    step_size = PARAMS["step_size"]
    LD = PARAMS["LD"]
    tmps = inputs+".tmp"
    outs = inputs +".bim"
    statement = '''awk '{print $1"\t"$1":"$4":"$5":"$6"\t"$3"\t"$4"\t"$5"\t"$6}' %(outs)s > %(tmps)s &&
                  rm %(outs)s && mv %(tmps)s %(outs)s &&
                  module load %(module)s && 
                  plink --bfile %(inputs)s
                         --indep-pairwise %(slinding)s %(step_size)s %(LD)s
                   --out %(output)s
                   --make-bed'''
    P.run (statement)

@transform (prunned, regex(r"(.+)_r8.bed"),r"\1_r9.bed")

def extract (infile,outfile):
    module = PARAMS["Plink_module"]
    output = P.snip(outfile,".bed")
    inputs = P.snip(infile,".bed")
    snp_list = inputs + ".prune.in"
    outs = inputs +".bim"
    tmps = inputs+".tmp"

    statement = '''module load %(module)s &&
                   plink --bfile %(inputs)s
                   --extract %(snp_list)s
                   --make-bed
                   --out %(output)s'''
    P.run (statement)


@transform (PARAMS["1000G"],formatter(), r"{basename[0]}_r1.bed")

def g1000 (infile,outfile):
    module = PARAMS["Plink_module"]
    output = P.snip(outfile,".bed")
    inputs = P.snip(infile,".bed")

    statement = '''module load %(module)s &&
                   plink --bfile %(inputs)s
                   --geno 0.2 
                   --mind 0.2
                   --make-bed
                   --out %(output)s'''
    P.run (statement)


@transform (g1000,regex(r"(.+)_r1.bed"), r"\1_r2.bed")

def g10001 (infile,outfile):
    out  = P.snip (outfile,".bed")
    inputs =  P.snip (infile,".bed")
    module = PARAMS["Plink_module"]

    statement = '''module load %(module)s &&
                   plink --bfile %(inputs)s
                   --geno 0.02 
                   --mind 0.02
                   --make-bed
                   --out %(out)s'''
    P.run(statement)


@transform (g1000,regex(r"(.+)_r1.bed"),
            add_inputs(prunned),
            r"\1_r3.bed")

def g1000_prune (infiles,outfile):
    bed, prunned = infiles
    out  = P.snip (outfile,".bed")
    inputs =  P.snip (bed,".bed")
    prunned = P.snip(prunned,".bed")
    prunin = prunned+".prune.in"
    module = PARAMS["Plink_module"]

    statement = '''module load %(module)s &&
                   plink --bfile %(inputs)s
                   --extract %(prunin)s
                   --make-bed
                   --out %(out)s'''
    P.run(statement)



@transform (g1000_prune,regex(r"(.+)_r3.bed"),
            add_inputs (extract),
            r"\1_r4.bed")

def g1000_verify (infiles,outfile):
    g1000,geno = infiles
    out  = P.snip (outfile,".bed")
    inputs =  P.snip (geno,".bed")
    ins = inputs+".bim"
    g100 = P.snip (g1000,".bed")
    extratcs = inputs+"extract"
    module = PARAMS["Plink_module"]
    statement = ''' awk '{print$2}' %(ins)s  > %(extratcs)s &&
                   module load %(module)s &&
                   plink --bfile %(g100)s
                   --extract %(extratcs)s
                   --make-bed
                   --out %(out)s'''
    P.run(statement)


@transform (g1000_verify,regex(r"(.+)_r4.bed"),
            add_inputs (extract),
            r"\1_geno_r5.bed")

def g1000_verify1 (infiles,outfile):
    g1000,geno = infiles
    out  = P.snip (outfile,".bed")
    inputs =  P.snip (g1000,".bed")
    ins = inputs+".bim"
    geno = P.snip (geno,".bed")
    extratcs = inputs+"extract"
    module = PARAMS["Plink_module"]
    statement = ''' awk '{print$2}' %(ins)s  > %(extratcs)s &&
                   module load %(module)s &&
                   plink --bfile %(geno)s
                   --extract %(extratcs)s
                   --make-bed
                   --out %(out)s'''
    P.run(statement)


@transform (g1000_verify1, regex(r"(.+)_geno_r5.bed"),
                        add_inputs (g1000_verify),
                             r"merge_genotype_1000G.bed")

def merge (infiles,outfile):
    genotype,g1000 = infiles
    g1000s = P.snip (g1000,".bed")
    gt = P.snip (genotype,".bed")
    module = PARAMS["Plink_module"]
    fam = g1000s+".fam"
    bim = g1000s+".bim"
    out = P.snip (outfile,".bed")
    statement = ''' module load %(module)s &&
                    plink --bfile %(gt)s --bmerge %(g1000)s %(bim)s %(fam)s
                    --make-bed
                    --out %(out)s'''
    P.run(statement)


@transform (merge, regex(r"(.+).bed"),
                        r"merge_genotype_1000G_filtred.bed")

def mergefilter (infile,outfile):
    g1000s = P.snip (infile,".bed")
    module = PARAMS["Plink_module"]
    out = P.snip (outfile,".bed")
    statement = ''' module load %(module)s &&
                    plink --bfile %(g1000s)s  
                    --maf 0.01
                    --hwe 0.000001
                    --mind 0.2
                    --geno 0.2
                    --make-bed
                    --out %(out)s'''
    P.run(statement)


@transform (mergefilter, regex(r"(.+).bed"),r"MDS_merge2")

def pca (infile,outfile):
    g1000s = P.snip (infile,".bed")
    genome = g1000s +".genome"
    module = PARAMS["Plink_module"]
    statement = ''' module load %(module)s &&
                    plink --bfile %(g1000s)s --genome --out %(g1000s)s &&
                    plink --bfile %(g1000s)s 
                    --read-genome %(genome)s 
                    --cluster 
                    --mds-plot 10 
                    --out %(outfile)s'''
    P.run(statement)




@follows(bedile,strandcorrect,missing1,missing2, missing3,maf,hwe,prunned,extract,g1000,g10001,g1000_prune,g1000_verify,g1000_verify1,merge,pca) 

def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))



