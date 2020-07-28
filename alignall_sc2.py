#!/usr/bin/env python3.6

import subprocess, os, string #command line calls, path showing, string things?
import multiprocessing #for parallelizing setuff
import csv #not sure where this is used
import math #not sure where this is used
import sys #program control
import argparse #for easy argument input
import re #for regex searching
import gzip #used for parsing fastq.gz files
import socket #this is used for determining the hostname
import json #this is used for counts (parsing kallisto run_info.json
from builtins import input
from datetime import datetime
from functools import partial
from collections import OrderedDict
import slack

exDir = os.getcwd()


VERSION = 1

npc = multiprocessing.cpu_count()
npc = int(npc*.40)
cur_time = datetime.now().strftime("%y%m%d_%H%M%S")

warningContinue = """WARNING: No arguments were specified, running with the following defaults:
genome=MOUSE
aligner=STAR
fastqDir=fastq
strand=3 (2nd strand)
steps=alignAndCount
indel=no
threads=40pct available
Are you sure you want to continue? (Y/N): """

parser = argparse.ArgumentParser(description="Align fastq reads.")
parser.add_argument("-g", "--genome",    default="MOUSE", type=str,                       help="-g [--genome]       Chose from STAR indexes in /Volumes/Projects/Genomes/ \n Examples: MOUSE, HUMAN (default: MOUSE")
parser.add_argument("-f", "--fastqDir",  default="fastq", type=str,                       help="-f [--fastqDir]     This is the directory in which to find the fastq files to align.\n default: ./fastq")
parser.add_argument("-p", "--steps",     default="all",   type=str, help="-p [--steps] These are the steps in the pipeline to run.\nChoices include all, alignAndQC, align, qcOnly, cluster, and projection.\n(default: alignAndQC)")
parser.add_argument("-t", "--threads",   default=npc,       type=int,                       help="-t [--threads]      Run this many samples concurrently.\n(default: " + str(npc) + ")")
parser.add_argument("--version",         action="store_true",                             help="[--version]         prints the version of the script and then exits") #will be very helpful for when I don't remember which version was used for a project, and what that means for the settings/changes
args = parser.parse_args()

if args.version == True:
        print("alignall main script, version: "+ str(VERSION))
        print("exiting...")
        sys.exit()

genome = args.genome
fastqDir = os.path.abspath(args.fastqDir)
steps = args.steps.lower()
threads = args.threads



runSettings = """########### RUNNING WITH OPTIONS: ###########
######### GENOME:\t\t"""+str(genome)+"""
######### FASTQDIR:\t\t"""+str(fastqDir)+"""
######### STEPS:\t\t"""+str(steps)+"""
######### THREADS:\t\t"""+str(threads)+'\n'
print(runSettings)

if len(sys.argv[1:]) == 0:
        usercontinue = input(warningContinue).upper()
        if usercontinue == 'Y':
                print("Running with defaults...")
        else:
                print("Exiting...")
                sys.exit()

print('Directory of fastq files is ' + fastqDir)
try:
        os.path.isdir(fastqDir)
except:
        print("Fastq directory does not exist. Check your files.")
        sys.exit()
print(steps)
if steps not in ["all", "alignandqc", "align", "cluster", "projection", "all", "qconly"]: # new idea for feature: definition that checks mismatch of 2 betwween entered step and list values, if close, askss user if they wante to continue with suggested step.
        print("The chosen step is not currently supported. Please specify the correct option for -p [--steps].")
        print("Options include all, alignAndQC, align, cluster, projection, qcOnly")
        sys.exit()

if genome == "HUMAN":
        ref = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        refIdx = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
        gtf = "Homo_sapiens.GRCh38.78.gtf"
        refFlat = "humanGRCh38.78_mRNA_refFlat.txt"
        starpath = "STARHUMAN"
        star_miRNApath = "STARsRNA/human/human_miRNA"
        star_tRNApath = "STARsRNA/human/human_tRNA"
        star_yRNApath = "STARsRNA/human/human_yRNA"
        miRNA_genelist = "STARsRNA/human/human_files/matureHUMAN_dna_genes.txt"
        kallistoIDX = "Kallisto_HG38/Homo_sapiens.GRCh38.ncrna.cdna.idx"
        kalmiRNAMature = "miRNA/humanKalisto/matureIndex"
        macsGenome = "hs" #need to skip this in the prefix adding loop
        exonBed = "Homo_sapiens.GRCh38.78.wholeExonicRegions.v3.bed"
        indels = "Mills_and_1000G_gold_standard.indels.hg38.Ens.vcf.gz"
        dbsnpHighCon = "1000G_phase1.snps.high_confidence.hg38.Ens.vcf.gz"
        dbsnp = "Homo_sapiens.fix2.vcf.gz"
        GenomeWindows = "Homo_sapiens.GRCh38.dna_sm.primary_assembly.100bp.bed"

elif genome == "MOUSE":
        ref = "Mus_musculus.GRCm38.dna.primary_assembly.fa"
        refIdx = "Mus_musculus.GRCm38.dna.primary_assembly.fa.fai"
        gtf = "Mus_musculus.GRCm38.78.gtf"
        refFlat = "mouseGRCm38.78_mRNA_refFlat.txt"
        starpath = "STARMOUSE"
        star_miRNApath = "STARsRNA/mouse/mouse_miRNA"
        star_tRNApath = "STARsRNA/mouse/mouse_tRNA"
        star_yRNApath = "STARsRNA/mouse/mouse_yRNA"
        miRNA_genelist = "STARsRNA/mouse/mouseFiles/mouse_miRNA_dna_genes.txt"
        kallistoIDX = "Kallisto_MM38/Mus_musculus.GRCm38.ncrna.cdna.idx"
        kalmiRNAMature = "miRNA/mouseKalisto/matureIndex"
        macsGenome = "mm"
        exonBed = ""
        indels = ""
        dbsnpHighCon = ""
        dbsnp = ""
        GenomeWindows = ""

else:
        print("#### Warning: May have difficulty with scRNA or SNP calling  ####")
        print("#### Genome is not currently fully supported. Good Luck!     ####")

fs_noprefix = [ref, refIdx, gtf, refFlat, starpath, star_miRNApath, star_tRNApath,
                star_yRNApath, miRNA_genelist, kallistoIDX, kalmiRNAMature,
                macsGenome, exonBed, indels, dbsnpHighCon, dbsnp, GenomeWindows]
if socket.gethostname() == "functionalgenomics":
        prefix = "/data/genomes/"
        [ref, refIdx, gtf, refFlat, starpath, star_miRNApath, star_tRNApath, star_yRNApath,
        miRNA_genelist, kallistoIDX, kalmiRNAMature, macsGenome,
        exonBed, indels, dbsnpHighCon, dbsnp, GenomeWindows] = [prefix+x for x in fs_noprefix]
        macsGenome = macsGenome.split("/")[-1]
        unzipInstruct = "zcat"
        whitelistBarcodes = "/data/barcodes/3M-february-2018.txt"

elif socket.gethostname() == "LBCJoshPollack2.ucsf.edu":
        prefix = "/Volumes/Pegasus/Genomes/"
        [ref, refIdx, gtf, refFlat, starpath, star_miRNApath,
        star_tRNApath, star_yRNApath, miRNA_genelist, kallistoIDX, kalmiRNAMature,
        macsGenome, exonBed, indels, dbsnpHighCon, dbsnp, GenomeWindows] = [prefix+x for x in fs_noprefix]
        macsGenome = macsGenome.split("/")[-1]
        unzipInstruct = "gzcat"
else:
        print("not sure what server this is...")
        sys.exit()

alignmentDir = "Alignments/STARsolo"+cur_time
alignmentQCDir = "AlignmentsQC/STARsolo"+cur_time




def main(argv):

    global alignmentDir, alignmentQCDir, threads, exDir

    try: #assign paired/single before checking lane #
            targets = os.listdir(fastqDir)
            targets = [i for i in targets if "fastq.gz" in i] 
            targetsFullPath = [fastqDir+"/"+x for x in targets]
    except OSError:
            print("didn't find a fastq directory...")
            print("cur fastqDir: " + fastqDir)
            sys.exit()
    targetsShort = []
    for f in targets:
        splitf = f.split("_")
        splitf = splitf[:-4]
        targetsFinal = ""
        for l in splitf:
            targetsFinal = targetsFinal+l+"_"
        targetsFinal = targetsFinal[:-1]
        targetsShort.append(targetsFinal)

    targetsShort = set(targetsShort)
    targetsCommaSep = []
    for h in targetsShort:
        fastqMatch = [i for i in targetsFullPath if h in i]
        targetsR1 = [j for j in fastqMatch if "_R1" in j]
        targetsR2 = [k for k in fastqMatch if "_R2" in k]
        targetsR1.sort()
        targetsR2.sort()
        R1s = ""
        R2s = ""
        indTarget = []
        for l in targetsR1:
            R1s = R1s+l+","
        for m in targetsR2:
        	R2s = R2s+m+","
        R1s = R1s[:-1]
        R2s = R2s[:-1]
        indTarget = [R1s,R2s]
        targetsCommaSep.append(indTarget)
    try:
            os.mkdir("AlignAllOut")
    except OSError:
            print("AlignAllOut already exists, starting here.")
            pass
    os.chdir("AlignAllOut") # make the other directory changing a lot easier
    
    runLog = open("runLog.txt", 'w')
    runLog.write(runSettings)

    #make the rest of the important dirs
    for d in ["Alignments", "AlignmentsQC", "SNPcalling", "Peakcalling",
            "SeqQC", "SeqQC/fastqc", "trimdfastqs", "cutadaptLogs"]:
            try:
                    os.mkdir(d)
            except OSError:
                    print("didn't make dir: "+d)
                    pass
    if steps in ["alignandqc", "align", "all"]:
        try:
            os.mkdir(alignmentDir)
            os.mkdir(alignmentQCDir)
        except:
            pass
        align_sc(targetsCommaSep)
        if steps in ["alignandqc", "all"]:
            alignDirs = os.listdir(alignmentDir)
            alignDirs = [i for i in alignDirs if "out" in i]
            alignDirs = [alignmentDir + "/" + i for i in alignDirs]
            print(alignDirs)
            print(alignDirs)
            pool = multiprocessing.Pool(int(threads/2))
            pool.map(cellFilter, alignDirs)
            pool.close()
            pool.join()

    if steps in ["qconly"]:
        try:
            print(os.getcwd())
            if len(os.listdir("Alignments")) > 1:
                user_continue = ""
                print("Multiple alignments have been run,\nwhich alignment dir do you want to perform QC on?\n["+" ".join(os.listdir("Alignments"))+"] ")
                while user_continue not in os.listdir("Alignments"):
                    user_continue = input("\nAwaiting response, ^C to exit: ")
                if user_continue in os.listdir("Alignments"):
                    alignmentDir = "Alignments/"+user_continue
                if user_continue in os.listdir("AlignmentsQC"):
                    alignmentQCDir = "AlignmentsQC/"+user_continue
            alignDirs = os.listdir(alignmentDir)
            alignDirs = [i for i in alignDirs if "out" in i]

            alignDirs = [alignmentDir + "/" + i for i in alignDirs]
            print(alignDirs)
            pool = multiprocessing.Pool(int(threads/2))
            pool.map(cellFilter, alignDirs)
            pool.close()
            pool.join()

        except OSError:
                print("couldn't find the Alignments dir...exiting...")
                sys.exit()





def align_sc(targetsI):
    for i in targetsI:
        if i[0].find("_R1") > 1:
            outfile = getOutfileName(i)
            cwd = os.getcwd()
            if outfile == "nosplitchar":
                print("Sample name formatting is required! (samples must contain _L00), exiting...")
                sys.exit()
            try:
                os.mkdir(alignmentDir+"/"+outfile)
                # if aligner == "kallisto":
                #     os.mkdir("KallistoOuts/"+outfile+"/logs.txt") #don't need this anymore -- just going to need a single alignment dir named for whatever aligner beig used
            except OSError:
                print("couldn't created " + alignmentDir+"/"+outfile)
                print("continuing")
                pass
        #print("running : " + str(target))
        R1fastq = i[0]
        R2fastq = i[1]
        print(outfile)
        print("running: \n" + R2fastq + "\n" + R1fastq + "\nGoing in : " + outfile)
        if R1fastq.endswith('.fastq.gz'):
            rfc = " --readFilesCommand "+unzipInstruct
        else:
            rfc = ""
        #handle paired 
        subprocess.call("STAR --runThreadN "+ str(threads) + " " + "--genomeDir "+starpath+rfc+"  --readFilesIn "+ R2fastq + " " + R1fastq + " \
         --outFileNamePrefix "+alignmentDir+"/"+outfile+"/" + " --outFilterType BySJout   --outFilterMultimapNmax 1   --outFilterMismatchNmax 999   --outFilterMismatchNoverLmax 0.04\
        --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   --alignSJDBoverhangMin 1   --soloType Droplet   --soloUMIlen 12\
        --clip3pAdapterSeq GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --outSAMtype BAM SortedByCoordinate --soloCBwhitelist "+whitelistBarcodes ,shell=True)

def cellFilter(alignDir):
	subprocess.call("Rscript --vanilla " + exDir + "/cellFiltering.R " + alignDir + "/Solo.out",shell=True)


def getOutfileName(name):
    #looking in reverse
    namei = name[0].split(",")[0]
    namei = namei.split("/")[-1]
    if len(namei.split("_L00")) > 1:
        namei = namei.split("_")
        namei = namei[:-4]
        finalName = ""
        for l in namei:
            finalName = finalName+l+"_"
        finalName = finalName[:-1]
        return  finalName + "out"
    else:
         return "nosplitchar"

if __name__ == "__main__":
        main(sys.argv[1:])
