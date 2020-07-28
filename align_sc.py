#!/usr/bin/env python3 

import subprocess, os, sys
import multiprocessing 
import scanpy as sc
import pandas as pd
from matplotlib import rcParams

def getOutfileName(name):
	#looking in reverse
	for i in range(len(name)-2):
		if name[-i-1:-i-4:-1]=="S__":
			return name[:-i-3] + "out"
			
	for i in range(len(name)-1):
		if name[-i-1:-i-3:-1]=="__":
			return name[:-i-2] + "out"
		elif name[-i-1:-i-3:-1]=="S_":
			return name[:-i-2] + "out"
	return "nosplitchar"

def align_sc(targets):
	for i in targets:
		if i.find("_R1") > 1:
			outfile = getOutfileName(i)
			if outfile == "nosplitchar":
				print "Sample name formatting is required! (samples must contain \"__\" or \"_S\"), exiting..."
				sys.exit()
			try:
				os.mkdir(alignmentDir+"/"+outfile)
				# if aligner == "kallisto":
				# 	os.mkdir("KallistoOuts/"+outfile+"/logs.txt") #don't need this anymore -- just going to need a single alignment dir named for whatever aligner beig used
			except OSError:
				pass

	for target in targets:
		if target.endswith('.fastq.gz'):
			rfc = " --readFilesCommand "+unzipInstruct
		else:
			rfc = ""

		#handle paired 
		subprocess.call("STAR --runThreadN 12 --genomeDir /data/genomes/STAR_HUMAN/   --readFilesIn CLH8_fastq/CLH8_S5_L001_R2_001.fastq.gz   CLH8_fastq/CLH8_S5_L001_R1_001.fastq.gz\
		--readFilesCommand zcat --outFileNamePrefix CLH8-starSolo/   --outFilterType BySJout   --outFilterMultimapNmax 1   --outFilterMismatchNmax 999   --outFilterMismatchNoverLmax 0.04\
		--alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   --alignSJDBoverhangMin 1   --soloType Droplet   --soloUMIlen 12\
		--soloCBwhitelist "+whitelistBarcodes)


		STAR --runThreadN 12 --genomeDir /data/genomes/STARHUMAN_272b/ --readFilesIn /data/projects/Jyoti_scRNA/fastq/gexl/061919Fibroid_GEXL_fastqs/061919Fibroid_GEXL_S24_L001_R1_001.fastq.gz /data/projects/Jyoti_scRNA/fastq/gexl/061919Fibroid_GEXL_fastqs/061919Fibroid_GEXL_S24_L001_R2_001.fastq.gz --readFilesCommand zcat --outFileNamePrefix Fibroid_061919/   --outFilterType BySJout   --outFilterMultimapNmax 1   --outFilterMismatchNmax 999   --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   --alignSJDBoverhangMin 1   --soloType Droplet   --soloUMIlen 12 --soloCBwhitelist /data/genomes/barcodes/3M-february-2018.txt



#get sequence QC data
#align sc data
	#do cliipping if necessary
	#get alignment QC and logs
#begin SC analysis 
	

