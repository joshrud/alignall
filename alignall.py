#!/usr/bin/env python3

import subprocess
import os
import string
import sys
import multiprocessing 
import argparse 
import re
import gzip
import socket
import json 
import slack
import PyRanges as pr
from builtins import input
from datetime import datetime
from functools import partial
from collections import OrderedDict


import alignall_snp as snp 
#import alignall_ChIP    #UNTESTED
#import alignall_ATAC    #UNTESTED

VERSION = 14.5 #finally begun new set of file creation, need a new repo and begin separating and testing the files as a single package

warningContinue = """ WARNING: No arguments were specified, running with the following defaults:
genome=MOUSE
aligner=STAR
fastqDir=fastq
strand=1 (forward strand)
steps=alignandcount
indel=no
threads=6
Are you sure you want to continue? (Y/N): """

genome = ''
fastqDir = ''
strand = ''
steps = ''
machine = socket.gethostname()
cur_time = datetime.now().strftime("%y%m%d_%H%M%S") #a string representing the time when program executes
EST_COUNTS = 3 #the column to extract counts from in abundance.tsv in kallisto alignments (4 would be TPM, may want this option later...)
GATKbase="~/bin/gatk-4.0.11.0/gatk " 
adapter = ""
whitelistBarcodes = "/data/barcodes/3M-february-2018.txt" 

parser = argparse.ArgumentParser(description="Align fastq reads.")
parser.add_argument("-g", "--genome",                   default="MOUSE",                type=str,                       help="-g [--genome]       Chose from STAR indexes in /Volumes/Projects/Genomes/ \n Examples: MOUSE, HUMAN (default: MOUSE")
parser.add_argument("-a", "--aligner",                  default="STAR",                 type=str,                       help="-a [--aligner]      Chose from kallisto (transcript level counts) or STAR (gene level counts) (default: STAR)")
parser.add_argument("-ad", "--adapter",                 default="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",     type=str,      help="-ad [--adapter]      adapter to trim, default is GATCGGAAGAGCACACGTCTGAACTCCAGTCAC, another common one is Nextera transposase: CTGTCTCTTATACACATCTCCGAGCCCACGAGAC")
parser.add_argument("-f", "--fastqDir",                 default="fastq",                type=str,                       help="-f [--fastqDir]     This is the directory in which to find the fastq files to align.\n default: ./fastq")
parser.add_argument("-s", "--strand",                   default=1,                      type=int,                       help="-s [--strand]       This is the strand of the RNA-seq library prep. \nChoices are 0(unstranded: Ovation/Smartseq), 1 (first strand: Lexogen FWD, and Nugen Universal) and 2 (second strand: TruSeq). (default: 3)")
parser.add_argument("-p", "--steps",                    default="alignandcount",        type=str,                       help="-p [--steps]        These are the steps in the pipeline to run.\nChoices include align, count, picard, alignandcount, callsnps, chip, atac, callpeaks, mirna, smallrna, scrna, scatac.(default: alignandcount)")
parser.add_argument("-t", "--threads",                  default=6,                      type=int,                       help="-t [--threads]      Run this many samples concurrently.\nThis sets the number of GATK/picard samples to run in parallel.\nRecommend setting this to 1/3 total thread count on a given machine.\nThese java commands can often occupy 2-4 threads each.(default: 6)")
parser.add_argument("--version",                                action="store_true",                           	 		help="[--version]         prints the version of the script and then exits") #will be very helpful for when I don't remember which version was used for a project, and what that means for the settings/changes
args = parser.parse_args()

if args.version == True:
	print("alignall main script, version: "+ str(VERSION))
	print("exiting...")
	sys.exit()

genome = args.genome
aligner = args.aligner
randbar = args.random_barcodes
adapterSeq = args.adapter
fastqDir = os.path.abspath(args.fastqDir)
strand = args.strand
steps = args.steps.lower()
threads = args.threads

runSettings = "########### RAN AT: "+cur_time+"\n"
"########### RUNNING WITH OPTIONS: ###########\n"
"######### GENOME:\t\t"+str(genome)+"\n"
"######### ALIGNER:\t\t"+str(aligner)+"\n"
"######### RANDOM_BARCODES:\t\t" +str(randbar)+" (unused for now)"+"\n"
"######### FASTQDIR:\t\t"+str(fastqDir)+"\n"
"######### STRAND:\t\t"+str(strand)+"\n"
"######### STEPS:\t\t"+str(steps)+"\n"
"######### INDEL:\t\t"+str(indel)+"\n"
"######### THREADS:\t\t"+str(threads)+'\n'
print(runSettings)

if len(sys.argv[1:]) == 0:
	usercontinue = input(warningContinue).upper()
	if usercontinue == 'Y':
		print("Running with defaults...")
	else:
		print("Exiting...")
		sys.exit()

######################################################################
## Assigning Directories for Reference Files and Indices            ##
##                                                                  ##
######################################################################

print('Directory of fastq files is ' + fastqDir)
try:
	os.path.isdir(fastqDir)
except:
	print("Fastq directory does not exist. Check your files.")
	sys.exit()

#assign strand constant info
strands = ["0","1","2"]
if str(strand) in strands:
	if float(strand) == float(0):
		htSeqStrand = "no"
		print("Strand is defined as 0, which implies nonstranded preps (e.g.Nugen Ovation or Illumina Truseq Unstranded)")
	if float(strand) == float(1):
		htSeqStrand = "yes"
		print("Strand is defined as 1, which implies first strand synthysis prep (probably lexogen forward")
	if float(strand) == float(2):
		htSeqStrand = "reverse"
		print("Strand is defined as 2, which implies second strand synthesis (e.g. Illumina Truseq Stranded)")
else:
	print("please choose from 0, 1, or 2.  Exiting...")
	sys.exit()

#make sure we have chosen a valid step
if steps not in ["align", "count", "alignandcount", "callsnps", "chip", "atac", "scrna", "scatac", "smallrna", "mirna", "mirnacount","picard"]: # new idea for feature: definition that checks mismatch of 2 betwween entered step and list values, if close, askss user if they wante to continue with suggested step.
	print("The chosen step is not currently supported. Please specify the correct option for -p [--steps].")
	print("Options include align, count, alignandcount, callsnps, chip, atac, scrna, smallrna, mirna, mirnacount, picard")
	sys.exit()
if steps == "scrna" and machine != "functionalgenomics":
	print("single cell rna alignment must be run on the functional genomics server! Exiting..")
	sys.exit()

#avoid running into problems for miRNA study
if steps=="mirna" and aligner != "STAR":
	print("This script only supports STAR for performing alignments on micro rna, \n Switching aligner to STAR...")
	aligner = "STAR"
if steps=="smallrna" and aligner != "STAR":
	print("This script only supports STAR for performing alignments on extracellular rna, \n Switching aligner to STAR...")
	aligner = "STAR"

#set all directory constants
if genome == "HUMAN":
	ref = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
	refIdx = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
	gtf = "Homo_sapiens.GRCh38.95.gtf"
	refFlat = "Homo_sapiens.GRCh38.95.refFlat.txt"
	riboIntervals = "Homo_sapiens.GRCh38.95_rRNAintervals.txt"
	starpath = "STARHUMAN"
	star_miRNApath = "STARsRNA/human/human_miRNA"
	star_tRNApath = "STARsRNA/human/human_tRNA"
	star_yRNApath = "STARsRNA/human/human_yRNA"
	miRNA_genelist = "STARsRNA/human/human_files/miRNA_files/matureHUMAN_dna_genes.txt"
	kallistoIDX = "Kallisto_HG38/Homo_sapiens.GRCh38.ncrna.cdna.idx"
	kalmiRNAMature = "miRNA/humanKalisto/matureIndex"
	macsGenome = "hs" #need to skip this in the prefix adding loop
	exonBed = "Homo_sapiens.GRCh38.95.wholeExonicRegions.v2.list"
	indels = "dbsnp_hsa/Mills_and_1000G_gold_standard.indels.hg38.Ens.vcf.gz" #not actually there
	dbsnpHighCon = "dbsnp_hsa/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
	dbsnp = "dbsnp_hsa/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
	GenomeWindows = "Homo_sapiens.GRCh38.dna_sm.primary_assembly.100bp.bed"

elif genome == "HUMANECOLI":
	ref = "Homo_sapiens.GRCh38.dna.primary_assembly_withEcoli.fa"
	refIdx = "Homo_sapiens.GRCh38.dna.primary_assembly_withEcoli.fa.fai"
	gtf = "Homo_sapiens.GRCh38.95_Ecoli.gtf"
	refFlat = "humanGRCh38.78_mRNA_refFlat.txt"
	riboIntervals = ""
	starpath = "STARHUMAN_ECOLI"
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

elif genome == "HUMANHG19":
	ref = "Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa"
	refIdx = ""
	gtf = "Homo_sapiens.GRCh37.75.gtf"
	refFlat = ""
	riboIntervals = ""
	starpath = "STARHUMANHG19"
	star_miRNApath = ""
	star_tRNApath = ""
	star_yRNApath = ""
	miRNA_genelist = ""
	kallistoIDX = ""
	kalmiRNAMature = ""
	macsGenome = ""
	exonBed = ""
	indels = ""
	dbsnpHighCon = ""
	dbsnp = ""
	GenomeWindows = ""

elif genome == "MOUSE":
	ref = "Mus_musculus.GRCm38.dna.primary_assembly.fa"
	refIdx = "Mus_musculus.GRCm38.dna.primary_assembly.fa.fai"
	gtf = "Mus_musculus.GRCm38.96.gtf"
	refFlat = "Mus_musculus.GRCm38.96.refflat.txt"
	riboIntervals = "" #need to generate a v96 riboIntervals
	starpath = "STARMOUSE"
	star_miRNApath = "STARsRNA/mouse/mouse_miRNA"
	star_tRNApath = "STARsRNA/mouse/mouse_tRNA"
	star_yRNApath = "STARsRNA/mouse/mouse_yRNA"
	miRNA_genelist = "STARsRNA/mouse/mouseFiles/mouse_miRNA_dna_genes.txt"
	kallistoIDX = "Kallisto_MM38/Mus_musculus.GRCm38.ncrna.cdna.idx"
	kalmiRNAMature = "miRNA/mouseKalisto/matureIndex"
	macsGenome = "mm"
	exonBed = "Mus_musculus.GRCm38.96wholeExonic.bed"
	indels = "GRCm38_snps/mgp.v3.indels.rsIDdbSNPv137.vcf.gz"
	dbsnpHighCon = "GRCm38_snps/mgp.v3.snps.rsIDdbSNPv137.vcf.gz"
	dbsnp = "GRCm38_snps/mgp.v3.snps.rsIDdbSNPv137.vcf.gz"
	GenomeWindows = ""

elif genome == "RAT":
	ref = "Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
	refIdx = ""
	gtf = "Rattus_norvegicus.Rnor_6.0.94.gtf"
	refFlat = "Rattus_norvegicus.Rnor_6.0.94.refFlat.txt"
	riboIntervals = ""
	starpath = "STARRAT"
	star_miRNApath = "STARsRNA/rat_miRNA"
	star_tRNApath = "STARsRNA/rat_tRNA"
	star_yRNApath =  "STARsRNA/rat_yRNA"
	miRNA_genelist = "STARsRNA/rat/ratFiles/ratGenes.txt"
	kallistoIDX = ""
	kalmiRNAMature = ""
	macsGenome = ""
	exonBed = ""
	indels = ""
	dbsnpHighCon = ""
	dbsnp = ""
	GenomeWindows = ""

elif genome == "ZEBRAFISH":
	ref = "Danio_rerio.GRCz10.dna.toplevel.fa"
	refIdx = ""
	gtf = "Danio_rerio.GRCz10.87.gtf"
	refFlat = "Danio_rerio.GRCz10.87.refFlat.txt"
	riboIntervals = ""
	starpath = "STARZEBRAFISH"
	star_miRNApath = ""
	star_tRNApath = ""
	star_yRNApath = ""
	miRNA_genelist = ""
	kallistoIDX = ""
	kalmiRNAMature = ""
	macsGenome = ""
	exonBed = ""
	indels = ""
	dbsnpHighCon = ""
	dbsnp = ""
	GenomeWindows = ""

elif genome == "COW":
	ref = "Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"
	refIdx = ""
	gtf = "Bos_taurus.ARS-UCD1.2.95.gtf"
	refFlat = "Bos_taurus.ARS-UCD1.2.95.refFlat_v3.txt"
	riboIntervals = ""
	starpath = "STARCOW"
	star_miRNApath = "STARsRNA/cow/cow_miRNA"
	star_tRNApath = ""
	star_yRNApath = ""
	miRNA_genelist = ""
	kallistoIDX = ""
	kalmiRNAMature = ""
	macsGenome = ""
	exonBed = ""
	indels = ""
	dbsnpHighCon = ""
	dbsnp = ""
	GenomeWindows = ""

elif genome == "CHICKEN":
	ref = ""
	refIdx = ""
	gtf = "Gallus_gallus.Gallus_gallus-5.0.90.gtf"
	refFlat = "Gallus_gallus.Gallus_gallus-5.0.90_refflat_v2.txt"
	riboIntervals = ""
	starpath = "STARCHICKEN"
	star_miRNApath = ""
	star_tRNApath = ""
	star_yRNApath = ""
	miRNA_genelist = ""
	kallistoIDX = ""
	kalmiRNAMature = "miRNA/chickenKalistok15/matureCHICKENk15"
	macsGenome = "gg"
	exonBed = ""
	indels = ""
	dbsnpHighCon = ""
	dbsnp = ""
	GenomeWindows = ""

elif genome == "DUCKEN":
	ref = ""
	refIdx = ""
	gtf = "duck_chicken_combined.gtf"
	refFlat = "duckchicken_refflatTest2.txt"
	riboIntervals = ""
	starpath = "STARDUCKEN"
	star_miRNApath = ""
	star_tRNApath = ""
	star_yRNApath = ""
	miRNA_genelist = ""
	kallistoIDX = ""
	kalmiRNAMature = ""
	macsGenome = ""
	exonBed = ""
	indels = ""
	dbsnpHighCon = ""
	dbsnp = ""
	GenomeWindows = ""

elif genome == "S_AUREUS":
	ref = ""
	refIdx = ""
	user_continue = input("S_AUREUS has a few genomes, which would you like to use? (TCH1516 or RIBO): ").upper()
	if user_continue == "TCH1516":
		starpath = "STARBACTERIA/S_AUREUS_tch1516"
		gtf = "Staphylococcus_aureus_subsp_aureus_usa300_tch1516.ASM1708v1.43.gtf"
		refFlat = "" #this needs to be added
		print("starpath is now tch1516")
	elif user_continue == "RIBO":
		print("using STAR_RIBO...")
		starpath = "STARBACTERIA/S_AUREUS_ribo"
		gtf = "Staphylococcus_aureus_GCF_001027105/GCF_001027105.1_ASM102710v1_genomic.v3.gtf"
		refflat = "Staphylococcus_aureus_GCF_001027105/GCF_001027105.1_ASM102710v1_genomic.refflat.txt"
	else:
		print("no subspecies specified, exiting...")
		sys.exit()
	riboIntervals = ""
	star_miRNApath = ""
	star_tRNApath = ""
	star_yRNApath = ""
	miRNA_genelist = ""
	kallistoIDX = ""
	kalmiRNAMature = ""
	macsGenome = ""
	exonBed = ""
	indels = ""
	dbsnpHighCon = ""
	dbsnp = ""
	GenomeWindows = ""

# if we don't have the genome already in here, another option is to just 
# create a text file with formatting that looks exactly like the reference 
# lists above, then check the file to fill in variables
# HINT: reference files with NO PREFIX are required, prefix for system 
# 	 is added after reference file is read
elif genome == "OTHER":
	reference_input = input("You have chosen a genome not in [HUMAN, MOUSE, "+
		"RAT, ZEBRAFISH, COW, HORSE, OR CHICKEN]\n Please specify a "+
		"file to use as reference: ")

	if os.path.exists(reference_input):
		lines = []
		with open(reference_input, 'rU') as refs:
			for line in refs:
				lines.append(line.rstrip("\n"))
		if len(lines) == 18:
			ref = lines[0].split("= ")[1].strip("\"")
			refIdx = lines[1].split("= ")[1].strip("\"")
			gtf = lines[2].split("= ")[1].strip("\"")
			refFlat = lines[3].split("= ")[1].strip("\"")
			riboIntervals = lines[4].split("= ")[1].strip("\"")
			starpath = lines[5].split("= ")[1].strip("\"")
			star_miRNApath = lines[6].split("= ")[1].strip("\"")
			star_tRNApath = lines[7].split("= ")[1].strip("\"")
			star_yRNApath = lines[8].split("= ")[1].strip("\"")
			miRNA_genelist = lines[9].split("= ")[1].strip("\"")
			kallistoIDX = lines[10].split("= ")[1].strip("\"")
			kalmiRNAMature = lines[11].split("= ")[1].strip("\"")
			macsGenome = lines[12].split("= ")[1].strip("\"")
			exonBed = lines[13].split("= ")[1].strip("\"")
			indels = lines[14].split("= ")[1].strip("\"")
			dbsnpHighCon = lines[15].split("= ")[1].strip("\"")
			dbsnp = lines[16].split("= ")[1].strip("\"")
			GenomeWindows = lines[17].split("= ")[1].strip("\"")
			for v in [ref, refIdx, gtf, refFlat, riboIntervals, starpath, star_miRNApath, 
		star_tRNApath, star_yRNApath, miRNA_genelist, kallistoIDX, kalmiRNAMature,
		macsGenome, exonBed, indels, dbsnpHighCon, dbsnp, GenomeWindows]:
				if v != "":
					print("Using: " + v)
		else:
			print("File formatting required, exiting...")
			sys.exit()

	else:
		print("Path to file not found, exiting...")
		sys.exit()
	

else:
	print("#### Warning: May have difficulty with scRNA or SNP calling  ####")
	print("#### Genome is not currently fully supported. Good Luck!     ####")

fs_noprefix = [ref, refIdx, gtf, refFlat, riboIntervals, starpath, star_miRNApath, 
		star_tRNApath, star_yRNApath, miRNA_genelist, kallistoIDX, kalmiRNAMature,
		macsGenome, exonBed, indels, dbsnpHighCon, dbsnp, GenomeWindows]
if machine == "functionalgenomics":
	prefix = "/data/genomes/"
	[ref, refIdx, gtf, refFlat, riboIntervals, starpath, star_miRNApath, star_tRNApath, 
	star_yRNApath, miRNA_genelist, kallistoIDX, kalmiRNAMature, macsGenome,
	exonBed, indels, dbsnpHighCon, dbsnp, GenomeWindows] = [prefix+x for x in fs_noprefix]
	macsGenome = macsGenome.split("/")[-1]
	unzipInstruct = "zcat"
elif machine == "LBCJoshPollack2.ucsf.edu":
	prefix = "/Volumes/Pegasus/Genomes/"
	[ref, refIdx, gtf, refFlat, riboIntervals, starpath, star_miRNApath,
	star_tRNApath, star_yRNApath, miRNA_genelist, kallistoIDX, kalmiRNAMature,
	macsGenome, exonBed, indels, dbsnpHighCon, dbsnp, GenomeWindows] = [prefix+x for x in fs_noprefix]
	macsGenome = macsGenome.split("/")[-1]
	unzipInstruct = "gzcat"
else:
	print("not sure what server this is...")
	sys.exit()

fs_withprefix = [ref, refIdx, gtf, refFlat, riboIntervals, starpath, 
	star_miRNApath, star_tRNApath, star_yRNApath, miRNA_genelist, kallistoIDX,
	kalmiRNAMature, macsGenome, exonBed, indels, dbsnpHighCon, dbsnp, 
	GenomeWindows]

alignmentDir = "Alignments/"+aligner+cur_time
alignmentQCDir = "AlignmentsQC/"+aligner+cur_time
trimdfolder = "trimdfastqs"

if machine == "functionalgenomics":
	picard_jarfile = "/usr/local/bin/picard.jar"
elif machine == "LBCJoshPollack2.ucsf.edu":
	picard_jarfile = "/Volumes/Pegasus/Projects/dependencies/picard.jar"


#last check, make sure all files for reference species are there 
isthere = list(map(lambda x:[x, os.path.exists(x)], fs_withprefix))
filesnotinplace = False
for i in isthere:
	if not i[1]:
		print(i[0]+" does not exist")
		filesnotinplace = True
if filesnotinplace:
	print("files aren't in place for running: "+steps+" on: "+genome)
	usercontinue = input("Do you still want to continue?(Y/n): ")
	if usercontinue.lower() != "y":
		sys.exit()




######################################################################
## Running the main programs                                        ##
##                                                                  ##
######################################################################

def main(argv):

	global alignmentDir, alignmentQCDir, threads, adapter
	alignallout = False

	print("working in: " + os.path.abspath("."))

	try: #assign paired/single before checking lane #
		targets = os.listdir(fastqDir)
		targets = [i for i in targets if "fastq.gz" in i]
		targetsFullPath = [fastqDir+"/"+x for x in targets]

	except OSError:
		print("didn't find a fastq directory...")
		print("cur fastqDir: " + fastqDir)
		sys.exit()

	print(targets)
	runType = runMode(targets)
	print(runType)

	#check if the samples are split into different lanes
	multipleLanes = False
	alreadyCombined = False
	lanedict = {} #can't just check if there is already a dict entry, becuase the R1/R2 gets cut off, won't work if paired end library
	for f in targets:
		if "L000" in f:
			alreadyCombined = True
			break
		splitf = f.split("_L00")
		if splitf[0] in lanedict.keys():
			lanedict[splitf[0]].append(splitf[1].split("_")[0])
		else: lanedict[splitf[0]] = [splitf[1].split("_")[0]]
	if alreadyCombined:
		print("lanes have already been combined! resetting targets..")
		targets = [x for x in os.listdir(fastqDir) if "L000" in x] #resets targets
	else:
		for k,v in lanedict.items():
			if len(set(v)) > 1:
				multipleLanes = True
				lanes = set(v)

		if multipleLanes:
			usercontinue = input("Multiple lanes of samples have been detected, would you like to combine?: ").lower()
			try:
				print("est. memory use: " + str(int(sum([os.stat(x).st_size for x in targetsFullPath])/1e6)) + "MB")
			except OSError:
				print("can't find any files in fastqDir...")
				sys.exit()
			if usercontinue in ['y', "yes"]:
				libraries = [x.split("_L00")[0] for x in targets]
				print(targets)
				libraries = set(libraries)
				pool = multiprocessing.Pool(int(threads)) # run this many threads concurrently
				func = partial(combineLanes, runType, targets)
				pool.map(func, libraries)
				pool.close()
				pool.join()	
			try:
				os.mkdir(fastqDir+"_original")
			except:
				print("Fastq directory does not exist. Check your files.")
			for i in targets:
				os.rename(fastqDir + "/" + i, fastqDir + "_original/" + i)

			targets = [x for x in os.listdir(fastqDir) if "L000" in x] #resets targets

	#cases= steps (kinda)
	# case1 (most common):  need to align and count bulk rna
		#subcase1: miRNA or smallRNA
		#subcase2: transcriptome instead of genome
	# case2: count ONLY previously aligned bulk rna
		#subcase1: miRNA or smallRNA
		#subcase2: transcriptome
	# case3: alignment and peak calling
		#subcase1: ChIP
		#subcase2: ATAC
	#case4: alignment and snp calling

	#first thing to do is check if an AlignAllOut basedir already exists...
	try:
		os.mkdir("AlignAllOut")
	except OSError:
		print("AlignAllOut already exists, starting here.")
		alignallout = True
		pass
	os.chdir("AlignAllOut") # make the other directory changing a lot easier `

	#really basic logfile for the run (what was run)
	runLog = open("runLog.txt", 'a') #open for appending, new log for ea. run
	runLog.write(runSettings)

	#make the rest of the important dirs
	if not alignallout:
		for d in ["Alignments", "AlignmentsQC", "SNPcalling", "Peakcalling",
			"SeqQC", "SeqQC/fastqc"]:
			try:
				os.mkdir(d)
			except OSError:
				print("didn't make dir: "+d)
				pass

	#now that AlignAllOut is known to exist, let's check Alignments for most recent alignmentDir
	R1_targets = list(filter(lambda x: "_R1" in x, targets))
	alignmentDir_new = alignmentAlreadyDone(R1_targets)
	if alignmentDir_new is not None: 
		print("found a completed alignment, switching alignment dir to: " + alignmentDir_new)
		alignmentDir = "Alignments/"+alignmentDir_new
		alignmentQCDir = alignmentDir.replace("Alignments", "AlignmentsQC")
	else:
		if steps in [ "callsnps", "count", "picard"]:
			print("A full alignment hasn't been completed yet. Run align before ("+steps+")")
			print("Exiting...")
			sys.exit()

	#maybe should just change dir to AlignAllOut right here...
	if steps in ["count", "picard", "mirnacount"]:
		try:
			if len(os.listdir("Alignments")) > 1:
				user_continue = ""
				print("Multiple alignments have been run,\nwhich alignment dir do you want to perform picard stats or counts for?\n["+" ".join(os.listdir("Alignments"))+"] ")
				while user_continue not in os.listdir("Alignments"):
					user_continue = input("\nAwaiting response, ^C to exit: ")
				if user_continue in os.listdir("Alignments"):
					alignmentDir = "Alignments/"+user_continue
				if user_continue in os.listdir("AlignmentsQC"):
					alignmentQCDir = "AlignmentsQC/"+user_continue
			elif len(os.listdir("Alignments")) == 0:
				print("No Alignment has been performed, exiting...")
				sys.exit()
			else: #a single alignment has been done (ideal case)
				print("Alignment has already been performed, using existing alignment dir... ")
				alignmentDir = "Alignments/"+os.listdir("Alignments")[0]
				alignmentQCDir = "AlignmentsQC/"+os.listdir("AlignmentsQC")[0]

			if steps == "picard":
				runPicardRNAMetrics()
				subprocess.call("multiqc -f -m star -m picard -d "+alignmentDir+" -o "+alignmentQCDir, shell=True)
			elif steps == "count":
				counts()
			elif steps == "mirnacount":
				counts_miRNA()

		except Exception as e:
			print(e)
			# print("couldn't find the Alignments dir...exiting...")
			sys.exit()


	if steps not in ["count", "picard", "mirnacount"]: ##["align", "count", "alignandcount", "callsnps", "chip", "atac", "scrna", "smallrna", "mirna"]

		if alignmentDir_new == None:
			try:
				os.mkdir(alignmentDir)
				os.mkdir(alignmentQCDir)
			except OSError:
				print("couldn't make a new alignment dir... exiting")
				sys.exit()
			
		if aligner == "kallisto" and steps in ["mirna", "smallrna"]: #because we need some extra directories
			try:
				os.mkdir(alignmentDir+"/Logfiles")
			except OSError:
				usercontinue = input("unable to make kallisto small/micro rna folders...\n do you want to continue? (Y/N): ").upper()
				if usercontinue != 'Y':
					print("Exiting...")
					sys.exit()

		if steps in ["align", "alignandcount", "smallrna", "mirna"]: #all the basic stuff
			
			doCutadapt = True
			if steps in ["smallrna", "mirna"]:
				try:
					os.mkdir("cutadapt_logs")
					os.mkdir("trimdfastqs")
				except OSError:
					print("not making cutadapt log folders or trimdfastqs folder...")

				trimdFilesInFastqDir = [f for f in targets if "trimd" in f] #sets a standard -> if we already trimmed, then files should be in fastqDir and program will detect if they're already trimmed
				if len(trimdFilesInFastqDir) == len(targets): #trimmed files are all present -> no need for cutadapt
					doCutadapt = False
					print("copying fastqDir to trimdfolder")
					subprocess.call("cp "+fastqDir+"/* "+trimdfolder, shell=True) #we'll be looking for cutadapt files in trimd folder later
					subprocess.call("touch cutadapt_logs/cutadapt_run_prior_NO_LOGS.txt", shell=True)

				if doCutadapt:
					cutadaptAll(targets) #makes trimdfolder fastqs
					targets = os.listdir(trimdfolder)
					print("done clipping...")

			print("In dir: "+os.getcwd())
			print("aligning...")
			align(fastqDir, targets, runType,threads) #should work not matter what biotype it is
			print("running fastqc...")
			subprocess.call("fastqc -t " + str(threads) + " -o SeqQC/fastqc " + fastqDir + "/*", shell=True)
			subprocess.call("multiqc -f -m fastqc -d SeqQC/fastqc -o SeqQC/", shell=True)

			if aligner == "STAR":
				if steps in ["smallrna", "mirna"]:
					subprocess.call("multiqc -f -m star -m cutadapt -d "+alignmentDir+"/ -o AlignmentsQC/", shell=True)
				else:
					print("running picard RNAseq Metrics")
					if refFlat != None:
						runPicardRNAMetrics()
						subprocess.call("multiqc -f -m star -m picard -d "+alignmentDir+" -o "+alignmentQCDir, shell=True)
					else:
						print("no refflat file, continuing with only star QC")
						subprocess.call("multiqc -f -m star -d "+alignmentDir+" -o "+alignmentQCDir, shell=True)
			elif aligner == "kallisto":
				if steps == "miRNA" or steps =="smallRNA":
					subprocess.call("multiqc -f -m kallisto -m cutadapt -d "+alignmentDir+"/Logfiles -d cutadapt_logs -o "+alignmentQCDir+"/", shell=True)
				else:
					subprocess.call("multiqc -f -m kallisto -d "+alignmentDir+"/Logfiles -o "+alignmentQCDir+"/", shell=True)

				#long term idea: do STAR alignment and kallisto alignment, then validate against eaech other, might also find novel seqs (basically code out mirbase, but make it actually work)

			if steps == "alignandcount":
				print("Gathering counts...")
				counts()
				generateMD5s()
			elif steps == "mirna":
				print("Gathering miRNA counts...")
				counts_miRNA() #this needs to be fixed I think, won't quite owrk just yet.
			elif steps == "smallrna":
				print("not yet configured, counting needs to be done manually")


		elif steps == "chip": #this needs a lot of work
			peakcall_dir = "peakcalling_"+cur_time
			chipgraphs_dir = "chipgraphs_"+cur_time
			try:
				os.mkdir(peakcall_dir)
				os.mkdir(chipgraphs_dir)
			except OSError:
				#need to add a check here to see if a project has already been run (or halfway through) so we can continue with it
				pass

			align(fastqDir, targets, runType,threads)

			subprocess.call("fastqc -t " + str(threads) + " -o SeqQC/fastqc " + fastqDir + "/*", shell=True)
			subprocess.call("multiqc -f -m fastqc -d SeqQC/fastqc -o SeqQC/", shell=True)

			subprocess.call("multiqc -f -m star -m picard -d "+alignmentDir+" -o "+alignmentQCDir, shell=True)


			#this should basically be the same as below, but the references to the files' methods should be slightly different....

		elif steps == "atac":
			#need to check (or ask for) an alignment dir to use first
			#also going to have to change a lot of the below directory paths - want to use alignmentDir, not Alignments

			numSamples = len([x for x in os.listdir(alignmentDir) if x.endswith("out")])
			if numSamples < int(threads):
				print("too many threads, setting it to maximum samples ("+str(numSamples)+")")
				threads = numSamples

			#need to set some names for the conditions in the experiment
			splitNames = [ n.split("_") for n in os.listdir(alignmentDir) ]#assuming this is AlignAllOut
			lengths = [ len(x) for x in splitNames ]
			if len(set(lengths)) != 1:
				#there is >1 type of sample splitting, so at least 1 sample doesn't fit the schema
				print("sample name formatting is required...")
				sys.exit() #maybe should return somerhing instead, since this will eventually be called by alignall
			splitlength = list(set(lengths))

			filecount=0
			alignmentSampleDirs = [x for x in os.listdir(alignmentDir) if os.path.isdir(os.path.join(alignmentDir, x))]
			for d in alignmentSampleDirs:
				if "Aligned.sortedByCoord.duplicateRemoved.out.bam" in os.listdir(alignmentDir+"/"+d):
					filecount+=1
			print("filecount: "+str(filecount))
			print("directories: "+str(len(alignmentSampleDirs)))
			if filecount == len(alignmentSampleDirs):
				print("files already exist")
			else:
				print("running runRmDups")
			snp.runRmDups(alignmentDir, threads, picard_jarfile)
			#print("finished")
			#print("Peak calling")
			#runCallPeaks_atac()
			#print("Running HTseq")
			#runChipHTseq_atac() #rename all of these in the file
			#print("Gathering htseq counts")
			#ChipCounts_atac()
			#print("Making ChIP-seq graphs")
			#runMakeChipGraphs_atac()

		elif steps == "scrna":
			print("Removing duplicates and adding read groups")
			snp.runRmDups(alignmentDir, threads, picard_jarfile)
			print("running HTseq")
			runHTseq()
			print("gathering htseq counts")
			scRnaCounts()

		elif steps == "callsnps": #working on this now...
			#reset alignmentDir if we have already done it 
			#each function needs to return a True/false and then success needs to equal all of them &'d together -only then can we delete all of the files

			print("Running SNP calling")
			# print("Removing duplicates")
			
			print("haven't figured out how to split files or do condition/sample specification, so just enter \"condition\" or \"sample\" ")
			user_continue = input("(condition/sample): ")
			if user_continue == "condition":
				readGroupOption = "condition"
			elif user_continue == "sample":
				print("ok, splitting files by sample...")
				readGroupOption = "sample"
			else:
				print(user_continue+" not recognized, exiting...")
				sys.exit()

			# if snp.checkReplicates(): #true = there are replicates; false = no replicates (checked by taking last split section using underscore)
			# 	readGroupOption = "condition"
			# else:
			# 	readGroupOption = "sample"

			runRmDupsCOMPLETE = False
			runAssignGroupsCOMPLETE = False
			runPreGATKCOMPLETE = False
			runMergeAndHaplotypeCallerCOMPLETE = False

			runRmDupsCOMPLETE = snp.runRmDups(alignmentDir,	 picard_jarfile, threads)
			print("adding read groups")
			runAssignGroupsCOMPLETE = snp.runAssignGroups(alignmentDir, picard_jarfile, fastqDir, threads, readGroupOption)
			print("Running GATK preprocessing (splitting bams)") #according to https://pubmed.ncbi.nlm.nih.gov/24974202/ effect is minimal when high coverage
			runPreGATKCOMPLETE = snp.runPreGATK(alignmentDir, threads, ref, indels, dbsnp, dbsnpHighCon, indel, GATKbase)
			print("Calling SNPs")
			runMergeAndHaplotypeCallerCOMPLETE = snp.runMergeAndHaplotypeCaller(exonBed, threads, ref, dbsnp, alignmentDir, GATKbase)

			if runRmDupsCOMPLETE and runAssignGroupsCOMPLETE and runPreGATKCOMPLETE and runMergeAndHaplotypeCallerCOMPLETE: 
				print("all processes succeeded. finishing up.")
			else:
				print("one or more processes failed...")
				print(["runRmDupsCOMPLETE", "runAssignGroupsCOMPLETE", "runPreGATKCOMPLETE", "runMergeAndHaplotypeCallerCOMPLETE"])
				print([runRmDupsCOMPLETE, runAssignGroupsCOMPLETE, runPreGATKCOMPLETE, runMergeAndHaplotypeCallerCOMPLETE])
				sys.exit()


	print("Cleaning up")
	print("Due to size constraints some files are removed")
	if steps in ["callsnps", "chip", "atac"]:
		if success:
			subprocess.call("rm "+alignmentDir+"/*/*bed", shell = True)
			subprocess.call("rm "+alignmentDir+"/*/*bed", shell = True)
			subprocess.call("rm "+alignmentDir+"/*/*bg", shell = True)
			# subprocess.call("rm "+alignmentDir+"/*/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam", shell = True)
			# subprocess.call("rm "+alignmentDir+"/*/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.bam", shell = True)
			# subprocess.call("rm "+alignmentDir+"/*/Aligned.sortedByCoord.duplicateRemoved.out.bam", shell = True)
		else:
			print("Didn't successfully finish the pipeline on step: ", steps)

	if steps in ["alignandcount", "align", "callsnps", "scrna", "scatac", "atac", "chip"]:
		print("finished alignment, messaging slack...")
		os.environ["SLACK_API_TOKEN"] = "xoxp-455361134131-456511378663-517030405335-44bef1524e555ddd40f58c8744dc98d9"
		slack_token = os.environ["SLACK_API_TOKEN"]
		sc = slack.WebClient(slack_token)
		message = "finished doing "+ steps +" in "+ os.path.abspath(".").split("/")[-2] #this should usually work for reporting the project (won't work when subprojects exist )

		#sc.api_call(
		#  "chat.postMessage",
		#  channel="bioinformatics",
		#  as_user=False,
		#  username="alignment bot",
		#  text=message
		#)

	runLog.close()


######################################################################
## align(fastqDir, runType)                                         ##
## Description: aligns mirna or bulk mrna using STAR or kallisto.   ##
######################################################################
def align(fastqDir, targets, runType, threads):
	fileNum = 0
	print("creating sample directories in Alignments/")
	for i in targets:
		if i.find("_R1") > 1:
			outfile = getOutfileName(i)
			if outfile == "nosplitchar":
				print("Sample name formatting is required! (samples must contain \"__\" or \"_S\"), exiting...")
				sys.exit()
			try:
				os.mkdir(alignmentDir+"/"+outfile)
				# if aligner == "kallisto":
				#       os.mkdir("KallistoOuts/"+outfile+"/logs.txt") #don't need this anymore -- just going to need a single alignment dir named for whatever aligner beig used
			except OSError:
				print("couldn't make dir for "+alignmentDir+"/"+outfile)
				pass

	targetsToRemove = []
	for i in targets:
		if i.find("_R2") > 1:
			targetsToRemove.append(i)
	targets = [x for x in targets if x not in targetsToRemove]

	if adapterSeq == "":
		adapterString = " "
	else:
		adapterString = " --clip3pAdapterSeq "+adapterSeq

	for i in targets:
		fileNum+=1

		if steps == "mirna":
			print("curr file path: " + os.path.abspath(trimdfolder+"/"+i))
		else:	
			print("curr file path: " + os.path.abspath(fastqDir+"/"+i))
		if i.find("_R1") > 1:
			outfile = getOutfileName(i)

			print(outfile)
			print("Loading and aligning", i)
			print("This is file "+ str(fileNum) + " of " + str(len(targets)) + " total files.")

			if "trimd" in i: #this is when cutadapt was alreeady run and the fastqDir is trimmedfastqs 
				trimd_i = i
				if runType == "paired":
					trimd_j = i.replace("_R1", "_R2")
			else:
				trimd_i = i.split(".")[0]+"_trimd."+i.split(".")[1]+"."+i.split(".")[2]
				if runType=="paired":
					j=i.replace("_R1", "_R2")
					trimd_j = trimd_i.replace("_R1", "_R2")

			if i.endswith('.fastq.gz'):
				rfc = " --readFilesCommand "+unzipInstruct
			else:
				rfc = ""

			if steps == "mirna": #181220: this block shouldn't even be here; one the exceRpt small-rna pipeline is merged into this pipeline, we will only be doing small RNA block, which shoujld have conditionals that continue if the organism actually has a [mi/t/r/y/pi/circular]RNA database

				if runType=="single": #full small rna truseq 3' adapter = NNNNTGGAATTCTCGGGTG ; 5' adapter = NNNNGUUCAGAGUUCUACA ;  3' full: TGGAATTCTCGGGTGCCAAGG
					#run cutadapt on single fastq
					# parallelCutadapt(i)
					subprocess.call("STAR --runMode alignReads --runThreadN " + str(threads) + " --genomeDir "+star_miRNApath+rfc+" --readFilesIn " + trimdfolder + "/" + trimd_i  + " --outFileNamePrefix "+alignmentDir+"/"+outfile+"/ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1 --outFilterMatchNmin 15 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 2 --alignIntronMax 1", shell=True)

				elif runType=="paired":
					print(i + " " + j)
					# parallelCutadapt(i)
					# parallelCutadapt(j) #run cutadapt on paired fastqs
					subprocess.call("STAR --runMode alignReads --runThreadN " + str(threads) + " --genomeDir "+star_miRNApath+rfc+" --readFilesIn  " + trimdfolder + "/" + trimd_i  + " " + trimdfolder + "/" + trimd_j  + " --outFileNamePrefix "+alignmentDir+"/"+outfile+"/ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1 --outFilterMatchNmin 15 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 2 --alignIntronMax 1", shell=True)

					# subprocess.call("STAR --runThreadN " + threads + "--genomeDir "+starpath+" --readFilesCommand gzcat --outFilterType BySJout --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --alignIntronMin 20 --alignIntronMax 1 --alignMatesGapMax 1000000  --readFilesIn  "  + "../" + fastqDir + "/" + i  + " " + "../" + fastqDir + "/" + j  +  " --clip3pAdapterSeq TGGAATTCTC --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix Alignments/"+outfile+"/", shell=True)
					# subprocess.call("samtools sort -o Alignments/"+outfile+"/Aligned.sortedByCoord.out.bam -O bam -T Alignments/"+outfile+"/Aligned.sortedByCoord -@" + threads + "" +" Alignments/"+outfile+"/Aligned.out.bam", shell=True)
					# subprocess.call("rm  Alignments/"+outfile+"/Aligned.out.bam", shell=True)

			if steps == "smallrna": #this is going to be a lengthy one when it's actually done
				if runType=="single":
					#miRNA
					subprocess.call("STAR --runMode alignReads --runThreadN " + str(threads) + " --genomeDir "+star_miRNApath+rfc+" --readFilesIn " + trimdfolder + "/" + trimd_i  + " --outFileNamePrefix "+alignmentDir+"/"+outfile+"/ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1 --outFilterMatchNmin 15 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 2 --alignIntronMax 1", shell=True)
					#tRNA
					subprocess.call("STAR --runMode alignReads --runThreadN " + str(threads) + " --genomeDir "+star_tRNApath+rfc+" --readFilesIn " + trimdfolder + "/" + trimd_i  + " --outFileNamePrefix "+alignmentDir+"/"+outfile+"/ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1 --outFilterMatchNmin 15 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 2 --alignIntronMax 1", shell=True)
					#yRNA
					subprocess.call("STAR --runMode alignReads --runThreadN " + str(threads) + " --genomeDir "+star_yRNApath+rfc+" --readFilesIn " + trimdfolder + "/" + trimd_i  + " --outFileNamePrefix "+alignmentDir+"/"+outfile+"/ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1 --outFilterMatchNmin 15 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 2 --alignIntronMax 1", shell=True)
					#entire genome
					subprocess.call("STAR --runMode alignReads --runThreadN " + str(threads) + " --genomeDir "+starpath+rfc+" --readFilesIn " + trimdfolder + "/" + trimd_i  + " --outFileNamePrefix "+alignmentDir+"/"+outfile+"/ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1 --outFilterMatchNmin 15 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 2 --alignIntronMax 1", shell=True)
				elif runType=="paired":
					print(i + " " + j)
					#miRNA
					subprocess.call("STAR --runMode alignReads --runThreadN " + str(threads) + " --genomeDir "+star_miRNApath+rfc+" --readFilesIn  " + trimdfolder + "/" + trimd_i  + " " + trimdfolder + "/" + trimd_j  + " --outFileNamePrefix "+alignmentDir+"/"+outfile+"/ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1 --outFilterMatchNmin 15 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 2 --alignIntronMax 1", shell=True)
					#tRNA
					subprocess.call("STAR --runMode alignReads --runThreadN " + str(threads) + " --genomeDir "+star_tRNApath+rfc+" --readFilesIn  " + trimdfolder + "/" + trimd_i  + " " + trimdfolder + "/" + trimd_j  + " --outFileNamePrefix "+alignmentDir+"/"+outfile+"/ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1 --outFilterMatchNmin 15 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 2 --alignIntronMax 1", shell=True)
					#yRNA
					subprocess.call("STAR --runMode alignReads --runThreadN " + str(threads) + " --genomeDir "+star_yRNApath+rfc+" --readFilesIn  " + trimdfolder + "/" + trimd_i  + " " + trimdfolder + "/" + trimd_j  + " --outFileNamePrefix "+alignmentDir+"/"+outfile+"/ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1 --outFilterMatchNmin 15 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 2 --alignIntronMax 1", shell=True)
					#entire genome
					subprocess.call("STAR --runMode alignReads --runThreadN " + str(threads) + " --genomeDir "+starpath+rfc+" --readFilesIn  " + trimdfolder + "/" + trimd_i  + " " + trimdfolder + "/" + trimd_j  + " --outFileNamePrefix "+alignmentDir+"/"+outfile+"/ --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outFilterMultimapNmax 1 --outFilterMatchNmin 15 --outFilterMismatchNoverLmax 0.05 --alignIntronMin 2 --alignIntronMax 1", shell=True)

			if steps == "chip" or steps == "atac":
				if runType == "single":
					print(i)
					subprocess.call("STAR --runThreadN " + str(threads) + " --genomeDir "+starpath+rfc+" --outFilterType BySJout --outFilterScoreMinOverLread 0.25 --outFilterMatchNminOverLread 0.25 --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --alignIntronMin 20 --alignIntronMax 1 --alignMatesGapMax 10  --readFilesIn  " + fastqDir + "/" + i + adapterString + " --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix "+alignmentDir+"/"+outfile+"/", shell=True)
				elif runType == "paired":
					print(i + " " + j)
					subprocess.call("STAR --runThreadN " + str(threads) + " --genomeDir "+starpath+rfc+" --outFilterType BySJout --outFilterScoreMinOverLread 0.25 --outFilterMatchNminOverLread 0.25 --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --alignIntronMin 20 --alignIntronMax 1 --alignMatesGapMax 10  --readFilesIn  " + fastqDir + "/" + i + " " + fastqDir + "/" + j  + adapterString + " --limitBAMsortRAM 96000000000 --outSAMtype BAM Unsorted --quantMode GeneCounts --outFileNamePrefix "+alignmentDir+"/"+outfile+"/", shell=True)
					subprocess.call("samtools sort -o "+alignmentDir+"/"+outfile+"/Aligned.sortedByCoord.out.bam -O bam -T "+alignmentDir+"/"+outfile+"/Aligned.sortedByCoord -@ " + str(threads) + " "+alignmentDir+"/"+outfile+"/Aligned.out.bam", shell=True)
					subprocess.call("rm  "+alignmentDir+"/"+outfile+"/Aligned.out.bam", shell=True)

			if aligner == "STAR" and steps not in ["chip", "atac", "mirna", "smallrna"]: #this is default alignandcount behavior
				if runType == "single":
					print(i)
					subprocess.call("STAR --runThreadN " + str(threads) + " --genomeDir "+starpath+rfc+" --outFilterType BySJout --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --readFilesIn  " + fastqDir + "/" + i  + adapterString + " --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix "+alignmentDir+"/"+outfile+"/", shell=True)
				elif runType == "paired":
					print(i + " " + j)
					subprocess.call("STAR --runThreadN " + str(threads) + " --genomeDir "+starpath+rfc+" --outFilterType BySJout --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --readFilesIn  " + fastqDir + "/" + i  + " " + fastqDir + "/" + j  + adapterString + " --limitBAMsortRAM 96000000000 --outSAMtype BAM Unsorted --quantMode GeneCounts --outFileNamePrefix "+alignmentDir+"/"+outfile+"/", shell=True)
					subprocess.call("samtools sort -o "+alignmentDir+"/"+outfile+"/Aligned.sortedByCoord.out.bam -O bam -T "+alignmentDir+"/"+outfile+"/Aligned.sortedByCoord -@ " + str(threads) + " "+alignmentDir+"/"+outfile+"/Aligned.out.bam", shell=True)
					subprocess.call("rm  "+alignmentDir+"/"+outfile+"/Aligned.out.bam", shell=True)

			if aligner == "kallisto" and steps not in ["chip", "atac", "mirna", "smallrna"]:
				if runType == "single":
					if strand == 0:
						subprocess.call("kallisto quant --single -l 250 -s 100 -t " + str(threads) + " -i "+ kallistoIDX +" -b 100  -o "+alignmentDir+"/"+outfile+"/ " + fastqDir + "/" + i  + " > "+alignmentDir+"/"+outfile+"/logs.txt 2>&1", shell=True)
					if strand == 1:
						subprocess.call("kallisto quant --single -l 250 -s 100 --fr-stranded -t " + str(threads) + " -i "+ kallistoIDX +" -b 100  -o "+alignmentDir+"/"+outfile+"/ " + fastqDir + "/" + i + " > "+alignmentDir+"/"+outfile+"/logs.txt 2>&1" , shell=True)
					if strand == 2:
						subprocess.call("kallisto quant --single -l 250 -s 100 --rf-stranded -t " + str(threads) + " -i "+ kallistoIDX +" -b 100  -o "+alignmentDir+"/"+outfile+"/ " + fastqDir + "/" + i  + " > "+alignmentDir+"/"+outfile+"/logs.txt 2>&1", shell=True)

				elif runType == "paired":
					print(i + " " + j)
					if strand == 0:
						subprocess.call("kallisto quant -t " + str(threads) + " -i "+ kallistoIDX +" -b 100  -o "+alignmentDir+"/"+outfile+"/ " + fastqDir + "/" + i  + " " + fastqDir + "/" + j + " > "+alignmentDir+"/"+outfile+"/logs.txt 2>&1", shell=True)
					if strand == 1:
						subprocess.call("kallisto quant --fr-stranded -t " + str(threads) + " -i "+ kallistoIDX +" -b 100  -o "+alignmentDir+"/"+outfile+"/ " + fastqDir + "/" + i  + " "+ fastqDir + "/" + j + " > "+alignmentDir+"/"+outfile+"/logs.txt 2>&1", shell=True)
					if strand == 2:
						subprocess.call("kallisto quant --rf-stranded -t " + str(threads) + " -i "+ kallistoIDX +" -b 100  -o "+alignmentDir+"/"+outfile+"/ " + fastqDir + "/" + i  + " " + fastqDir + "/" + j + " > "+alignmentDir+"/"+outfile+"/logs.txt 2>&1", shell=True)
		else:
			print("didn't find any files with \'_R1\' in name") #it's probably just paired end (or none of the files have R1/R2, which needs to be fixed)


######################################################################
##counts()                                                              ##
##Description: Collects gene counts for STAR or kallisto alignments ##
######################################################################

## Nonstranded preps (e.g.Nugen Ovation or Illumina Truseq Unstranded): 0
## 1st strand synthesis: 1
## 2nd strand synthesis (e.g. Illumina Truseq Stranded): 2

def counts():
	#assuming we're already in AlignAllOut/ ... (this should probably be fixed)

	logs_counts=open("counts."+aligner+"."+genome+".txt",'w')
	logs_alignstats=open("Log."+aligner+"."+genome+".txt",'w')
	if aligner == "STAR":
		logs_minstats=open("stats."+aligner+"."+genome+".txt",'w')

	d={}
	qd=OrderedDict()
	samples=[]

	name = ""
	## loop through individual counts files and store them all in a single
	#should be for both kallisto and STAR, but kallisto doesn't really need it
	if aligner == "STAR":
		for i in os.listdir(alignmentDir):
			smalld={}
			if i.endswith("out"):
				name=i.rstrip('out')
				samples.append(name)
				x=open(alignmentDir+"/"+i+"/ReadsPerGene.out.tab", 'rU').readlines()
				for j in x:
					j=j.strip().split('\t')
					id=j[0]
					count=j[int(strand)+1]
					smalld[id]=count
				d[name]=smalld
				y=open(alignmentDir+"/"+i+"/Log.final.out", 'rU').readlines()
				for j in y:
					if len(j) > 1:
						j=j.strip(" ").strip("\n")
						j=j.split('\t')
						if len(j) > 1:
							value=j[1]
							id=j[0].strip(" |")
							if id not in qd:
								qd.setdefault(id, [])
								# value=j[1]
								qd[id].append(value)
							else:
								# value=j[1]
								qd[id].append(value)

		##write headers to counts and stats files
		logs_counts.write('Ensembl_Gene_ID\t'+'\t'.join(samples)+'\n')
		logs_alignstats.write('Stat\t'+'\t'.join(samples)+'\n')
		logs_minstats.write('Stat\t'+'\t'.join(samples)+'\n')

	elif aligner == "kallisto":
		for i in os.listdir(alignmentDir):
			smalld={}
			if i.endswith("out"):
				name=i.rstrip('out')
				samples.append(name)
				x=open(alignmentDir+"/"+i+"/abundance.tsv", 'rU').readlines()
				for j in x:
					j=j.strip().split('\t')
					id=j[0]
					if id == "target_id":
						continue
					count=j[EST_COUNTS]
					smalld[id]=count
				d[name]=smalld
				y=open(alignmentDir+"/"+i+"/run_info.json", 'rU')
				y_par=json.load(y)
				for key,value in y_par.items():
					if key == "call": #this is the log of what was run on commandline, don't need this in an excel spreadsheet...
						continue
					if key not in qd:
						qd.setdefault(key,[])
						qd[key].append(str(value))
					else:
						qd[key].append(str(value))

		logs_counts.write('Ensembl_Transcript_ID\t'+'\t'.join(samples)+'\n')
		logs_alignstats.write('Stat\t'+'\t'.join(samples)+'\n')

	## loop through dictionaries in ensembl ID order and write to files
	if name == "":
		usercontinue = input("WARNING: Alignment dir empty in folder: " + os.path.dirname("../Alignments") + " \nNote: this warning usually results in fatal program error.\nAre you sure you want to continue? (Y/N):").upper()
		if usercontinue == 'Y':
			print("Running...")
		else:
			print("Exiting...")
			sys.exit()

	for key in sorted(d[name].keys()):
		line=[key]
		for i in samples:
			line.append(d[i][key])
		if key.startswith("N_") and aligner == "STAR":
			logs_minstats.write('\t'.join(line)+'\n')
		else:
			logs_counts.write('\t'.join(line)+'\n')
	for key, values in qd.items():
		line = [key]
		line.extend(qd[key])
		logs_alignstats.write('\t'.join(line)+'\n')

	logs_counts.close()
	logs_alignstats.close()
	if aligner == "STAR":
		logs_minstats.close()


def counts_miRNA():

	all_counts_file = "counts.miRNA."+genome+".txt"
	miRNA_counts_file = alignmentDir+"/"+all_counts_file
	samples = list(filter(lambda x: x != all_counts_file, os.listdir(alignmentDir)))

	if all_counts_file in os.listdir(alignmentDir):
		usercontinue = input(all_counts_file+" already exists, do you want to overwrite? (y/n): ")
		if usercontinue.lower() ==  'y':
			print("ok, removing and creating new miRNA counts file...")
			subprocess.call("rm "+alignmentDir+"/"+all_counts_file, shell=True)
			for s in samples:
				subprocess.call("rm "+alignmentDir+"/"+s+"/miRNA_counts.txt", shell=True)
		else:
			print("ok, avoiding overwrite.  exiting...")
			sys.exit()

	sortedGenes = []
	with open(miRNA_genelist, 'rU') as rgf:
		for l in rgf:
			sortedGenes.append(l.rstrip("\n"))

	if aligner != "STAR":
		print("aligner must be STAR, re-run this step after generating STAR alignments")
		sys.exit()

	#STEP 1: can't just extract out the gene alignments, so we ahve to make a list of them in a file before we create a variable -> takes a while, need parallelization
	pool = multiprocessing.Pool(int(threads))
	pool.map(extract_miRNA_alignments, samples)
	pool.close()
	pool.join()

	#STEP2: collect all the gene names/counts into a dictionary of dictionaries
	print("creating dictionary for recorded reads..")
	d = {}
	qd = OrderedDict()
	for s in samples:
		smalld = {}
		with open(alignmentDir+"/"+s+"/miRNA_counts.txt", 'rU') as readsFile:
			for line in readsFile:
				gene = line.rstrip("\n").split("\t")
				smalld[gene[0]] = gene[1] 
			d[s] = smalld
		y=open(alignmentDir+"/"+s+"/Log.final.out", 'rU').readlines()
		for j in y:
			if len(j) > 1:
				j=j.strip(" ").strip("\n")
				j=j.split('\t')
				if len(j) > 1:
					value=j[1]
					id=j[0].strip(" |")
					if id not in qd:
						qd.setdefault(id, [])
						value=j[1]
						qd[id].append(value)
					else:
						value=j[1]
						qd[id].append(value)


	#STEP3: iterate through all sortec miRNA and if there are any counts in dictionary for each sample, output to file
	print("writing counts file...")
	header = "\t"+"\t".join(samples)+"\n"
	with open(miRNA_counts_file, 'w') as mcf: #miRNA_counts_file
		mcf.write(header)
		for g in sortedGenes:
			temp = []
			temp.append(g)
			for s in samples:
				if g in d[s].keys():
					temp.append(str(d[s][g]))
				else:
					temp.append('0')
			mcf.write("\t".join(temp) +"\n")

	#last step: create summarized log file
	with open("Log."+aligner+"."+genome+".txt",'w') as logfile:
		logfile.write('Stat\t'+'\t'.join(samples)+'\n')
		for key, values in qd.items():
			line = [key]
			line.extend(qd[key])
			logfile.write('\t'.join(line)+'\n')


def extract_miRNA_alignments(s):
	print("extracting aligned genes from bam of sample: "+s) #260 flag isn't really necessary I think beacuse 256 is multi-aligned and 4 is paired or something
	subprocess.call(""" samtools view -F 260 """+alignmentDir+"""/"""+s+"""/Aligned.out.bam | cut -f 3 | sort | uniq -c | awk '{printf "%s\\t%s\\n", $2, $1}' > """+alignmentDir+"""/"""+s+"""/miRNA_counts.txt """, shell=True)
	# subprocess.call("samtools view -F 260 "+alignmentDir+"/"+s+"/Aligned.out.bam | cut -f 3 | sort | uniq -c | awk '{printf(\"%s\t%s\n\", $2, $1)}' > "+alignmentDir+"/"+s+"/miRNA_counts.txt", shell=True)

# execute a single cutadapt thread
def exCutadapt(processAndFile):
	process = processAndFile.split(" ")[0]
	file = processAndFile.split(" ")[1]
	trimmed_file = file.split(".")[0]+"_trimd."+file.split(".")[1]+"."+file.split(".")[2]
	outfile = getOutfileName(file)
	adapterSeq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" 	

	try: 
		print("trimming file process:  "+process) ### need to add checking of known addaptor sequences as a prereq to this code... illumina: NNNNTGGAATTCTCGGGTGCCAAGG; Lexogen: CCGCCGGGCAGCTTCCGGGAA
		print("cutadapt -m 10 -a "+adapterSeq+" -o trimdfastqs/"+trimmed_file+" "+fastqDir+"/"+file+" > cutadapt_logs/"+outfile+"/cutadapt_logs.txt")
		subprocess.call("cutadapt -m 10 -a "+adapterSeq+" -o trimdfastqs/"+trimmed_file+" "+fastqDir+"/"+file+" > cutadapt_logs/"+outfile+"/cutadapt_logs.txt", shell=True)

	except OSError: 
		print("well that didn't work")
	except: 
		print("it wasn't an OSError")

# execute all threads of a cutadapt run using exCutadapt
def cutadaptAll(files, t=6):
	print("using Illumina TruSeq small rna adapter for trimming...")
	outfiles = [getOutfileName(file) for file in files]
	for i in outfiles:
		try:
			os.mkdir("cutadapt_logs/"+i)
		except OSError:
			pass
	files_enum = zip(range(1,len(files)+1), files)
	files_enum = [" ".join([str(x[0]), x[1]]) for x in files_enum] #this is literally just for fancy process numbers
	pool = multiprocessing.Pool(t) # let's just split all fastqs into 4 chunks (n=4)
	pool.map(exCutadapt, files_enum)
	pool.close()
	pool.join()

########################
# picard rna seq metrics
########################

def runPicardRNAMetrics():
	iList = []
	for i in os.listdir(alignmentDir):
		if i.endswith("out"):
			iList.append(i)
	if strand == 0:
		picardStrand = "NONE"
	elif strand == 1:
		picardStrand = "FIRST_READ_TRANSCRIPTION_STRAND"
	elif strand == 2:
		picardStrand = "SECOND_READ_TRANSCRIPTION_STRAND"
	else:
		print("You have an error in your strand ID.")
		print("Choose 0, 1, or 2.")
	pool = multiprocessing.Pool(int(threads)) # run this many threads concurrently
	func = partial(exPicardMetrics, refFlat, riboIntervals, picardStrand)
	pool.map(func, iList)
	pool.close()
	pool.join()

def exPicardMetrics(refFlat, riboIntervals, strand, i):
	name = i.strip("out")
	if riboIntervals != "":
		subprocess.call("samtools view -H "+alignmentDir+"/"+i+"/Aligned.sortedByCoord.out.bam > cur_ribo_intervals.txt; cat "+riboIntervals+" >> cur_ribo_intervals.txt", shell=True)
		riboCode = "RIBOSOMAL_INTERVALS=cur_ribo_intervals.txt"
	else:
		riboCode = ""
	subprocess.call("java -Xmx10g -jar "+picard_jarfile+" CollectRnaSeqMetrics I="+alignmentDir+"/"+i+"/Aligned.sortedByCoord.out.bam \
		TMP_DIR="+alignmentDir+"/"+i+" REF_FLAT="+refFlat+" STRAND_SPECIFICITY="+strand+" "+riboCode+" QUIET=true OUTPUT="+alignmentDir+
		"/"+i+"/picard.txt 2>"+alignmentDir+"/"+i+"/picard.err", shell=True)

################################
#### utility definitions #######
################################

def runMode(targetsList):
	runType = "single"
	for i in targetsList:
		if i.find("_R2") > 1:
			runType = "paired"
			#targetsList.remove(i)
		continue
	return runType

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

def combineLanes(runtype, targets, s):
	if runtype == "paired":
		paired = True
	else: paired = False

	L000file1 = fastqDir+"/"+s+"_L000_R1_001.fastq.gz"
	if paired:
		L000file2 = fastqDir+"/"+s+"_L000_R2_001.fastq.gz"
	
	sampleList = [fastqDir+"/"+t for t in targets if s in t]
	sampleList_r1 = [t for t in sampleList if "_R1_" in t]
	sampleList_r2 = [t for t in sampleList if "_R2_" in t]

	if len(sampleList_r1) > 0:
		catList1 = " ".join(sampleList_r1)
		print("cat "+catList1+" > "+L000file1)
		subprocess.call("cat "+catList1+" > "+L000file1, shell=True)

	if len(sampleList_r2) > 0:
		catList2 = " ".join(sampleList_r2)
		print("cat "+catList2+" > "+L000file2)
		subprocess.call("cat "+catList2+" > "+L000file2, shell=True)
		
# alignedSamples = [s for s in ad if "Aligned.sortedByCoord.out.bam" in os.listdir(s)]		
def alignmentAlreadyDone(r1Targets):
	for ad in os.listdir("Alignments/"):
		alignedSamples = []
		for s in os.listdir("Alignments/"+ad):	
			if "Aligned.sortedByCoord.out.bam" in os.listdir("Alignments/"+ad+"/"+s):
				alignedSamples.append(s)
		if len(alignedSamples) == len(r1Targets):
			return ad
	return None


# def multiprocessMD5s(threads):
# 	total = len(os.listdir(fastqDir))
# 	names = ["md5sum"]
# 	pool = multiprocessing.Pool(int(threads)) # run this many threads concurrently
# 	pool.map(generateMD5s, names)
# 	pool.close()
# 	pool.join()	

#makes a text file in AlignAllOut that has all of the md5 hashes for fastq files
def generateMD5s():
	md5call = "md5" if machine == "LBCJoshPollack2.ucsf.edu" else "md5sum"
	for f in os.listdir(fastqDir):
		subprocess.call(md5call+" "+fastqDir+"/"+f+" >> md5sums.txt", shell=True)



# this was mostly written down in order to track what I did to make the rRNA interval for hg38
# convert a gtf to a ribosomal rna interval file and format the new file
def getrRNAIntervals(input_gtf, output_bed):
	if not os.path.isfile(input_gtf):
		print("GTF file not found...")
		return False 
	gtf = pr.get_gtf(input_gtf)
	gtf_rrna = gtf[(gtf.gene_biotype == "rRNA") | (gtf.gene_biotype == "Mt_rRNA")]
	gtf_rrna_g = gtf_rrna[(gtf_rrna.Feature == "gene")]
	gtf_rrna_g.to_bed(output_bed)
	temp_bed = output_bed.replace(".bed", "_ts.bed")
	final_bed = output_bed.replace("ts.bed", "_f.bed")
	final_out = "final_out.txt" 
	subprocess.call("tr -s '\t' '\t' < "+output_bed+" > "+temp_bed, shell=True)
	subprocess.call("cut -f 1,2,3,6,9 "+temp_bed+" > "+final_bed, shell=True)
	subprocess.call("""samtools view -H """+sample_bam+""" | awk 'BEGIN{FS=OFS="\t"}; $1 ~ /^@SQ/ {print $0}') """+
		"""> """+final_out, shell=True)
	return True

# object for organizing fastq files that need to be run in alignment
# one target can have an R1, R2, I1, and I2 file, but for SE50 will probably only have a R1 file
# directory is absolute path to that fastq file, so should just be 'fastqDir' 
class Target: 
	def __init__(self, sname, directory = None, R1 = None, R2 = None, I1 = None, I2 = None):
		self.sname = sname
		if directory is not None:
			self.directory = directory
		if R1 is not None:
			self.R1 = R1
		if R2 is not None:
			self.R2 = R2
		if I1 is not None:
			self.I1 = I1
		if I2 is not None:
			self.I2 = I2		
	def assignDir(self, dir_):
		self.directory = dir_
	def assignR1(self, R1):
		self.R1 = R1
	def assignR2(self, R2):
		self.R2 = R2
	def assignI1(self, I1):
		self.I1 = I1
	def assignI2(self, I2):
		self.I2 = I2

if __name__ == "__main__":
	main(sys.argv[1:])


