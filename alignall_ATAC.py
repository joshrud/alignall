import subprocess, os, string
import multiprocessing
import csv
import math
import getopt
import sys
import os.path
import argparse
import re
import gzip
from functools import partial
from collections import OrderedDict
from slackclient import SlackClient

#now import auxiliary files
# import alignall_atacOnly	#UNTESTED
# import mergeVcfs			#UNTESTED
# import callSNPs				#UNTESTED
# import alignall_chipOnly	#UNTESTED

# from parseST import parseST
# import xlrd

#from normalize_seq_depth_and_signal_to_noise.defs import CountAlignedReads, BedReads, BedSlop, makeBedGraph, makeNormBedGraph, makeNormBigWig

#  CHANGE LOG: 
# <date> | <change_in_detail>
# ------------------------------
# 190214 | Edited with intention of making it a standalone script
# 190628 | Edited with intention of making this a plugin to alignall
#
#
#
#
#
VERSION = 1 

###### plug in for performing ATAC-seq peak calling in alignall 
###### basically just need some functions 

## intended to be run within project area, just like AlignAllOut
## ASSUMPTIONS: 1. aligned reads are in AlignAllOut/Alignments/* 
## 				2. bam files are all named Aligned.sortedByCoord.out.bam; and were aligned with STAR
	

####   Function Definitiions   ################################################                                 
###############################################################################


GATKbase="java -jar ~/tools/GATK/GenomeAnalysisTK.jar -T "

#### getConditionGroups
# For 

def getConditionGroups(alignmentDir):
	#if duplicate names in a single sample, last position is recorded only
	finalGrouping = {}
	allSamples = [x.rstrip("out") for x in os.listdir(alignmentDir) if bool(re.search("_[0-9]out$", x))]
	splitNames = {} #example: {0: ['ATAC'], 1: ['CP101', 'JD130'], 2: ['Q111', 'Q20'], 3: ['1', '2', '3']} #for more complex conditions
	for s in allSamples:

		#add to splitNames
		length = len(s.split("_"))
		for i in range(length):
			if i in splitNames.keys():
				if s.split("_")[i] not in splitNames[i]: #works like dictionary
					splitNames[i].append(s.split("_")[i])
			else:
				splitNames[i] = [s.split("_")[i]]

	#remove non-conditions
	for i in splitNames.keys():
		if len(splitNames[i]) == 1:
			del splitNames[i] #remove the redundant sequence entry
			continue
		try:
			temp = [int(x) for x in splitNames[i]]
		except ValueError:
			pass
		else:
			del splitNames[i]

	#now all complex groupings (across conditions) ex: {CP101: [CP101_Q20_1out, CP101_Q20_2out, CP101_Q111_2out, CP101_Q175_1out, ......], .....}
	for k in splitNames.keys():
		for v in splitNames[k]:
			finalGrouping[v] = [x+"out" for x in allSamples if v in x]

	return finalGrouping

def getBasenames(alignmentDir):
	#if duplicate names in a single sample, last position is recorded only
	finalGrouping = {}
	allSamples = [x.rstrip("out") for x in os.listdir("Align") if bool(re.search("_[0-9]out$", x))]
	basenames = {}
	for s in allSamples: 
		#add to basenames
		cur_basename = re.sub("_[0-9]$","",s)
		if cur_basename in basenames:
			basenames[cur_basename] +=1
		else:
			basenames[cur_basename] = 1 

	#add all groupings of sample names for basenames ex: {CP101_Q20: [CP101_Q20_1out, CP101_Q20_2out, CP101_Q20_3out], ...}
	basenames_group = []
	for k in basenames.keys():
		basenames_group.append(k)
		finalGrouping[k] = [x+"out" for x in allSamples if k in x] #add directory-confirmed samplenames back into dictionary as values of basename key

	return finalGrouping

def runRmDups(alignmentDir, threads):
	iList = []
	for i in os.listdir(alignmentDir):
		if i.endswith("out"):
			iList.append(i)
	pool = multiprocessing.Pool(int(threads)) # on 3 processors
	pool.map(exRmDups, iList)
	pool.close()
	pool.join() 

def exRmDups(j):
	name = j.strip("out")
	subprocess.call("java -jar "+picard_jarfile+" MarkDuplicates I="+alignmentDir++j+"/Aligned.sortedByCoord.out.bam O="+alignmentDir+j+"/Aligned.sortedByCoord.duplicateRemoved.out.bam M="+alignmentDir+j+"/marked_dup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true > "+alignmentDir+j+"/picard.log.txt 2> "+alignmentDir+j+"/picard.err.txt", shell=True)

def runCallPeaks(): 
	samples = [x.rstrip("out") for x in os.listdir(alignmentDir) if bool(re.search("_[0-9]out$", x))]
	if len(samples) == 0:
		print "aligned samples don't seem to follow formatting requirements, please re-format. "
		sys.exit()

	for i in samples:
		tempname = i.rstrip("out")
		if "_" not in tempname:
			print("need to choose a different character/string, or just reformat samples; exiting...")
			sys.exit()

	#make all the main sample directories in preparation for merge macs
	basenames = getBasenames()
	if "peakCalling" not in os.listdir("."):
		os.mkdir("peakCalling")
	else:
		user_continue = raw_input("peakCalling already exists; do you want to overwrite?(type 'y' to proceed): ").lower()
		if user_continue != 'y':
			print "ok, quitting..."
			sys.exit()
	for k in basenames.keys():
		if k+"_rep0" not in [x for x in os.listdir(alignmentDir) if os.path.isdir(os.path.join(alignmentDir, x))]:
			os.mkdir(os.path.join(alignmentDir,k+"_rep0"))
		if k+"_rep0" not in [x for x in os.listdir("peakCalling") if os.path.isdir(x)]:
			os.mkdir(os.path.join("peakCalling",k))

	chipMasterList = []
	for key in basenames.keys():
		groupString = key+" "+" ".join([""+alignmentDir+x+"/Aligned.sortedByCoord.duplicateRemoved.out.bam" for x in basenames[key]]) #need to make the key the 1st in the list and then take it for recognizing dir. when passed into merged macs
		chipMasterList.append(groupString)

	print ""
	print "Running ChIP-seq peak calling on merged replicates"
	macsMergedThreads = int(threads)
	pool = multiprocessing.Pool(int(macsMergedThreads)) # run this many threads concurrently
	func = partial(exMergedMacs, macsGenome) #don't need the inputs for atac 
	pool.map(func, chipMasterList) #master list is all non-control samples 							## if we can feed it a key,value pair, then this would work a lot better; (we'd actually be able to find the directory)
	pool.close()
	pool.join()
	
	print ""
	print "Merging replicates on condition"
	try:
		os.mkdir("peakCalling/FinalMergedPeaks")
	except OSError:
		print "unable to make dir: peakCalling/FinalMergedPeaks, continuing..." 
		pass

	print "merging peak beds on condition.."
	for i in conditionGroups.keys():
		print "cat peakCalling/*"+i+"*/Merged_reps_pass_idr.broadPeak | sort -k1,1 -k2,2n | bedtools merge -i - > peakCalling/FinalMergedPeaks/"+i+".merged.peaks.bed"
		subprocess.call("cat peakCalling/*"+i+"*/Merged_reps_pass_idr.broadPeak | sort -k1,1 -k2,2n | bedtools merge -i - > peakCalling/FinalMergedPeaks/"+i+".merged.peaks.bed", shell=True)
		
	print ""
	print "Partitioning peaks between conditions"
	subprocess.call("bedops --partition peakCalling/FinalMergedPeaks/*.merged.peaks.bed | awk -v OFS=$'\t' '{print $1,$2,$3,\"peak_\"FNR}' > peakCalling/FinalMergedPeaks/Partitioned_peaks.bed",shell=True)
	print ""
	print "Converting bed to gtf"
	subprocess.call("awk -v OFS=$'\t' '{print $1,\"ChIP\",\"peak\",$2,$3,\".\",\".\",\".\",\"gene_id \\\"\",$4,\"\\\";\"}' peakCalling/FinalMergedPeaks/Partitioned_peaks.bed | awk -v OFS=$'\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9\" \"$10$11$12}' > peakCalling/FinalMergedPeaks/Partitioned_peaks.gtf", shell = True)		
	
	interCommands = []
	for j in os.listdir("peakCalling"):
		if j.endswith("rep0"):
			interCommand = " | intersectBed -a - -b peakCalling/"+j+"/Merged_reps_pass_idr.sort.broadPeak -wa -u -f .5 "
			if interCommand not in interCommands:
				interCommands.append(interCommand)
	interCommandsString = " ".join(interCommands)

	subprocess.call("cat peakCalling/FinalMergedPeaks/Partitioned_peaks.bed " + interCommandsString + " > peakCalling/FinalMergedPeaks/Partitioned_peaks.common.bed", shell = True)
	
#####################
## runs macs2 callpeak on an individual chip sample
## macsGenome: the 2-character representation of the genome of interest (hg, mm, gg, etc...)
## INPUTS: the string corresonding to all control/input files for the chip project
## CHIPS: the replicates that all correspond to the same sample
##
##
#####################
def exMergedMacs(macsGenome, CHIPS): #from alignall v10

 	groupName = CHIPS.split(" ")[0]
 	chipsFiles = " ".join(CHIPS.split(" ")[1:])
 	chipsNames = [x.split("/")[1] for x in chipsFiles.split(" ")]

	groupDIR = groupName+"_rep0"
	mergeBam = groupName+"_rep0.bam"
	
	print "merging all bams of condition: " +groupName
	if os.path.exists(""+alignmentDir++groupDIR+"/"+mergeBam) and os.path.getsize()>0:
		print "looks like this merged bam already exists, skipping this step..."
	else:
		subprocess.call("samtools merge -@ 4 "+alignmentDir+ + groupDIR + "/" + mergeBam + " " + chipsFiles, shell = True)
	
	pseudo1DIR = groupName+"_pseudo1"
	pseudo2DIR = groupName+"_pseudo2"
	pseudo1Bam = pseudo1DIR + ".bam"
	pseudo2Bam = pseudo2DIR + ".bam"
	try: #make the directories
		os.mkdir(""+alignmentDir++pseudo1DIR)
		os.mkdir(""+alignmentDir++pseudo2DIR)
	except OSError:
		print "didn't make pseudo dirs for samples in group: " + groupName
		pass
		# user_continue = raw_input("ERROR: wasn't able to make the directories for pseudo-samples in alignments.\nIf these directories exist in Alignments/ , hit enter to continue. -> ")
		# if user_continue != "": 
		# 	print "Exiting..."
		# 	sys.exit()
	
	#1.5 and 2.5 are seeds used to generate the merged pseudo bam files
	subprocess.call("samtools view -@ 4 -b -s 1.5 -o "+alignmentDir+ +pseudo1DIR+ "/" + pseudo1Bam + " "+alignmentDir+ + groupDIR + "/" + mergeBam, shell = True)
	subprocess.call("samtools view -@ 4 -b -s 2.5 -o "+alignmentDir+ +pseudo2DIR+ "/" + pseudo2Bam + " "+alignmentDir+ + groupDIR + "/" + mergeBam, shell = True)
		
	MACSOUTBASE="peakCalling/"+groupDIR
	MACSOUT="peakCalling/"+groupDIR+"/Merged_peaks_for_idr" #this name might need an appended "_peaks" 
	MACSOUTBASE_pseudo1="peakCalling/"+pseudo1DIR
	MACSOUT_pseudo1="peakCalling/"+pseudo1DIR+"/pseudo_rep_peaks_for_idr"
	MACSOUTBASE_pseudo2="peakCalling/"+pseudo2DIR
	MACSOUT_pseudo2="peakCalling/"+pseudo2DIR+"/pseudo_rep_peaks_for_idr"
		
	try:
		os.mkdir(MACSOUTBASE)
		os.mkdir(MACSOUTBASE_pseudo1)
		os.mkdir(MACSOUTBASE_pseudo2)
	except OSError:
		print "unable to make peak calling dirs for samples in group: " + groupName
		pass

	#call peaks on pseudo replicates
	subprocess.call("macs2 callpeak -t "+alignmentDir+ +pseudo1DIR+ "/" + pseudo1Bam + " -f BAM -g " + macsGenome + " --keep-dup 1 --bw 500 -n " + MACSOUT_pseudo1 + " --nomodel --extsize 300 --slocal 2000 --llocal 20000 -p 0.001 --broad > "+ MACSOUTBASE_pseudo1+"/macs2.callpeak.log.txt 2>"+ MACSOUTBASE_pseudo1+"/macs2.callpeak.err.txt", shell = True)
	subprocess.call("macs2 callpeak -t "+alignmentDir+ +pseudo2DIR+ "/" + pseudo2Bam + " -f BAM -g " + macsGenome + " --keep-dup 1 --bw 500 -n " + MACSOUT_pseudo2 + " --nomodel --extsize 300 --slocal 2000 --llocal 20000 -p 0.001 --broad > "+ MACSOUTBASE_pseudo2+"/macs2.callpeak.log.txt 2>"+ MACSOUTBASE_pseudo2+"/macs2.callpeak.err.txt", shell = True)
	#call peaks on actual replicates
	subprocess.call("macs2 callpeak -t " + chipsFiles + " -f BAM -g " + macsGenome + " --keep-dup 1 --bw 500 -n " + MACSOUT + " --nomodel --extsize 300 --slocal 2000 --llocal 20000 -p 0.001 --broad > "+ MACSOUTBASE+"/macs2.callpeak.log.txt 2>"+ MACSOUTBASE+"/macs2.callpeak.err.txt", shell = True)
	
	#perform IDR statistics for called peaks
	print "Running idr pipeline on : " + groupName
	print "peakCalling/"+groupDIR+"/Merged_reps_pass_idr.broadPeak"	
	subprocess.call("idr --samples "+MACSOUT_pseudo1+"_peaks.broadPeak "+MACSOUT_pseudo2+"_peaks.broadPeak --peak-list peakCalling/"+groupDIR+"/Merged_peaks_for_idr_peaks.broadPeak --plot --idr-threshold 0.01 --input-file-type broadPeak --output-file peakCalling/"+groupDIR+"/Merged_reps_pass_idr.broadPeak > peakCalling/"+groupDIR+"/idr.log.txt 2> peakCalling/"+groupDIR+"/idr.err.txt", shell = True)
	print "idr --samples "+MACSOUT_pseudo1+".broadPeak "+MACSOUT_pseudo2+".broadPeak --peak-list peakCalling/"+groupDIR+"/Merged_peaks_for_idr_peaks.broadPeak --plot  --idr-threshold 0.1 --input-file-type broadPeak --output-file peakCalling/"+groupDIR+"/Merged_reps_pass_idr.broadPeak > peakCalling/"+groupDIR+"/idr.log.txt 2> peakCalling/"+groupDIR+"/idr.err.txt"
	subprocess.call("sort -k1,1 -k2,2n peakCalling/"+groupDIR+"/Merged_reps_pass_idr.broadPeak > peakCalling/"+groupDIR+"/Merged_reps_pass_idr.sort.broadPeak", shell = True)
    
def runChipHTseq():
	#fix any possible GTF errors (happens for non-human I think)
	subprocess.call("awk -v OFS=\"\t\" '$4==0{$4=1};1' peakCalling/FinalMergedPeaks/Partitioned_peaks.gtf > peakCalling/FinalMergedPeaks/Partitioned_peaks_fixed.gtf")
	#closestBed -D b -a Partitioned_peaks_fixed.bed -b /Volumes/Pegasus/Genomes/Mus_musculus.GRCm38.93.tss.bed > peaks_near_tss.txt
	#bedtools window -a Partitioned_peaks_fixed.bed -b /Volumes/Pegasus/Genomes/Homo_sapiens.GRCh38.93.gtf -w 0 > features_by_peak.txt #this one is for overlapping genes/features, but we want the above( closeset TSS)
	#awk '$7=="gene" {gsub("[;\"]", "", $14); gsub("[;\"]", "", $18); print $1,$2,$3,$4,$7,$8,$9,$11,$14,$18}' features_by_peak.txt > genes_overlapping_peaks.txt
	#cut -f10 -d " " genes_overlapping_peaks.txt | sort | uniq -c | sort -k1,1 -nr | head 

	chipGTF = "peakCalling/FinalMergedPeaks/Partitioned_peaks_fixed.gtf" #if not human genome, 0s can be made as location i.e. start: 0 end: 1000
	iList = []
	### check if >100 mapped reads before running HTseq
	for i in os.listdir(alignmentDir):
		if i.endswith("out"):
			z=open(""+alignmentDir++i+"/Log.final.out", 'rU').readlines()
			for k in z:
				if len(k) > 1:
					k=k.strip(" ").strip("\n")
					k=k.split('\t')
					id = k[0]
					id = id.strip(" |")
					if id == "Uniquely mapped reads number":
						value = k[1]
			if value >= 100:
				iList.append(i)
	HTseqThreads = int(threads)*2
	pool = multiprocessing.Pool(int(HTseqThreads)) # run this many threads concurrently
	func = partial(exChipHTseq, htSeqStrand, chipGTF)
	pool.map(func, iList)
	pool.close()
	pool.join()

def exChipHTseq(htSeqStrand, gtf, j):
	subprocess.call("htseq-count -f bam -r pos -s " + htSeqStrand + " -t peak -i gene_id -m intersection-nonempty  "+alignmentDir++j+"/Aligned.sortedByCoord.duplicateRemoved.out.bam " + gtf + " > "+alignmentDir++j+"/ATAC.duplicateRemoved.counts.txt 2>  "+alignmentDir++j+"/ATAC.HTseq.duplicateRemoved.log.txt", shell=True)

def ChipCounts():
	
	newfile4=open("counts.ATAC."+genome+".txt",'w')
	d=OrderedDict()
	samples=[]

	name = ""
	for i in os.listdir(alignmentDir):
		smalld=OrderedDict()
		if i.endswith("out"):
			name=i.strip('out')
			HTseqFile = ""+alignmentDir++i+"/ATAC.duplicateRemoved.counts.txt"
			if os.path.isfile(HTseqFile):
				if os.stat(HTseqFile).st_size > 0:
					samples.append(name)
					x=open(HTseqFile, 'rU').readlines()
					for j in x:
						j=j.strip().split('\t')
						id=j[0]
						count=j[1]
						if id.find("__") != 0:
							smalld[id]=count
						d[name]=smalld
	
	newfile4.write('Peak_ID\t'+'\t'.join(samples)+'\n')

	for key in sorted(d[name].iterkeys()):
		line=[key]
		for i in samples:
			line.append(d[i][key])
		newfile4.write('\t'.join(line)+'\n')

def runMakeChipGraphs():
	chipList = []
	for i in os.listdir(alignmentDir):
		if i.endswith("out"):
			bam = ""+alignmentDir++i+"/Aligned.sortedByCoord.duplicateRemoved.out.bam"
			chipList.append(bam)
	
	pool = multiprocessing.Pool(int(threads)) # run this many threads concurrently
	func = partial(exMakeChipGraphs, refIdx, GenomeWindows)
	pool.map(func, chipList)
	pool.close()
	pool.join()

def exMakeChipGraphs(refIdx, GenomeWindows, inBam):
	base=inBam.split("/")[1]
	baseName=base.replace("out","")
	BedGraphNormFile=""+alignmentDir++base+"/"+baseName+".bg"
	BigWigNormFile="graphs/"+baseName+".bw"
	AlignedReads = CountAlignedReads(inBam)
	BedReads(inBam)
	BedSlop(inBam, refIdx)
	makeBedGraph(inBam, GenomeWindows)
	makeNormBedGraph(inBam, AlignedReads, BedGraphNormFile)
	makeNormBigWig(inBam, refIdx, BedGraphNormFile, BigWigNormFile)

def CountAlignedReads(inBam):
	inBamBase = inBam.split(".bam")[0]
	CountAlignedReadsFile=inBamBase+".AlignedReads.txt"
	if os.path.isfile(CountAlignedReadsFile):
		print "Using aligned read counts from previous file for " + inBam + " \n"
	else:
		CountAlignedReadsCommand = "samtools view -c -F 4 "+ inBam +" > " + CountAlignedReadsFile
		subprocess.call(CountAlignedReadsCommand, shell=True)
	inputReadA = open(CountAlignedReadsFile, 'rb')
	inputReadCountsA = list(csv.reader(inputReadA, delimiter = ' '))
	AllReads = float(inputReadCountsA[0][0])
	return AllReads

def BedReads(inBam):
	inBamBase = inBam.split(".bam")[0]
	BedReadFile=inBamBase+".bed"
	if os.path.isfile(BedReadFile):
		print "Using previous file for bam to bed conversion: "+ BedReadFile + " \n"
	else:
		BamToBedCommand = "bamToBed -i "+ inBam + " >" + BedReadFile
		subprocess.call(BamToBedCommand, shell=True)

def BedSlop(inBam, GenomeIndex):
	inBamBase = inBam.split(".bam")[0]
	BedSloppedFile=inBamBase+".slopped.bed"
	inBed=inBamBase+".bed"
	if os.path.isfile(BedSloppedFile):
		print "Using previous file for slopped bed: "+ BedSloppedFile + " \n"
	else:
		SlopBedCommand = "slopBed -i " + inBed +" -g " + GenomeIndex + " -l 0 -r 300  -s > "+ BedSloppedFile
		subprocess.call(SlopBedCommand, shell=True)

def makeBedGraph(inBam, GenomeWindows):
	inBamBase = inBam.split(".bam")[0]
	BedGraphFile=inBamBase+".slopped.bg"
	BedSlopped=inBamBase+".slopped.bed" 
	if os.path.isfile(BedGraphFile):
		print "Using previous bed graph file: ", BedGraphFile + "  \n"
	else:
		BedGraphCommand="coverageBed -counts -b " + BedSlopped + " -a "+ GenomeWindows + " | sort -k1,1 -k2,2n > "+ BedGraphFile
		subprocess.call(BedGraphCommand, shell=True)

def makeNormBedGraph(inBam, AllReads, BedGraphNormFile):
	inBamBase = inBam.split(".bam")[0]
	BedGraphFile= inBamBase+".slopped.bg"
	#BedGraphNormFile=inBamBase+".slopped.normDepthNoise.bg"
	if os.path.isfile(BedGraphNormFile):
		print "Output file " + BedGraphNormFile + " already exists  \n"
	else:
		row=0
		BedGraphOpen = open(BedGraphFile, 'rb')
		BedGraphList = csv.reader(BedGraphOpen, delimiter = '\t')
		covList = [((float(row[3]))/ float(AllReads))*(1000000.0)*(1000.0/float(float(row[2]) - float(row[1]))) for row in BedGraphList]
		covListPerc99P99 = str(np.percentile(covList, 99.99))
		covListPerc99P9 = str(np.percentile(covList, 99.9))
		covListPerc99 = str(np.percentile(covList, 99))
		covListPerc90 = str(np.percentile(covList, 90))
		covListPerc50 = str(np.percentile(covList, 50))
		BedGraphOpen = open(BedGraphFile, 'rb')
		BedGraphList = csv.reader(BedGraphOpen, delimiter = '\t')
		expfileOut = open(BedGraphNormFile, 'wb')
		writeExpfile = csv.writer(expfileOut, lineterminator='\n', delimiter = '\t')
		for row in BedGraphList:
			WindowSize = float(row[2]) - float(row[1])
			NormedReadDepth = (((float(row[3]))/ float(AllReads))*(1000000.0)*(1000.0/float(WindowSize)) - float(covListPerc50))/(float(covListPerc99P9) - float(covListPerc50))
			newrow = [ row[0], row[1], row[2] ]
			if NormedReadDepth > 0:
				writeExpfile.writerow(newrow + [round(NormedReadDepth,2)])
			else:
				writeExpfile.writerow(newrow + [0])
		BedGraphOpen.close()

def makeNormBigWig(inBam, GenomeIndex, BedGraphNormFile, BigWigNormFile):
	inBamBase = inBam.split(".bam")[0]
	#BedGraphNormFile=inBamBase+".slopped.normDepthNoise.bg"
	#BigWigNormFile=inBamBase+".slopped.normDepthNoise.bw"
	if os.path.isfile(BigWigNormFile):
		print "Output file " + BigWigNormFile + " already exists  \n"
	else:
		BigWigCommand="bedGraphToBigWig " + BedGraphNormFile + " "+ GenomeIndex + " " + BigWigNormFile
		subprocess.call(BigWigCommand, shell=True)


if __name__=="__main__":
	main()



