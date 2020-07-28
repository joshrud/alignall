#import some stuff that's probably required 
import subprocess, os, shutil, string #command line calls, path showing, string things, deleting directories
import multiprocessing #for parallelizing stuff
from functools import partial
import sys #program control 
import re #for regex searching 
import gzip #used for parsing fastq.gz files
from datetime import datetime
from collections import OrderedDict


def runRmDups(alignmentDir, picard_jarfile, threads):
	now = datetime.now().strftime("%H:%M:%S")
	print("running rm dupes " + now)
	iList = []
	for i in os.listdir(alignmentDir):
		if "Aligned.sortedByCoord.duplicateRemoved.out.bam" not in os.listdir(alignmentDir+"/"+i) and \
		"Aligned.sortedByCoord.out.bam" in os.listdir(alignmentDir+"/"+i):
			iList.append(i)
		else:
			print("dupes already removed for: "+i)
	pool = multiprocessing.Pool(int(threads)) 
	func = partial(exRmDups, alignmentDir, picard_jarfile)
	pool.map(func, iList)
	pool.close()
	pool.join() 
	print("removed all duplicates. ")
	

def exRmDups(alignmentDir, picard_jarfile, j):
	now = datetime.now().strftime("%H:%M:%S")
	print(now + " marking duplicates in file: "+j)
	name = j.strip("out")
	subprocess.call("java -Xmx10g -jar "+picard_jarfile+" MarkDuplicates I="+alignmentDir+"/"+
		j+"/Aligned.sortedByCoord.out.bam O="+alignmentDir+"/"+j+
		"/Aligned.sortedByCoord.duplicateRemoved.out.bam TMP_DIR="+alignmentDir+"/"+j+
		" M="+alignmentDir+"/"+j+"/marked_dup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true > "+
		alignmentDir+"/"+j+"/picard.dupeRem.log.txt 2> "+alignmentDir+"/"+j+"/picard.dupeRem.err.txt", shell=True)


def runAssignGroups(alignmentDir, picard_jarfile, fastqDir, threads, readGroupOption):
	now = datetime.now().strftime("%H:%M:%S")
	print(now + " assigning groups...")

	iList = []
	for i in os.listdir(alignmentDir):
		name = i.strip("out")
		print("adding sample " + i)
		if "Aligned.sortedByCoord.duplicateRemoved.withReadGroups.out.bam" not in os.listdir(alignmentDir+"/"+i):
			iList.append(i)
	print(iList)
	pool = multiprocessing.Pool(int(threads)) # on 3 processors
	func = partial(exAssignGroups, alignmentDir, picard_jarfile)
	pool.map(func, iList)
	pool.close()
	pool.join() 


def exAssignGroups(alignmentDir, picard_jarfile, j):
	now = datetime.now().strftime("%H:%M:%S")
	print(now + " assigning groups in file: "+j)
	name = j.strip("out")
	rep = name.split('_')[-1] #when a file standard is present, this will be helpful, but it currently just adds giberish to the RG
	repLong = "Rep" + str(rep)	
	subprocess.call("java -Xmx10g -jar "+picard_jarfile+" AddOrReplaceReadGroups I="+alignmentDir+"/"+j+
		"/Aligned.sortedByCoord.duplicateRemoved.out.bam O="+alignmentDir+"/"+j+
		"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.out.bam TMP_DIR="+alignmentDir+"/"+j+
		" SORT_ORDER=coordinate RGID="+repLong+" RGLB="+repLong+" RGPL=illumina RGSM="+repLong+
		" RGPU=HiSeq4000 > "+alignmentDir+"/"+j+"/picard.assiGroup.log.txt 2> "+alignmentDir+"/"+j+"/picard.assiGroup.err.txt", shell=True)	


def runPreGATK(alignmentDir, threads, ref, indels, dbsnp, dbsnpHighCon, indel, GATKbase):
	now = datetime.now().strftime("%H:%M:%S")
	print(now + "running pre gatk...")
	iList = []
	for i in os.listdir(alignmentDir):
		if i.endswith("out"):
			iList.append(i)
	# exPreThreads = int(3 * int(threads))
	pool = multiprocessing.Pool(int(threads)) # run this many threads concurrently
	func = partial(exPreGATK, ref, indels, dbsnp, dbsnpHighCon, indel, alignmentDir, GATKbase)
	pool.map(func, iList)
	pool.close()
	pool.join()	


def exPreGATK(ref, indels, dbsnp, dbsnpHighCon, indel, alignmentDir, GATKbase, i):
	now = datetime.now().strftime("%H:%M:%S")
	print(now+" Running PreGATK on ... " + i + "\n")
	if not os.path.isfile(alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.out.bam"): #another check we need so we don't remake the .bai file
		subprocess.call("samtools index "+alignmentDir+"/"+i+
			"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.out.bam > "+alignmentDir+"/"+i+
			"/samtools.log.txt 2> "+alignmentDir+"/"+i+"/samtools.err.txt", shell=True)

	print("This step is slow, so checking if file exists and is non-empty")
	print("WARNING: Partial files will cause down stream errors")
	print("")
	outFile = alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam"
	if os.path.isfile(outFile):
		if os.stat(outFile).st_size > 0:
			print("file " + outFile + " exists, skipping")
		else: #-RF ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS removed in version 4.0 of gatk
			#https://gatkforums.broadinstitute.org/gatk/discussion/11567/splitncigarreads-unrecognized-options
			subprocess.call(GATKbase + " SplitNCigarReads -R " + ref + " -I "+alignmentDir+"/"+i+
				"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.out.bam -O "+alignmentDir+"/"+
				i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam  --tmp-dir "+
		 		alignmentDir+"/"+i+" > "+alignmentDir+
				"/"+i+"/split.log.txt 2> "+alignmentDir+"/"+i+"/split.err.txt", shell=True)
	else:
		subprocess.call(GATKbase + " SplitNCigarReads -R " + ref +
		 " -I "+alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.out.bam -O "+
		 alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam  --tmp-dir "+
		 alignmentDir+"/"+i+"> "+
		 alignmentDir+"/"+i+"/split.log.txt 2> "+alignmentDir+"/"+i+"/split.err.txt", shell=True)
		
	print("skipping BQSR steps...")
	# if indel == "yes": #let's just get rid of this...
	# 	subprocess.call(GATKbase + " RealignerTargetCreator -R "+ ref + " -I "+alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam  --known " + indels + " -O "+alignmentDir+"/"+i+"/forIndelRealigner.intervals > "+alignmentDir+"/"+i+"/realigner.log.txt 2> "+alignmentDir+"/"+i+"/realigner.err.txt", shell=True)
	# 	subprocess.call(GATKbase + " IndelRealigner -R " + ref + " -I "+alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam -known "+ indels + " -targetIntervals "+alignmentDir+"/"+i+"/forIndelRealigner.intervals -o "+alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.Realign.bam > "+alignmentDir+"/"+i+"/realigner.log.txt 2> "+alignmentDir+"/"+i+"/realigner.err.txt", shell=True)
	# 	subprocess.call(GATKbase + " BaseRecalibrator -R " + ref + " -nct 2 -I "+alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.Realign.bam  -knownSites " + dbsnpHighCon + " -knownSites " + indels + " -O "+alignmentDir+"/"+i+"/recal_data.table > "+alignmentDir+"/"+i+"/recalibrator.log.txt 2> "+alignmentDir+"/"+i+"/recalibrator.err.txt ", shell=True)
	# 	subprocess.call(GATKbase + " BaseRecalibrator -R " + ref + " -nct 2 -I "+alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.Realign.bam  -knownSites " + dbsnpHighCon + " -knownSites " + indels + " -BQSR "+alignmentDir+"/"+i+"/recal_data.table -O "+alignmentDir+"/"+i+"/post_recal_data.table > "+alignmentDir+"/"+i+"/recalibrator.log.txt 2> "+alignmentDir+"/"+i+"/recalibrator.err.txt", shell=True)
	# 	subprocess.call(GATKbase + " AnalyzeCovariates -R " + ref + " -before "+alignmentDir+"/"+i+"/recal_data.table -after "+alignmentDir+"/"+i+"/post_recal_data.table -plots "+alignmentDir+"/"+i+"/recalibration_plots.pdf > "+alignmentDir+"/"+i+"/recalibrator.log.txt 2> "+alignmentDir+"/"+i+"/recalibrator.err.txt", shell=True)
	# 	subprocess.call(GATKbase + " PrintReads -R " + ref + " -I "+alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.Realign.bam  -BQSR "+alignmentDir+"/"+i+"/recal_data.table -O "+alignmentDir+"/"+i+"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.Recal.bam > "+alignmentDir+"/"+i+"/recalibrator.log.txt 2> "+alignmentDir+"/"+i+"/recalibrator.err.txt", shell=True)
	# 	print(i + " completed")
	# 	print("")
	
	# elif indel == "no":
	# 	print("skipping indel realignment")
	# 	subprocess.call(GATKbase + " BaseRecalibrator -R " + ref + " -I "+alignmentDir+"/"+i+
	# 		"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam --known-sites " + 
	# 		dbsnpHighCon + " -O "+alignmentDir+"/"+i+"/recal_data.table > "+
	# 		alignmentDir+"/"+i+"/recalibrator.1.log.txt 2> "+alignmentDir+"/"+i+"/recalibrator.1.err.txt", shell=True)
	# 	subprocess.call(GATKbase + " BaseRecalibrator -R " + ref + " -I "+alignmentDir+"/"+i+
	# 		"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam --known-sites " +
	# 		 dbsnpHighCon + " -bqsr "+alignmentDir+"/"+i+"/recal_data.table -O "+
	# 		 alignmentDir+"/"+i+"/post_recal_data.table > "+alignmentDir+"/"+i+"/recalibrator.2.log.txt 2> "+
	# 		 alignmentDir+"/"+i+"/recalibrator.2.err.txt", shell=True)
	# 	subprocess.call(GATKbase + " AnalyzeCovariates -R " + ref + "  -before "+alignmentDir+"/"+i+
	# 		"/recal_data.table -after "+alignmentDir+"/"+i+"/post_recal_data.table -plots "+alignmentDir+"/"+i+
	# 		"/recalibration_plots.pdf > "+alignmentDir+"/"+i+"/recalibrator.3.log.txt 2> "+alignmentDir+"/"+i+
	# 		"/recalibrator.3.err.txt", shell=True)
	# 	subprocess.call(GATKbase + " PrintReads -R " + ref + " -I "+alignmentDir+"/"+i+
	# 		"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam -bqsr "+
	# 		alignmentDir+"/"+i+"/recal_data.table -O "+alignmentDir+"/"+i+
	# 		"/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.Recal.bam > "+
	# 		alignmentDir+"/"+i+"/recalibrator.4.log.txt 2> "+alignmentDir+"/"+i+"/recalibrator.4.err.txt", shell=True)
	# 	print(i + " completed")
	# 	print("")


def runMergeAndHaplotypeCaller(exonBed, threads, ref, dbsnp, alignmentDir, GATKbase):
	
	hapDir = alignmentDir+"/Merged_HaplotypeCalling"
	alignmentfolders = os.listdir(alignmentDir)
	if "Merged_HaplotypeCalling" in alignmentfolders:
		alignmentfolders = alignmentfolders.remove("Merged_HaplotypeCalling")
		user_continue = input("Merged_HaplotypeCalling/ has already been created, \
			\ncontinuing will overwrite current files. Continue? (Y/N): ")
		if user_continue.lower() != "y":
			sys.exit()
		else:
			shutil.rmtree(hapDir)
			os.mkdir(hapDir)
	else: 
		os.mkdir(hapDir)

	splitfolders = [f for f in alignmentfolders if "Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam" in os.listdir(alignmentDir+"/"+f)]
	if len(splitfolders) < len(alignmentfolders):
		print("something went wrong in run pre gatk, rerun script.")
		sys.exit()

	# reps = [] #what if there aren't any replicates
	# repNumber = []
	# for i in os.listdir(alignmentDir):
	# 	if i.endswith("out"):
	# 		name=i.strip('out')
	# 		rep=name.split('_')[1]
	# 		repLong = "Rep" + str(rep)
	# 		if rep not in repNumber:
	# 			repNumber.append(rep)

	masterList = []
	splitCount = 0			
	bedLines = []
	with open(exonBed,'r') as text_file:
		for row in text_file:
			bedLines.append(row)
			splitCount = splitCount + 1
	
	split = round(splitCount / 50)
	outer_count = 1
	line_count = 0
	sorting = True
	while sorting:
		count = 0
		increment = (outer_count-1) * split
		left = len(bedLines) - increment
		file_name = "subset_exons_" + str(outer_count * split) + ".list" #gatk - style somehow works, when bed doesn't (bed = (<chr>\t<start>\t<stop>); list=(<chr>:<start>-<stop>)) 
		hold_new_lines = []
		if left < split:
			while count < left:
				hold_new_lines.append(bedLines[line_count])
				count += 1
				line_count += 1
			sorting = False
		else:
			while count < split:
				hold_new_lines.append(bedLines[line_count])
				count += 1
				line_count += 1
		outer_count += 1
		with open(alignmentDir+"/Merged_HaplotypeCalling/"+file_name,'w') as next_file:
			for row in hold_new_lines:
				next_file.write(row)
	
	bedList = []
	for l in os.listdir(alignmentDir+"/Merged_HaplotypeCalling"):
		if "subset_exons" in l and "raw" not in l and "filtered" not in l: # and l.find("subset") > 1
			bedList.append(l)
	
	outBamList = alignmentDir+"/Merged_HaplotypeCalling/bamsIn.list"
	
	bamList = []
	for k in os.listdir(alignmentDir):
		if k.find("out") > 1:
			bam = alignmentDir +"/"+k + "/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.bam"
			bamList.append(bam)

	with open(outBamList,'w') as bamListOutFile:
		for n in bamList:
			bamListOutFile.write("%s\n" % n)
	
	#miniList.append(arg)
	#masterList.append(miniList)
	print("running haplotype caller...")
	pool = multiprocessing.Pool(int(threads)) 
	func = partial(exHaplotypeCaller, ref, dbsnp, alignmentDir, GATKbase, outBamList)
	print(bedList)
	pool.map(func, bedList)
	pool.close()
	pool.join()	

	#now merge using samtools 
	# subprocess.call("samtools merge ")


def exHaplotypeCaller(ref, dbsnp, alignmentDir, GATKbase, bamListInput, bedListInput):
	print("working on ", bedListInput, " while in directory: ", os.path.abspath("."))

	outBase = bedListInput.rstrip('.list')
	#rep = str(input[0])
	#arg = str(input[1])
	#print( "Running merge of all Rep" + rep + " ... "
	#subprocess.call(GATKbase + " printReads -R " + ref  + " " + arg + " --read_filter MappingQualityZero -o Alignments/Merge_Rep" + str(rep) + "/Aligned.sortedByCoord.duplicateRemoved.withReadGroups.Split.Realign.Recal.merge.bam > Alignments/Merge_Rep" + str(rep) + "/merge.log.txt 2> Alignments/Merge_Rep" + str(rep) + "/merge.err.txt", shell=True)
	# subprocess.call("touch "+alignmentDir+"/Merged_HaplotypeCalling/output." + outBase + ".raw.log.txt", shell=True)
	# subprocess.call("touch "+alignmentDir+"/Merged_HaplotypeCalling/output." + outBase + ".raw.err.txt", shell=True)
	# subprocess.call("touch "+alignmentDir+"/Merged_HaplotypeCalling/output." + outBase + ".filtered.log.txt", shell=True)
	# subprocess.call("touch "+alignmentDir+"/Merged_HaplotypeCalling/output." + outBase + ".filtered.err.txt", shell=True)
	if not os.path.isfile(alignmentDir+"/Merged_HaplotypeCalling/output."+outBase+".raw.vcf"):
		print("Calling SNPs in exonic region: "+outBase+" using merge of all samples...") #no longer an option: --stand_emit_conf 20.0
		subprocess.call(GATKbase + " HaplotypeCaller -L " + alignmentDir + "/Merged_HaplotypeCalling/" + 
			bedListInput + " -R " + ref + " --dbsnp " + dbsnp + " -I " + bamListInput + 
			" --dont-use-soft-clipped-bases -stand-call-conf 20.0 -O "+alignmentDir+"/Merged_HaplotypeCalling/output." + 
			outBase + ".raw.vcf > "+alignmentDir+"/Merged_HaplotypeCalling/output." + outBase + ".raw.log.txt 2> "+
			alignmentDir+"/Merged_HaplotypeCalling/output." + outBase + ".raw.err.txt", shell=True) 
	print("Filtering SNPs on merge of all Rep" + outBase + " ... ")
	subprocess.call(GATKbase + " VariantFiltration -R " + ref + " -L " + alignmentDir + "/Merged_HaplotypeCalling/" +
		bedListInput + " -V "+alignmentDir+"/Merged_HaplotypeCalling/output." + outBase + 
		".raw.vcf -window 35 -cluster 3 --filter-name FS -filter \"FS > 30.0\" --filter-name QD -filter \"QD < 2.0\" -O "+
		alignmentDir+"/Merged_HaplotypeCalling/output." + outBase + ".filtered.vcf > "+alignmentDir+
		"/Merged_HaplotypeCalling/output." + outBase + ".filtered.log.txt 2> "+alignmentDir+"/Merged_HaplotypeCalling/output." + 
		outBase + ".filtered.err.txt", shell=True)


