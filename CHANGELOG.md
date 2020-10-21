# Change Log for All files in Alignall repo 
## Doesn't need to be edited every commit, but keeping track of changes when they are made is helpful. 

#### alignall.py 
<date> | <change_in_detail>
------------------------------
180426 | Added exRNA options in steps and in main before setup(), as well as in setup().  Added exRNA code in align() under miRNA code.  exRNA won't use any hairpin files.  Will change for miRNA in future update.
180618 | Minor change: added library creation protocols in help notes
180924 | Major organization
182031 | Added reference to parseST() (function to make simple text files by parsing sample tracker)
190103 | Added version number, manually, to beginning, and option to print version and exit; also somewhat recently added random barcodes option for small RNA, which still needs a lot of work (generally just need to improve small RNA stuff)
190207 | Added import calls to all files that will eventually be used for SNP, ATAC, and ChIP analyses -> enough to change version to 12 (also going to be adding miRNA to this, as well as random barcode functionality)
190213 | Added parallel cutadapt prototype, added and cleaned up align() with trimming and cleaner miRNA parts, and updated version to 13 (atac part worked on myriam's, but haven't added accessibility and not sure if whole thing will work)
190321 | Changed the settings print out to print from variable which will also print to a "runLog.txt" file directly after setup
190424 | Added S_AUREUS genome detection, but only starpath, gtf, and refflat (and refflat might not even be working)
190430 | Added another S_Aureus genome, and a user-question to differentiate between S_AUREUS genomes
190509 | Added logfile output for kallisto run (outputs both stderr and stdout to log.txt)
190510 | Added hostname determination so we can establish the filesystem after knowing which machine we're on
190522 | Fixed hostname determination so it's simpler also added functionality if only want to do counts and several alignmentdirs exist
190603 | Created the new script (this one) to package up necessary components
190724 | Now works for single end rna-seq alignnment, continuing testing and fixing for atac, chip, snp, etc...
190726 | Added combineLanes, a function which combines files of the different lane, same sample
191014 | Python3-ized much of this script and added multiprocessing of combining lanes.
191205 | Been adding functionality for a lot of the callsnps option stuff, and added option countmirna for only doing mirna counts 
191206 | Changed strand options from 1,2,3 to 0,1,2
191211 | Standard for cutadapt: f we already trimmed, then files should be in fastqDir and program will detect if they're already trimmed so won't run cutadapt again
200127 | Polishing the variant calling pipeline (uses return values from multiprocessing pools to determine success/failure of command)
200420 | No longer going to work on smallRNA seq in this pipeline because docker image (genboree) is working on Linux machine for this
200512 | Adding dbsnp for mouse 
200522 | Added whole exonic file for mouse
200605 | Added an md5 generation step in case people wanted to submit data to GEO 
200707 | Added ability to select "OTHER" for genome, so that we can do wacky mouse-rat-human hybrids without adding to this script 



