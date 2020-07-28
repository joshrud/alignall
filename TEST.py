# script designed to test variations of alignall.py

#tests: 
#1. rna seq
#2. variant calling
#3. chip with input
#5. atac (behaves the same as chip without input)
#6. single cell
#7. single cell variant calling (demuxlet & freemuxlet)
#8. rna seq -> variant calling -> single cell -> demux (if time)
#9. adapter clipping
#10. kallisto (transcriptome alignment)
#11. picard and qc 

#features:
# - timing (need a @time decorator function to time each test)
# - logs (write a central log file for general processes and an individual one for each process run)
# - more tests (there's probably lots of edge cases that won't get tested by above tests)



REFDIR = "" #root directory of all above test files


