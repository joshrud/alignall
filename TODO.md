# The To do list for the alignall repo
## If you think of any ideas, please add to the bottom of the list and each idea will be triaged by importance 


### High Priority Todos 
- standardize how inputs are given: make really simple to start using, but can supply a .alignconfig file or something
- make TEST script and sample files with ~1000 reads each, with maybe 8 samples in each test experiment, including single cell test also; find random dataset online for this
- check adapter before trimming
- should print some sort of log file that can be checked by alignall if run again on same files and used to set checkpoints	


### Medium Priority Todos
- standardize how reference files are downloaded and isntalled (maybe make an install option for this script to get all files and put them into a reference directory (e.g. --get MOUSE ; --get HUMAN ))
- remove all of the miRNA and smallRNA parts of this script or link to docker image (maybe just in 'help' say to run docker image)
- integrate alignall_sc2 and test
- accumulate getReadLengths.py and cutadaptAll.py and integrate into this script (this should just be referencing it)


### Low Priority Todos
- test and get atac/chip working (low priority because I don't need to do this anymore)
- clean up code so that it doesn't go past 80 characters in line length


