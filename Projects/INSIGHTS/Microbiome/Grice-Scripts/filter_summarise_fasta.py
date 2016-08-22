#!/usr/bin/python

# Arguments following script should be as follows: fasta, min, max

# Import libraries
import sys
import os.path

# Read in command line arguments
fasta = str(sys.argv[1])
outName = os.path.splitext(fasta)[0]
minLen = int(sys.argv[2])
maxLen = int(sys.argv[3])

# Make a dictionary containing all of the lengths and their counts
readLengths = {}

# Open the text file for output
keep = open(outName + "_size_filtered.fna", "w")
discard = open(outName + "_size_discarded.fna", "w")

# Read in file one line at a time
with open(fasta) as infile:
    for line in infile:
        # Make sure you aren't counting the identifier line
        if line.startswith('>'):
        	# store header in case sequence gets past filtering
        	header = line
        	continue
    	else:
    		# Get length
        	length = len(line)
    		
			# Add 1 to length to dictionary, if doesn't exist, create it
	        if length in readLengths:
	        	readLengths[length] += 1
	    	else:
	    		readLengths[length] = 1

    		# size filter
    		if length >= minLen and length <= maxLen:
    			# Write sequence to keep file
    			keep.write(header)
    			keep.write(line)
	    	else:
    			# Write sequences to discard file
    			discard.write(header)
    			discard.write(line)

# Close text files
keep.close()
discard.close()

# print out the sorted lengths        
with open(outName + "_length_summary.csv", "w") as length_summary:
	for key in sorted(readLengths):
		line = "%s,%s\n" % (key,readLengths[key])
		length_summary.write(line)
