#
# 
#

# include some funny functions
import os
from os import listdir
from os.path import isfile, join

# get the list of files
mypath = 'shortinfo/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

i = 0;
with open('shortinfo_final.txt', 'w') as outfile:
    for fname in onlyfiles:
        i = i + 1
        if i==1:
            with open(mypath+fname) as infile:
                for line in infile:
                    outfile.write(line)
                    
        else:
            with open(mypath+fname) as infile:
                j = 0
                for line in infile:
                    j = j + 1
                    if j >= 2:
                        outfile.write(line)
 
 
