# importing libraries
import os
import sys
from itertools import islice

# defining location of parent folder
BASE_DIRECTORY = 'DirectoryFolder'
output_file = open('output.txt', 'w')
output = {}
file_list = []

# read all txt files in the folder and store it in file_list
for (dirpath, dirnames, filenames) in os.walk(BASE_DIRECTORY):
    for f in filenames:
        if 'txt' in str(f):
            e = os.path.join(str(dirpath), str(f))
            file_list.append(e)

# find all lines start with key string and append it to output. 
for f in file_list:
    output[f] = []
    with open(f, "r") as ifile:
        for line in ifile:
            if line.startswith("Keyword"):            
                output[f].append(line)
                output[f].append(next(ifile, '').strip())


# if need, we can also output several lines around the key string line:
for f in file_list:
    output[f] = []
    with open(f, "r") as ifile:
        linesBefore = list()
        for line in ifile:
            linesBefore.append(line.rstrip())
            if len(linesBefore) > 4: #Adding up to 4 lines
                linesBefore.pop(0)
            if "Description" in line:
                #if len(linesBefore) == 4: # if there are at least 3 lines before the match
                for i in range(3):
                    output[f].append(linesBefore[i])
                output[f].append(("".join(line.rstrip())))
                #print ("".join(islice(f,3)))
                output[f].append((next(ifile, '').rstrip()))

# now output the results to a txt file. 
tabs = []
for tab in output:
    tabs.append(tab)

tabs.sort()
for tab in tabs:
    output_file.write(tab + '\n')
    output_file.write('\n')
    for row in output[tab]:
        output_file.write(row + '')
    output_file.write('\n')
    output_file.write('----------------------------------------------------------\n')

