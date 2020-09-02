import os, sys
import numpy as np

filelist = "filelist_30_v15_test.txt" 

with open(filelist) as f:
    lines = f.readlines()
    # print(data)
    for line in lines:
        if (line != lines[0]):
            words = line.split()
            value = words[1].find("LFI")
            # getting rid of path to file leaving only the filename intact
            words[1] = words[1][value:]
            pid = words[0]
            tod = words[1]#[4:-4]
            words[0] = tod
            words[1] = pid
            new_words = list()
            new_words.append(words[0])
            new_words.append(words[0][4:-4])
            new_words[1] = new_words[1][4:]
            new_words.append(words[1])
            new_words.append(words[2])
            new_words.append(words[3])
            new_words.append(words[4])
            print(f"{new_words}")
