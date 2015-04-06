#!/usr/bin/env python

'''
Author: Ian Tibbetts
Co-authors: Prof. Elizabeth McGrath
Last Edited: 4/5/2015
Colby College Astrophysics Research
'''

import sys
import os
from collections import OrderedDict

def mergeSummaries(summary1, summary2, delim=","):
    '''
    merge the summary information from summary 2 into summary 1
    '''
    merged = OrderedDict()
    gnameCol = -1
    avalCol = -1
    for col, header in enumerate(summary1[0].split(delim)):
        if header == "galaxyID":
            gnameCol = col
        elif header == "timeStep":
            avalCol = col
    if gnameCol == -1 or avalCol == -1:
        print("summary 1 does not have expected headers")
        exit()
    for row in summary1[2:]:
        fields = row.split(delim) + ["0", "0", "0"]
        gname = fields[gnameCol]
        aval = float("0."+fields[avalCol])
        key = "%s%.6f"%(gname,aval)
        if key in merged:
            merged[key].append(fields)
        else:
            merged[key] = [fields]
            
    sfrCol = -1
    massCol = -1
    gnameCol = -1
    avalCol = -1
    for col, header in enumerate(summary2[0].split(delim)):
        if header == "sim":
            gnameCol = col
        elif header == "aexp":
            avalCol = col
        elif header == "mass_sfr":
            sfrCol = col
        elif header == "mass_stars":
            massCol = col
    if sfrCol == -1 or massCol == -1:
        print("summary 2 does not have expected headers")
        exit()
    for row in summary2[1:]:
        fields = row.split(delim)
        gname = fields[gnameCol]
        aval = float(fields[avalCol])
        sfr = fields[sfrCol]
        mass = fields[massCol]
        ssfr = str(float(sfr)*(10**9)/float(mass))
        key = "%s%.6f"%(gname,aval)
        try:
            for dup in merged[key]:
                dup[-3:] = [sfr, ssfr, mass]
            print("found key "+key)
        except KeyError:
            print("no key " + key)
            
    summary1[0] = delim.join([summary1[0], "mass_sfr", "mass_ssfr", "mass_stars"])
    mergedLines = summary1[:2]
    for key in merged:
        for dup in merged[key]:
            mergedLines.append(delim.join(dup))
    return mergedLines
        
if __name__ == "__main__":

    # define the command line interface with simUtility.py
    usage = ("USAGE: python mergeSummaries.py <summary1> <summary2>")

    # must have exactly one positional command line argument
    if (len(sys.argv) != 3 
        or not os.path.isfile(sys.argv[1])
        or not os.path.isfile(sys.argv[2])):
        print(usage)
        exit()
        
    # read summary files
    with open(sys.argv[1], 'r') as sumFile1:
        sum1 = [line.rstrip("\n") for line in sumFile1.readlines()]
    with open(sys.argv[2], 'r') as sumFile2:
        sum2 = [line.rstrip("\n") for line in sumFile2.readlines()]
    
    # attempt to merge the summary files
    mergeSumLines = mergeSummaries(sum1, sum2)
    
    # write the merged summary back to summary 1
    fnSplit = sys.argv[1].split(".")
    with open(".".join(fnSplit[:-1])+"_merged."+fnSplit[-1], 'w') as sumMerge:
        sumMerge.write("\n".join(mergeSumLines))