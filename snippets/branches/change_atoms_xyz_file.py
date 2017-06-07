#usr/bin/env python

import os
import argparse
import numpy as np
atomcount = 53


parser=argparse.ArgumentParser(description='This script it for converting the methyl group to TFM')

parser.add_argument('molfile_name', help='required input; /your/molfile/location/and/name.xyz; no default')

args = parser.parse_args()
molfile = args.molfile_name #This line sets the file that is being converted to the variable molfile

def switch_xyz(input1):
    f = open('%s.xyz' %input1, 'a')
    f.write('53 \n')
    f.write('\n')
    atoms=np.genfromtxt(fname=input1,dtype="S8", skip_header=2)
    atoms[51,0] = "F"
    atoms[52,0] = "F"
    atoms[50,0] = "F" 
    for lines in atoms:
        cleanline = "     ".join(lines)
	f = open('%s.xyz' %input1, 'a')
	f.write('%s \n ' %cleanline) 

switch_xyz(molfile)
