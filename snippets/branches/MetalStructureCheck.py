#usr/bin/env python
##module list###
import sys
import os
import argparse
import numpy as np
import pybel as pb
#import openbabel
##smallchange
#### variables ####
#operating directory#
directory= '/export/zimmerman/ipendlet/4-CobaltChemistry/2_Combinations/4-OptOctahedral_H2/zzPreFinalXYZSTRUCTURES_needconfirMzzz/'
#total bonded thresh
tb = 6

### File line arguments ####
#parser=argparse.ArgumentParser(description='This script is intended for converting xyz molfiles into useable molfiles for gamess input. The basis set and header can be specified in the molfile "header"')

#parser.add_argument('molfile_name', help='required input; /your/molfile/location/and/name.xyz; no default')
#parser.add_argument('Length_name')
#parser.add_argument('Angle_name')

##demoline## parser.add_argument('-ECP', type=int, default=0, help='sets whether the pro')

#args = parser.parse_args()
#Sets the inputfiles to the variable molfile for each run#
#molfile = args.molfile_name 
#lengthfile = args.Length_name

########File Handling - Program operates on all available XYZ files in directory#########
def ddcheck(molfile):
    if os.path.exists(molfile):
        print "File Found"
    else:
        print 'File not located'
### 

##output bonded atoms ###
def OBread(mol1):
    for atom in mol1:
        if atom.atomicnum == 27 :
            obatom = atom.OBAtom
            for bonded in pb.ob.OBAtomAtomIter(obatom):
                print bonded.GetAtomicNum()

## output number of connected atoms ##
def OBconnect(mol1, file):
    for atom in mol1:
        if atom.atomicnum == 27:
            obatom = atom.OBAtom
#            print len(list(pb.ob.OBAtomAtomIter(obatom)))
            if len(list(pb.ob.OBAtomAtomIter(obatom))) < 6:
                file2 = file.strip("'")
                print file2
                #, len(list(pb.ob.OBAtomAtomIter(obatom)))
            else:
                pass

#handles arrays based on input from XXX 
def CalcBondLengths(array, X, Y):
    xyza = np.genfromtxt(array, skip_header=2, usecols=(1,2,3))
    a = xyza[X,:]
    print 'a'
    print a
    b = xyza[Y,:]
    print 'b'
    print b
    dist = np.linalg.norm(a-b)
    return dist

#Takes input from the user defined length array file and analyzes each of the bonds
def ProcessNeededLengths(xyzfile, lnag):
    Array = np.genfromtxt(lnag)
    mA = Array.shape
    m = np.array(mA)
    n = m[0]
    i=0
    while i<n:
        atom1 = Array[i,0]
        atom2 = Array[i,1]
        X = atom1-1
        Y = atom2-1
        Solve = CalcBondLengths(xyzfile, X, Y)
        i=i+1
        f = open('out.txt', 'ab')
        print >> f, Solve, atom1, atom2
        f.close()
        #print 

for file in os.listdir(directory):
#    mol1 = next(pb.readfile("xyz", molfile))
    if file.endswith(".xyz"):
        mol1 = next(pb.readfile("xyz", file))
#        OBread(mol1)
        OBconnect(mol1, file)
#        print(os.path.join("/mydir", file))
#        print(file)

#OBread(mol1)
