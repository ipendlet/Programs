#usr/bin/env python

import sys
import argparse
import numpy as np
import openbabel

parser=argparse.ArgumentParser(description='This script is intended for converting xyz molfiles into useable molfiles for gamess input. The basis set and header can be specified in the molfile "header"')

parser.add_argument('molfile_name', help='required input; /your/molfile/location/and/name.xyz; no default')
parser.add_argument('Length_name')
parser.add_argument('Angle_name')

#parser.add_argument('-ECP', type=int, default=0, help='sets whether the program will incorporate ECP information from the basis sets. Default is 0 (off), set to 1 if needed (on)')
#parser.add_argument('-Basis', type=int, default=0, help='sets whether the program will incorporate manually creates basis sets. Default is 0 (off), meaning that the user needs to specify a basis through the inline option in gamess.  The option should be set to 1  (on) if manual entry is desired')

args = parser.parse_args()
molfile = args.molfile_name #This line sets the file that is being converted to the variable molfile
lengthfile = args.Length_name
anglefile = args.Angle_name



########CODE BEGINNING#########
#Just a check to make sure nothing died at the start...
d={}
def ddcheck(molfile):
    if os.path.exists(molfile):
        print "File Found"
    else:
        print 'File not located'

#handles arrays based on input from XXX 
def CalcBondLengths(array, X, Y):
    xyza = np.genfromtxt(array, skip_header=2, usecols=(1,2,3))
    a = xyza[X,:]
    b = xyza[Y,:]
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
        #print i
        
#handles arrays based on input from XXX 
def CalcBondAngles(array, X, Y, Z):
    xyza = np.genfromtxt(array, skip_header=2, usecols=(1,2,3))
    a = xyza[X,:]
    b = xyza[Y,:]
    c = xyza[Z,:]
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)

#Takes input from the user defined length array file and analyzes each of the bonds
def ProcessNeededAngles(xyzfile, lnag):
    Array = np.genfromtxt(lnag)
    mA = Array.shape
    m = np.array(mA)
    n = m[0]
    i=0
    while i<n:
        atom1 = Array[i,0]
        atom2 = Array[i,1]
        atom3 = Array[i,2]
        X = atom1-1
        Y = atom2-1
        Z = atom3-1
        Solve = CalcBondAngles(xyzfile, X, Y, Z)
        i=i+1
        f = open('out.txt', 'ab')
        print >> f, Solve, atom1, atom2, atom3
        f.close()
        #print i
        


#sorts out which of user requests are for lengths and which are for angles (i.e. just trying to get how many are arrayed)
#def Sort(xyzfile, lnag, angfil) 
#def Sort(xyzfile, lnag)
    

#Final Calls
#BondLengths(molfile, lengthfile, anglefile)
#f = open('out.txt', 'w')
#print >> f, ProcessNeededLengths(molfile, lengthfile)
#f.close()

f = open('out.txt', 'ab')
print >> f, molfile
f.close()
ProcessNeededLengths(molfile, lengthfile)
ProcessNeededAngles(molfile, anglefile)
