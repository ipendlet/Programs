#usr/bin/env python

#writing this code to python standard.  This is designed to function so far with all molfile types.
#requirements to get this working.
#1. Needs to have the basis for each of the molecules in question labeled as element.basis (ex. H.basis for hydrogen)
#2. Basis molfiles need to be saved to a folder in the directory of the python script and input molfile named basis_sets
#3. The molfile needs to be specified in the argument of the python script
#4. The name of the output molfile will be nameofinput.xyz.inp
#5. If ECP is inluded in the basis set add -ECP to the python command
#5.1 If you are not sure whether you basis needs ECP, then I would recommend double checking this before going futher!!!!!!
#6. Test on one molecule molfile before running on all!!!!
#7. If you do not follow these steps correctly do not be surprised if gamess hates you, breaks, or the world comes to a screeching halt.
#8 BACKUP ALL FILES BEFORE RUNNING THIS SCRIPT! (Hopefully you read all the way through to this.. if not.. sucks to be you.)
import os
import argparse
import numpy as np
import openbabel

parser=argparse.ArgumentParser(description='This script is intended for converting xyz molfiles into useable molfiles for gamess input. The basis set and header can be specified in the molfile "header"')

parser.add_argument('molfile_name', help='required input; /your/molfile/location/and/name.xyz; no default')
parser.add_argument('-ECP', type=int, default=0, help='sets whether the program will incorporate ECP information from the basis sets. Default is 0 (off), set to 1 if needed (on)')
parser.add_argument('-Basis', type=int, default=0, help='sets whether the program will incorporate manually creates basis sets. Default is 0 (off), meaning that the user needs to specify a basis through the inline option in gamess.  The option should be set to 1  (on) if manual entry is desired')

args = parser.parse_args()
molfile = args.molfile_name #This line sets the file that is being converted to the variable molfile
ECP = args.ECP
Basis = args.Basis

########CODE BEGINNING#########
#Just a check to make sure nothing died at the start...
d={}
def warmup():
    if os.path.exists(molfile):
        print "File Found"
    else:
        print 'File not located'

# The following commands are used in combination with openbabel to convert the file to the gamess input for further work
def convertxyztogamess(input1): 
    conv = openbabel.OBConversion()
    conv.OpenInAndOutFiles(input1, "%s_converted" %molfile)
    conv.SetInAndOutFormats("xyz", "inp")
    conv.Convert()
    confile = open("%s_converted" %molfile)
    lines = confile.read()
    lines = lines[:-8]
    f = open("%s_converted" %molfile, 'w')
    f.write("%s" %lines)
    return "%s_converted" %molfile

def Autobasis_finalsteps_configuration(input1): #This section is responsible for add the basis information following the respective row for each of the elements
    confile = convertxyztogamess(input1)
    rowgrab = np.genfromtxt(fname=confile, dtype="S8", skip_header=5)
    x = 0
    header = open("header")
    headerread = header.read()
    f = open("%s.inp" %confile, 'w')
    f.write("%s" %headerread)
    for lines in rowgrab:
	cleanreadline = "      ".join(lines)
        f.write("%s \n" %cleanreadline)
#	with open("%s.inp" %confile, "a") as myfile:
#	    myfile.write("%s \n" %cleanreadline)
        x+=1
    f.write(" $END \n")

def basis_set_configuration(input1): #This section is responsible for add the basis information following the respective row for each of the elements
    confile = convertxyztogamess(input1)
    rowgrab = np.genfromtxt(fname=confile, dtype="S8", skip_header=5)
    x = 0
    header = open("header")
    headerread = header.read()
    f = open("%s.inp" %confile, 'w')
    f.write("%s" %headerread)
    for lines in rowgrab:
	cleanreadline = "      ".join(lines)
	elementname = rowgrab[x,0]
	f = open('basis_sets/%s.basis' %elementname)
        fread = f.read()
	#print rowprint
	with open("%s.inp" %confile, "a") as myfile:
	    myfile.write("%s \n" %cleanreadline)
	    myfile.write("%s \n" %fread)
        x+=1
    with open("%s.inp" %confile, "a") as myfile:
	myfile.write(" $END \n")

def ECP_configuration(input1):
    confile = convertxyztogamess(input1)
    rowgrab = np.genfromtxt(fname=confile, dtype="S8", skip_header=5)
    x = 0
    header = open("header")
    headerread = header.read()
    f = open("%s.inp" %confile, 'w')
    f.write("%s" %headerread)
    for lines in rowgrab:
	cleanreadline = "      ".join(lines)
	elementname = rowgrab[x,0]
	f = open('basis_sets/%s.basis' %elementname)
        fread = f.read()
	#print rowprint
	with open("%s.inp" %confile, "a") as myfile:
	    myfile.write("%s \n" %cleanreadline)
	    myfile.write("%s \n" %fread)
        x+=1
    with open("%s.inp" %confile, "a") as myfile:
	myfile.write(" $END \n")
    x = 0
    with open("%s.inp" %confile, "a") as myfile:
	myfile.write(" $ECP \n")
        for lines in rowgrab:
	    elementname = rowgrab[x,0]
            #try:
	     #    (d['%s.elementcounter' %elementname])
	#	 print '%s.fail' %elementname
	#    except KeyError:
	 #        d['%s.elementcounter' %elementname]="hello"
	#	 print '%s.success' %elementname
	    try:
	         with open('basis_sets/%s.basis.ecp' %elementname):
                    try:
	               (d['%s.elementcounter' %elementname])
		       #print '%s.fail' %elementname
		       myfile.write("%s NONE \n" %elementname)
	            except KeyError:
	               d['%s.elementcounter' %elementname]="hello"
	 	       #print '%s.success' %elementname
	               f = open("basis_sets/%s.basis.ecp" %elementname)
		       fread = f.read()
		       myfile.write("%s" %fread)
	    except IOError:
		    myfile.write("%s NONE \n" %elementname)
	    x+=1	
	myfile.write(" $END \n")


def RunAllFinalFunction(input1): 
    if Basis == 0:
        Autobasis_finalsteps_configuration(input1)
    else:
        if ECP == 1:
	     ECP_configuration(input1)
        else:
	    print "ECP was not included!! Make sure you don't need this. See help for how to toggle!"
            basis_set_configuration(input1)

RunAllFinalFunction(molfile)
