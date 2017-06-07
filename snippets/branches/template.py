#usr/bin/env python


##### NOT up to python standard.  All lines that need to be fixed are identified by:
# MAC SYSTEM CALL 

###prerequisite information for importing specific file, no modification of file
import os
import argparse #separates information from linux command
import re
import numpy as np
#import pybel
import openbabel
from StringIO import StringIO
parser = argparse.ArgumentParser(description='Converts inputfile.xyz to ElementaryStep_moleculeX.cdxml where X represents the number of times this script executed')
parser.add_argument('file_name', help='required input; /your/file/location/and/name.xyz; no default')
parser.add_argument('-pages_num', type=int, default=4, help='increases the pages occupied by output molecules; increase for more space: default= 4pages')
args = parser.parse_args()
files = args.file_name #files will remain an argument throughout this script!!
pages = args.pages_num

#variables
ytotalunits = 0 #global var. written in pagecal()
ytotal = 0 #global var. written in pagecal()
scratchvar = 1 #global var. written in scratch()
molvar = 0 #global var. written in moleculeinfo() 
ytop = 0 #global var. written in pagespace()
ybottom = 0 #global var. written in pagespace()
firstrun = 0 #global var. written in moleculeinfo()
molfoldvar = 1 #global var. written in moeculeinfo()
molcount = 1 #global var. written in xyzsplitcount() 
atct = 0 #global var. written in atomcount()
Xpagescale = int(10) #global variable written in xadjustment()

###beginning of code
def filessearchconfirm():
    if os.path.exists(files): # MAC SYSTEM CALL 
        print "File Found"
    else:
        print 'Files not located'

#creates scratch files for processing without overright
def scratch():
    global scratchvar
    x = int(scratchvar)
    if not os.path.exists('scratch/'): # MAC SYSTEM CALL 
        os.makedirs("scratch/") # MAC SYSTEM CALL 
    file = open("scratch/scratch_%d.txt" % x, 'w+')
    scratchvar +=1
    return "scratch/scratch_%d.txt" % x 

#creates directory and files for final molecule output
def moleculeinfo():
    global molvar
    global firstrun
    global molfoldvar
    molvar +=1
    m = int(molvar)
    f = int(firstrun)
    while firstrun == 0:
        if os.path.exists('molecule_%d/' % molfoldvar): # MAC SYSTEM CALL 
            molfoldvar +=1
        else:
            os.makedirs("molecule_%d/" % molfoldvar)
            firstrun +=1
    file = open("molecule_%d/molecule_%d.txt" % (molfoldvar, m), 'w+')
    return "molecule_%d/molecule_%d.txt" % (molfoldvar, m) 

#read atom header
def atomheader(inputfile):
    f=open(inputfile)
    lines = f.readlines()
    atoms = lines[0] 
    return atoms

#reads number of atoms
def atomcount(inputfile):
    f=open(inputfile)
    global atct
    lines = f.readlines()
    atoms = lines[0]
    a = atoms.strip() 
    atoms = int(a)
    atct = atoms
    return atoms


#adds energy to the comment line w/o disturbing other text (if it is not already present)  Mostly error prevention
def addenergy(inputfile):
    temp = scratch() 
    f = open(temp, 'w+')
    array = open(inputfile)
    top = atomheader(inputfile)
    for line in array.readlines():
        if re.search(top, line):
            x = [value for value in line.split()]
            g =('[%s]' % '  '.join(map(str, x)))
            y = g[1 : -1]
            print>>f, y
        elif re.search(r"[A-Za-z]", line):
            x = [value for value in line.split()]
            g =('[%s]' % '  '.join(map(str, x)))
            y = g[1 : -1]
            print>>f, y
        else:
            x = [value for value in line.split()]
            g =('[%s]' % ' '.join(map(str, x)))
            y = g[1 : -1]
            print>>f, 'Energy: %s' % y
    return temp

#xy array template, cuts the text in line indicated within input file from argument    
def xyzsplitcount(inputfile): #input file is coming from args labeled as 'files'
    if not os.path.exists('numpyexport.npy'): # MAC SYSTEM CALL  
        array = open(inputfile)
        global molcount
        for lines in array.readlines():
                if 'Energy' in lines: # change countvar at top to change count variable
                    molcount +=1
                else:
                    pass
        ep = int(molcount-1)
        xytable = np.zeros(shape=(4, ep)) #forms array template [molecule, xcoordinate, ycoordinate, energy of molecule,]
        xytable[0] = range(1, molcount) 
        fname = 'numpyexport'
        np.save(fname, xytable)
        return fname+'.npy' #output saved molecule array as "numpyexport.npy"
    else: 
        return 'numpyexport.npy' #output saved molecule array as "numpyexport.npy"

#calibrates energy to page length depending on scale
def pagecal():
    global ytotalunits
    global ytotal
    ytotal = ytop-ybottom
    totalunits = 24*pages
    ytotalunits = float(ytotal/totalunits)
    
#reads in y min and max for page setting#
def pagespace(): 
    global molfoldvar
    global molvar
    atomcount(files)
    global atct
    global ytop
    global ybottom 
    yminfile = 2 
    if molvar > yminfile:
        y1 = "molecule_%d/molecule_%d.txt" % (molfoldvar, yminfile) 
        y2 = "molecule_%d/molecule_%d.txt"  % (molfoldvar, molvar)
        yfind1 = np.genfromtxt(y1, skiprows=1, delimiter=",", usecols=(1),  skip_footer=atct)
        yfind2 = np.genfromtxt(y2, skiprows=1, delimiter=",", usecols=(1),  skip_footer=atct)
        ytop = yfind2
        ybottom = yfind1
    else:
        y2 = "molecule_%d/molecule_%d.txt" % (molfoldvar, 1) 
        yfind2 = np.genfromtxt(y2, skiprows=1, delimiter=", ", usecols=(1),  skip_footer=atct)
        ytop = abs(yfind2)
    pagecal()
    
#creates scratch/molecule, output array n[scratch, file#]
def xyzsplit(inputfile):
    splitfilenames = open('splitfiles.txt', "w+")
    f = open(inputfile)
    atoms = atomheader(inputfile) 
    for line in f.readlines():
        if line.startswith(atoms):        
            x = [value for value in line.split()]
            g =('[%s]' % ', '.join(map(str, x)))
            y = g[1 : -1]
            mkfile = "%s" %(moleculeinfo())
            fp = open(mkfile, "w+")
            print>> fp, y
            print>> splitfilenames,  mkfile
        else: 
            x = [value for value in line.split()]
            g =('[%s]' % ', '.join(map(str, x)))
            y = g[1 : -1]
            print>> fp, y
    return 'splitfiles.txt'

#export reformatted array
def gbarray(inputfile):
    temp = scratch()
    fh = open(inputfile)
    f = open(temp, "w+")
    for lines in fh.readlines():
         x = [value for value in lines.split()]
         g =('[%s]' % ' '.join(map(str, x))) 
         y = g[1 : -1]
         print>> f, y
    f = open(temp, "r")
    return f.read()

#reads xyz file header and returns header only
def Format(inputfile):
    headerfile = open(inputfile)
    count = 0
    temp = scratch()
    writefile = open(temp, "w+")
    for lines in headerfile.readlines(): #formatting loop
        if count <2:
            x = [value for value in lines.split()]
            g = ('[%s]' % '  '.join(map(str, x)))
            y = g[1 : -1]
            print>> writefile, y
            count +=1
        else:
            break
    return temp
#energy reading and send to Y for adjusting

 
#yadjustment calculation
def yadjust(inputfile1, inputfile2, n):
    ymax = max(inputfile1)
    ymin = min(inputfile1)
    ydiff = ymax-ymin
    nparray = np.load('numpyexport.npy')
    yvalueadjusted = ytop - inputfile2+abs(ydiff)
    ymidpset = yvalueadjusted/ytotalunits
    nparray[2,n] = ymidpset
    fname = 'numpyexport'
    np.save(fname, nparray)
    return ydiff

#xadjustment calculation
def xadjust(inputfile1, yadjust, input1):
    xmax = max(inputfile1)
    xmin = min(inputfile1)
    xmid =  (xmax+xmin)/2
    global Xpagescale
    n = input1
    array = np.load('numpyexport.npy')
    neededydiff = yadjust+2
    xmidpoint = xmid
    xdiff = xmax-xmin
    ymid = array[2, n]
    ypastmid = array[2, n-1]
    presentydiff = abs(ymid) - abs(ypastmid)
    while xmidpoint < 10: 
        xmidpoint = xmidpoint+1
    if n == 0:
        xmidpoint = xmidpoint-4
    if n == 1:
        pass
    elif abs(presentydiff) < abs(neededydiff):
        xmidpoint = array[1, n-1] + xdiff + 2 
    else:
        pass 
    array[1, n] = xmidpoint
    fname = 'numpyexport'
    if xmidpoint > Xpagescale:
        Xpagescale = xmidpoint
    np.save(fname, array)
    return xmidpoint


#for dividing xyz file, removes atom column for handling by "numpy"
def divide(inputfile1, input1):
    n = input1
    temp22 = gbarray(inputfile1)     
    yfind = np.genfromtxt(StringIO(temp22), skiprows=2, delimiter=", ", usecols=(2))
    yfind2 = np.genfromtxt(StringIO(temp22), skiprows=1, delimiter=", ", usecols=(1),  skip_footer=atct)
    ydiff = yadjust(yfind, yfind2, n)
    A = np.genfromtxt(StringIO(temp22), skiprows=2, delimiter=", ", usecols=(1,2,3))
    xfind = np.genfromtxt(StringIO(temp22), skiprows=2, delimiter=", ", usecols=(1)) 
    xadjustment = xadjust(xfind, ydiff, n)
    moleculetrackingarray = np.load('numpyexport.npy')
    moleculetrackingarray[3, n] = yfind2
    yadjustment = moleculetrackingarray[2, input1]
    fname = 'numpyexport'
    np.save(fname, moleculetrackingarray)
    A[:,0] +=xadjustment
    A[:,1] +=yadjustment
    dt = np.dtype(str) 
    C = np.swapaxes(A,0,1)
    B = np.genfromtxt(StringIO(temp22), dtype="S8", skiprows=2, delimiter=", ", usecols=(0))
    D = np.vstack((B, C)) 
    E = np.swapaxes(D,0,1)
    saveE = scratch()
    np.savetxt(saveE, E, fmt="%s", delimiter="       ")
    return saveE

#combining files for final product
def combine(inputfile1, input1):
    n = input1
    file1 = Format(inputfile1)
    file2 = divide(inputfile1, n) 
    babelready = "%s.xyz" %scratch()
    destination = open(babelready, "wb")
    babelconverted = "%s.cdxml" %scratch()
    os.system('cat %s %s > %s' %(file1, file2, babelready)) # MAC SYSTEM CALL 
    conv =openbabel.OBConversion()
    conv.OpenInAndOutFiles("%s" %babelready,"%s" %babelconverted)
    conv.SetInAndOutFormats("xyz","cdxml")
    conv.Convert() 
    return "%s" %babelconverted

#Locates position of chemdraw molecule bounding box
def energyxyfind(array): 
    tf = scratch()
    conv = openbabel.OBConversion()
    conv.OpenInAndOutFiles(array, "%s" %tf)
    conv.SetInAndOutFormats("cdxml", "xyz")
    conv.Convert()
    return tf

#writes cdxml code for molecule in question
def energylabel(array, num):
    tf = energyxyfind(array)
    mNum = num+1
    xrow = np.genfromtxt(tf, skip_header=2, usecols=(1))
    yrow = np.genfromtxt(tf, skip_header=2, usecols=(2))
    ymax1 = max(yrow)
    xmin = min(xrow)
    ymax = ymax1+10
    #textbox = open('%s' %tf2, "w")
    #print tf2
    with open("scratch/tfile", "a") as myfile:
        myfile.write('<t p="%d %d"><s font="47136" size="12">Molecule %d</s></t>\n' %(xmin, ymax, mNum))
    #textbox = ('<t p="%d %d"><s font="47136" size="12">Molecule %d</s></t>/n' %(xmin, ymax, mNum))
    #os.system("cat %s >> scratch/tfile" %tf2) 


#running divided script through the process!
def runscript(inputfile1):
    m = addenergy(inputfile1) 
    xyzsplitcount(m) #responsible for the initial generation of the numpy array
    num = 0
    linenum = 0
    oxyz = open(xyzsplit(m), "r") #currently writes all of the split molecule file names to a single output
    pagespace() 
    for lines in oxyz.readlines():
        x = [value for value in lines.split()]
        g = ('[%s]' % '  '.join(map(str, x)))
        y = g[1 : -1]
        array = combine(y, num)
        oarray = open("%s" %array, "r")
        energylabel(array, num)
        for lines in oarray.readlines():
            if "<" in lines: 
                linenum +=1
            else:
                pass 
        if num == 0:
            crparray = open("%s" %array, "rw+")  
            os.system("head -n +%d %s  >> %s" %(linenum-1, array, "scratch/mfile.cdxml")) # MAC SYSTEM CALL 
        else:
            crparray = open("%s" %array, "rw+")  
            os.system("head -n +%d %s | tail -n %d >> %s" %(linenum-1, array, linenum-6, "scratch/mfile.cdxml")) # MAC SYSTEM CALL 
       # else:
       #     crparray = open("%s" %array, "rw+")  
       #     os.system("tail -n %d %s  >> %s" %(linenum-4, array,"scratch/mfile.cdxml")) # MAC SYSTEM CALL 
        num +=1
        linenum = 0
    os.remove('splitfiles.txt') # MAC SYSTEM CALL 
    os.remove('numpyexport.npy') # MAC SYSTEM CALL 
    return "scratch/mfile.cdxml"

#Program to add chemdraw header to the final output file (replaced the babelconvertdoubles program)
def finishingtouch(inputfile1):
    chemdrawfile = runscript(inputfile1)
    global molfoldvar
    cfile = "ElementaryStep_%d.cdxml" % molfoldvar 
    pageswide = int(Xpagescale/17.8)+1
    with open("%s" %(chemdrawfile)) as infile:
        with open(cfile, 'w') as outfile:
            for i, line in enumerate(infile):
                if i==0:
                    outfile.write("<?xml version='1.0'?>\n")
                elif i==1:
                    outfile.write("<CDXML BondLength='30'>\n")
                elif i==2:
                    outfile.write("<page\n")
                elif i==3:
                    outfile.write("HeightPages='%d'\n" % pages)
                elif i==4:
                    outfile.write("WidthPages='%d'>\n" % pageswide)
                else:
                    outfile.write(line)
    a = "scratch/tfile"
    os.system('cat %s >> %s' % (a, cfile))
    with open("%s" %cfile, "a") as infile:
         infile.write('</page></CDXML>')
    os.remove("%s" %chemdrawfile) # MAC SYSTEM CALL 
    os.system("rm -r -f scratch/") # MAC SYSTEM CALL 
    os.system("rm -r -f molecule_*") # MAC SYSTEM CALL 
    # os.system("open /Applications/CS\ ChemOffice\ 2010/ %s" % cfile)

finishingtouch(files)

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

