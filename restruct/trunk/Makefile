# Use this Makefile with make

# Executable name
CMD = restruct.exe

# -------- description of DFLAGS ---------------


# -------- Define environmental variable C_COMPILER -----------
# Make sure it is defined
#          ifeq ($(strip$(FORTRAN_COMPILER)),)
# Otherwise you can define it here also by uncommenting next line
 FC = g++ 
DFLAGS = -DGEOMETRIC #-Define the cpp flags to be used
#DFLAGS = -DGEOMETRIC #-Define the cpp flags to be used
OFLAGS =  # optimization
LINKERFLAGS = -lm # link math libraries



#
# Implicit rules for handling src files
#  ( uses static pattern rules, see info for make )
.c.o:
	$(FC) -c -g $(DFLAGS) -Wimplicit $<
.cpp.o:
	$(FC) -c -g $(DFLAGS) $<

OBJECTS = zstruct.o main.o icoord.o tm.o pTable.o stringtools.o  mm_grad.o utils.o mem.o opt.o print.o geombasis.o atom.o iso.o mopac.o paths.o pgsm.o dft.o

$(CMD) : $(OBJECTS)
	$(FC) $(DEBUG_FLAGS) $(OFLAGS) $(LINKERFLAGS) $(OBJECTS)  -o ./$(CMD)

clean:
	/bin/rm -f *.o *.i *.mod *.exe make.log

cleano:
	rm -f *.o *.i

depend :
	g++ -MM *.cpp *.c >> Makefile 

# DO NOT DELETE created with g++ -MM *.cpp *.c
main.o: main.cpp zstruct.h icoord.h
paths.o: paths.cpp zstruct.h icoord.h
pgsm.o: pgsm.cpp pgsm.h constants.h
dft.o: dft.cpp dft.h constants.h
iso.o: iso.cpp zstruct.h icoord.h geombasis.h constants.h utils.h mopac.h
atom.o: atom.cpp atom.h icoord.h
zstruct.o: zstruct.cpp zstruct.h constants.h print.h
icoord.o: icoord.cpp icoord.h zstruct.h
mm_grad.o: mm_grad.cpp icoord.h constants.h
mopac.o: mopac.cpp mopac.h constants.h
geombasis.o: geombasis.cpp geombasis.h icoord.h
mem.o: mem.cpp icoord.h
opt.o: opt.cpp icoord.h 
pTable.o: pTable.cpp pTable.h
print.o: print.cpp icoord.h
stringtools.o: stringtools.cpp stringtools.h
tm.o: tm.cpp icoord.h utils.h
utils.o: utils.cpp utils.h icoord.h

