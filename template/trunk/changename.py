#!/usr/bin/python

import re,os,sys,commands

list_string_files = commands.getoutput('ls stringfile*')


#length = len(list_string_files.split('\n'))/3
#j = range(0,length)

#list_directories = commands.getoutput('ls -d q*.xyz')
#zs=''
#list_directories = list_directories.split('\n')
#for i in list_directories:
#	os.chdir('/export/zimmerman/jnmetz/Zstruct_Files/' + i)
#	outp = str(commands.getoutput('./status0'))
#	outp = outp.split('\n')
#	new=''
#	m = 0001
#	for i in outp:
#		new += str( "%03d" % (m,)) + ' ' + i + '\n'
#		m += 1
#	new = new.split('\n')
#	new = new[0:-1]
#	zs=''
#	for i in new:
#		if 'XTS' in i:
#			zs += new[int(i[0:4])-3] + '\n'
#zs = zs.split('\n')
#zs = zs[0:-1]
#y=y2=y3=y4=''
#for i in zs:
#	if int(i[4:8]) == 0001:
#		y2 += i[4:8] + '\n'
#		y3 += i[24:] + '\n'
#	else:
#		y += i[4:8] +'\n'
#		y4 += i[24:] + '\n'
#print y3
#print y4
#list_directories2 = str('0000') + '\n' + str('0001')
#list_directories2 = list_directories2.split('\n')
#for i in list_directories2:
#	for x in list_string_files:
	#	os.chdir('/export/zimmerman/jnmetz/Zstruct_Files/' + x +'/'+ str(i))
#		k2=k3=''
#		if int(i) == 0001:
#			for i in (y3.split('\n'))[0:-1]:
#				k2 += 'stringfile' + i + '.xyz' + '\n'
#				k3 += 'stringfile.xyz' + i + '\n'
#		else:
#			for i in (y4.split('\n'))[0:-1]:
#				k2 += 'stringfile' + i + '.xyz' + '\n'
#				k3 += 'stringfile.xyz' + i  + '\n'
#
#		k2 = k2.split('\n')
#		k2 = k2[0:-1]
#		k3 = k3.split('\n')
#		k3 = k3[0:-1]
#		print k3
#		list_string_files = commands.getoutput('ls stringfile*')
#		k4 = ''
#		print list_string_files
#		for i in k3:
#			if i in list_string_files:
#				with open(i,'r') as inf:
#					k4 += 'end'
#					for num, line in enumerate(inf):
#						k4 += line
#		k4 = k4.split('end')
#		k4 = k4[1:]	
#		print k4[0]
#		m = 0
#		for x in k2:
#	
#			if m > len(k4):
#				break
#			else:
#				with open(x,'w') as outf:
#					print >> outf, k4[m]
#					m += 1
# 
##
