# Script to annotate the TE content of a genome using combined RepeatMasker & RepeatModeler

import os
import sys
import subprocess
import random as random

# input genome fasta file
GenomeFile = 'dmel-X.fas'
# number of cores to use
Cores = 8
# species to use in RepeatMasker
Species = 'Metazoa'

# function to run RepeatMasker on genome
def RepeatMaskerRunner(infile,cores,species):
	# run RepeatMasker (NB: E.coli insert check turned off)
	cmd = 'RepeatMasker -pa ' + str(cores) + ' -no_is -species \"' + species + '\" -x -gff -a ' + infile
	subprocess.call(cmd,shell=True)
	# remove all output files apart from the gff
	os.remove(infile.replace('.fas','.fas.align'))
	os.remove(infile.replace('.fas','.fas.cat.gz'))
	os.remove(infile.replace('.fas','.fas.masked'))
	os.remove(infile.replace('.fas','.fas.out'))
	os.remove(infile.replace('.fas','.fas.tbl'))
	# in the gff, replace the TE RepeatMasker label ("similarity") with "exon" (so that gffread can read it)
	cmd = 'sed s/similarity/exon/ <' + infile.replace('.fas','.fas.out.gff') + ' > ' + infile.replace('.fas','_TE.gff')
	subprocess.call(cmd,shell=True)
	# remove the (now-redundant) initially-outputted gff
	os.remove(infile.replace('.fas','.fas.out.gff'))
	# screen out any annotations of rRNA
	cmd = 'grep -v \"rRNA\" ' + infile.replace('.fas','_TE.gff') + ' > ' + infile.replace('.fas','_TE_clean.gff')
	subprocess.call(cmd,shell=True)
	os.remove(infile.replace('.fas','_TE.gff'))


# function to combine multiple gffs, and output non-redundant gff and fasta files
	# # RepeatMasker will give same elements identical names - therefore append random numbers to end of all element names to allow parsing
	# numberer(infile.replace('.fas','_TE_clean.gff'))
	# os.remove(infile.replace('.fas','_TE_clean.gff'))
	# os.rename(infile.replace('.fas','_TE_clean_renamed.gff'),infile.replace('.fas','_TE_RepeatMasker.gff'))
	# # write fasta file of TEs
	# cmd = 'gffread ' + infile.replace('.fas','_TE.gff') + ' -g ' + infile + ' -w ' + infile.replace('.fas','_TE.fas')
	# subprocess.call(cmd,shell=True)

# function to append a random number to the end of identical entries in a gff file, and output to file
def numberer(infile):
	for i in open(infile,'r'):
		name = ''
		newname = ''
		line = i.split('\t')
		for j in line:
			if j.startswith('Target'):
				temp = j.split(' ')
				name = temp[0]
				newname = name + '_' + str(random.randint(1,1000000000000000000000000000000000000))
		newline = i.replace(name,newname)
		outfile = open(infile.replace('.gff','_renamed.gff'),'a')
		outfile.write(newline)
		outfile.close()

RepeatMaskerRunner(GenomeFile,Cores,Species)

#test