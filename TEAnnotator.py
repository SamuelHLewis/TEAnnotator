# Script to annotate the TE content of a genome using combined RepeatMasker & RepeatModeler
# Tested on program versions: RepeatMasker open-4.0.6, RepeatModeler open-1.0.8, RepeatScout 1.0.5, TRF 4.0.4

import os
import sys
import subprocess
import random as random
import shutil

# input genome fasta file
GenomeFile = 'dmel-X.fas'
if GenomeFile.endswith('.fasta'):
	os.rename(GenomeFile,GenomeFile.replace('.fasta','.fas'))
# number of cores to use
Cores = 8
# species to use in RepeatMasker
Species = 'Metazoa'

# function to run RepeatMasker on genome
def RepeatMaskerRunner(infile,cores=1,species='Metazoa',customlib=''):
	# run RepeatMasker (NB: E.coli insert check turned off)
	# if customlib not specified, find TEs according to species
	if customlib=='':
		# set outfile name early to reflect its origin from species
		outfile = infile.replace('.fas','_' + Species + '_TE.gff')
		cmd = 'RepeatMasker -pa ' + str(cores) + ' -no_is -species ' + species + ' -x -gff -a ' + infile
		print('Running RepeatMasker with species library ' + species)
		subprocess.call(cmd,shell=True)
	# if a custom repeat library is specified, find TEs according to the custom library
	else:
		# set outfile name early to reflect its origin from customlib
		outfile = infile.replace('.fas','_customlib_TE.gff')
		cmd = 'RepeatMasker -pa ' + str(cores) + ' -no_is -lib ' + customlib + ' -x -gff -a ' + infile
		print('Running RepeatMasker with custom library ' + customlib)
		subprocess.call(cmd,shell=True)
	# remove all output files apart from the gff
	os.remove(infile.replace('.fas','.fas.align'))
	os.remove(infile.replace('.fas','.fas.cat.gz'))
	os.remove(infile.replace('.fas','.fas.masked'))
	os.remove(infile.replace('.fas','.fas.out'))
	os.remove(infile.replace('.fas','.fas.tbl'))
	# this removes the RepeatMasker temp dir
	dirs = os.listdir(path='.')
	for i in dirs:
		if i.startswith('RM_'):
			shutil.rmtree(i)
	# in the gff, replace the TE RepeatMasker label ("similarity") with "exon" (so that gffread can read it)
	cmd = 'sed s/similarity/exon/ <' + infile.replace('.fas','.fas.out.gff') + ' > ' + infile.replace('.fas','_TEwithrib.gff')
	subprocess.call(cmd,shell=True)
	# remove the (now-redundant) initially-outputted gff
	os.remove(infile.replace('.fas','.fas.out.gff'))
	# screen out any annotations of rRNA
	cmd = 'grep -v \"rRNA\" ' + infile.replace('.fas','_TEwithrib.gff') + ' > ' + outfile
	subprocess.call(cmd,shell=True)
	os.remove(infile.replace('.fas','_TEwithrib.gff'))

# function to run RepeatModeler on genome
def RepeatModelerRunner(infile,cores=1):
	# build database
	cmd = 'BuildDatabase -name ' + infile.strip('.fas') + ' -engine ncbi ' + infile
	subprocess.call(cmd,shell=True)
	# run RepeatModeler
	cmd = 'RepeatModeler -engine ncbi -pa ' + str(cores) + ' -database ' + infile.strip('.fas')
	subprocess.call(cmd,shell=True)
	# move RepeatModeler output "consensi.fa" (will be in randomly-named subdir) to current directory
	dirs = os.listdir(path='.')
	for i in dirs:
		if i.startswith('RM_'):
			RepModDir = i
	cmd = "mv ./" + RepModDir + '/' + 'consensi.fa.classified ' + infile.replace('.fas','_consensi.fa.classified')
	subprocess.call(cmd,shell=True)
	# remove database files & RepeatModeler temp dir	
	os.remove(infile.replace('.fas','.nhr'))
	os.remove(infile.replace('.fas','.nin'))
	os.remove(infile.replace('.fas','.nnd'))
	os.remove(infile.replace('.fas','.nni'))
	os.remove(infile.replace('.fas','.nog'))
	os.remove(infile.replace('.fas','.nsq'))
	os.remove(infile.replace('.fas','.translation'))
	os.remove('unaligned.fa')
	shutil.rmtree(RepModDir)
	# run RepeatMasker on RepeatModeler output
	RepeatMaskerRunner(infile=GenomeFile,cores=Cores,customlib=infile.replace('.fas','_consensi.fa.classified'))


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

RepeatMaskerRunner(infile=GenomeFile,cores=Cores,species=Species)
# RepeatModelerRunner(infile=GenomeFile,cores=Cores)