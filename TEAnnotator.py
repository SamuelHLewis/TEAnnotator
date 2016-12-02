#!/usr/bin/env python3

# This annotates the TE content of a genome using a combination of RepeatMasker & RepeatModeler methods
# Dependencies (all of which need to be in PATH): 
	# RepeatMasker open-4.0.6
	# RepeatModeler open-1.0.8
	# RepeatScout 1.0.5
	# TRF 4.0.4
	# bedtools v2.17.0
	# cufflinks v2.2.1 

import os
import sys
import subprocess
import random as random
import shutil
import re

import argparse

#############################
## DEFAULT ARGUMENT VALUES ##
#############################
# input genome fasta file
GenomeFile = ''
# number of cores to use
Cores = 1
# species to use in RepeatMasker
Species = 'Metazoa'
# default cutoff value for RepeatMasker
CutOff = 250
# minimum length of repeat annotation
MinLength = 500
# low complexity repeats excluded (True) or included (False)
NoLow = True

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-i', '--infile', type=str, help='input fasta file')
parser.add_argument('-c', '--cores',  type=int, help='number of cores')
parser.add_argument('-s', '--species',  type=str, help='species to use for RepeatMasker')
parser.add_argument('-u', '--cutoff',  type=int, help='cutoff value for RepeatMasker')
parser.add_argument('-m', '--minlength',  type=int, help='minimum length of repeat annotation')
parser.add_argument('--lowcomplex', help='include low complexity repeats', action="store_true")
args = parser.parse_args()
# input fasta file parsing
GenomeFile = args.infile
if GenomeFile is not None:
	if GenomeFile.endswith('.fas'):
		print('Input file is ' + GenomeFile)
	elif GenomeFile.endswith('.fasta'):
#		os.rename(GenomeFile,GenomeFile.replace('.fasta','.fas'))
		GenomeFile = GenomeFile.replace('.fasta','.fas')
		print('Input file renamed to ' + GenomeFile)
	elif GenomeFile.endswith('.fa'):
#		os.rename(GenomeFile,GenomeFile.replace('.fa','.fas'))
		GenomeFile = GenomeFile.replace('.fa','.fas')
		print('Input file renamed to ' + GenomeFile)
	else:
		print('Input file does not have normal fasta filename ending - please check file format')
		sys.exit(0)
else:
	print('ERROR: no input file specified')
# cores parsing
if args.cores is None:
	print('Using default number of cores (' + str(Cores) + ')')
else:
	Cores = args.cores
	if Cores > 0:
		print(str(Cores) + ' cores specified')
	else:
		print('ERROR: 1 or more cores (-c) required')
		sys.exit(0)
# species parsing
if args.species is None:
	print('Using default species (' + Species + ')')
else:
	Species = args.species
	print('RepeatMasker species set to ' + Species)
# cutoff value parsing
if args.cutoff is None:
	print('Using default RepeatMasker cutoff value (' + str(CutOff) + ')')
else:
	if args.cutoff > 0:
		CutOff = args.cutoff
		print('RepeatMasker cutoff value set to ' + str(CutOff))
	else:
		print('ERROR: cutoff value (-u) must be higher than 0')
		sys.exit(0)
# minimum annotation length parsing
if args.minlength is None:
	print('Using default minimum annotation length (' + str(MinLength) + ')')
else:
	if args.minlength > 0:
		MinLength = args.minlength
		print('Minimum annotation length set to ' + str(MinLength))
	else:
		print('ERROR: minimum annotation length (-m) must be more than 0')
		sys.exit(0)
# low complexity repeat handling parsing
if args.lowcomplex:
	NoLow = False
	print('Including low complexity repeats')
else:
	print('Excluding low complexity repeats')

#####################################################
##                    CORE CODE                    ##
## (Don't touch unless you know what you're doing) ##
#####################################################

# function to run RepeatMasker on genome
def RepeatMaskerRunner(infile,cores=1,species='Metazoa',customlib='',cutoff=225,nolow=True):
	# run RepeatMasker (NB: E.coli insert check turned off)
	# if customlib not specified, find TEs according to species
	if customlib=='':
		# set outfile name early to reflect its origin from species
		outfile = infile.replace('.fas','_' + Species + '_TE.gff')
		if nolow is False:
			cmd = 'RepeatMasker -pa ' + str(cores) + ' -norna -cutoff ' + str(cutoff) + ' -no_is -species ' + species + ' -x -gff -a ' + infile
		elif nolow is True:
			cmd = 'RepeatMasker -pa ' + str(cores) + ' -norna -cutoff ' + str(cutoff) + ' -nolow -no_is -species ' + species + ' -x -gff -a ' + infile
		print('Running RepeatMasker with species library ' + species)
		subprocess.call(cmd,shell=True)
		# make a new dir in TEAnnotator to put the tempfiles into
		tempdir = 'TEAnnotator/'+Species+'/'
		if os.path.isdir(tempdir) is True:
			shutil.rmtree(tempdir)
			print('Old ' + tempdir + ' removed')
		os.makedirs(tempdir)
		# move all output files apart from gff to species-specific dir in TEAnnotator/
		dirs = os.listdir(path='.')
		for i in dirs:
			if i.startswith('RM_'):
				# RepeatMaskerDir = i.split('/')[-2]
				# shutil.move(RepeatMaskerDir,tempdir+RepeatMaskerDir)
				None
			else:
				filename = os.path.split(i)[1]
				if filename.endswith('.fas.align'):
					shutil.move(filename,tempdir+filename)
				elif filename.endswith('.fas.cat'):
					shutil.move(filename,tempdir+filename)
				elif filename.endswith('.fas.cat.gz'):
					shutil.move(filename,tempdir+filename)
				elif filename.endswith('.fas.masked'):
					shutil.move(filename,tempdir+filename)
				elif filename.endswith('.fas.out'):
					shutil.move(filename,tempdir+filename)
				elif filename.endswith('.fas.tbl'):
					shutil.move(filename,tempdir+filename)
	# if a custom repeat library is specified, find TEs according to the custom library
	else:
		# set outfile name early to reflect its origin from customlib
		outfile = infile.replace('.fas','_customlib_TE.gff')
		if nolow is False:
			cmd = 'RepeatMasker -pa ' + str(cores) + ' -norna -cutoff ' + str(cutoff) + ' -no_is -lib ' + customlib + ' -x -gff -a ' + infile
		elif nolow is True:
			cmd = 'RepeatMasker -pa ' + str(cores) + ' -norna -cutoff ' + str(cutoff) + ' -nolow -no_is -lib ' + customlib + ' -x -gff -a ' + infile
		print('Running RepeatMasker with custom library ' + customlib)
		subprocess.call(cmd,shell=True)
		# make a new dir in TEAnnotator to put the tempfiles into
		tempdir = 'TEAnnotator/Custom/'
		if os.path.isdir(tempdir) is True:
			shutil.rmtree(tempdir)
			print('Old ' + tempdir + ' removed')
		os.makedirs(tempdir)
		# move all output files apart from gff to Custom dir in TEAnnotator/
		dirs = os.listdir(path='.')
		for i in dirs:
			if i.startswith('RM_'):
				RepeatMaskerDir = i.split('/')[-2]
				shutil.move(RepeatMaskerDir,tempdir+RepeatMaskerDir)
			else:
				filename = os.path.split(i)[1]
				if filename.endswith('.fas.align'):
					shutil.move(filename,tempdir+filename)
				elif filename.endswith('.fas.cat'):
					shutil.move(filename,tempdir+filename)
				elif filename.endswith('.fas.cat.gz'):
					shutil.move(filename,tempdir+filename)
				elif filename.endswith('.fas.masked'):
					shutil.move(filename,tempdir+filename)
				elif filename.endswith('.fas.out'):
					shutil.move(filename,tempdir+filename)
				elif filename.endswith('.fas.tbl'):
					shutil.move(filename,tempdir+filename)
	# in the gff, give each annotation a unique name based on start-end positions, and replace the TE RepeatMasker label ("similarity") with "CDS" (so that gffread can extract coding seqs)
	newgffcontents = ''
	for line in open(infile.replace('.fas','.fas.out.gff')):	
		if line.startswith('#'):
			None
		else:
			linesplit = line.split('\t')
			TEname = linesplit[8].replace('Target \"Motif:','ID=').replace('\"','').replace(' ','_').strip('\n')
			if re.search('rRNA',TEname) == None:
				newline = linesplit[0] + '\t' + linesplit[1] + '\t' + 'CDS' + '\t' + linesplit[3] + '\t' + linesplit[4] + '\t' + linesplit[5] + '\t' + linesplit[6] + '\t' + linesplit[7] + '\t' + TEname + '\n'
				newgffcontents += newline
	newgff = open(outfile,'wt')
	newgff.write(newgffcontents)
	newgff.close()
	return(outfile)

# function to run RepeatModeler on genome
def RepeatModelerRunner(infile,cores=1):
	# build database
	cmd = 'BuildDatabase -name ' + infile.strip('.fas') + ' -engine ncbi ' + infile
	subprocess.call(cmd,shell=True)
	# run RepeatModeler
	cmd = 'RepeatModeler -engine ncbi -pa ' + str(cores) + ' -database ' + infile.strip('.fas')
	subprocess.call(cmd,shell=True)
	# make a new dir in TEAnnotator to put the tempfiles into
	tempdir = 'TEAnnotator/RepeatModeler/'
	if os.path.isdir(tempdir) is True:
		shutil.rmtree(tempdir)
		print('Old ' + tempdir + ' removed')
	os.makedirs(tempdir)
	dirs = os.listdir(path='.')
	for i in dirs:
		# identify the RepeatModeler dir
		if i.startswith('RM_'):
			RepModDir = i
		# move database indexing files to TEAnnotator/RepeatModeler/ dir 
		else:
			filename = os.path.split(i)[1]
			if filename.endswith('.nhr'):
				shutil.move(filename,tempdir+filename)
			elif filename.endswith('.nin'):
				shutil.move(filename,tempdir+filename)
			elif filename.endswith('.nnd'):
				shutil.move(filename,tempdir+filename)
			elif filename.endswith('.nni'):
				shutil.move(filename,tempdir+filename)
			elif filename.endswith('.nog'):
				shutil.move(filename,tempdir+filename)
			elif filename.endswith('.nsq'):
				shutil.move(filename,tempdir+filename)
			elif filename.endswith('.translation'):
				shutil.move(filename,tempdir+filename)
	# if RepeatModeler has identified some TEs, move the RepBase-formatted file ("consensi.fa" in randomly-named subdir) into the working directory and rename it, remove the other intermediate files & temp dir, and return 'TEs found'
	if os.path.isfile('./' + RepModDir + '/consensi.fa.classified') is True:
		cmd = "mv ./" + RepModDir + '/consensi.fa.classified ' + infile.replace('.fas','_consensi.fa.classified')
		subprocess.call(cmd,shell=True)
		dirs = os.listdir(path='.')
		for i in dirs:
			filename = os.path.split(i)[1]
			if filename.endswith('unaligned.fa'):
				shutil.move(filename,tempdir+filename)
			elif filename.endswith('tmpConsensi.fa.cat.all'):
				shutil.move(filename,tempdir+filename)
		# move the RepeatModeler temp dir into the overall tempfile dir
		shutil.move(RepModDir,tempdir+RepModDir)
		return(infile.replace('.fas','_consensi.fa.classified'))
	# if RepeatModeler hasn't identified any TEs, remove the temp dir and return 'No TEs found'
	else:
		# move the RepeatModeler temp dir into the overall tempfile dir
		shutil.move(RepModDir,tempdir+RepModDir)
		return('No TEs found')

# function to combine 2 gffs, and output non-redundant gff file
def GffCombiner(fastafile,gff1,gff2):
	combinedgff = fastafile.replace('.fas','_TEcombined.gff')
	# combine gff1 & gff2 into one gff file
	cmd = 'cat ' + gff1 + ' >> ' + combinedgff
	subprocess.call(cmd,shell=True)
	cmd = 'cat ' + gff2 + ' >> ' + combinedgff
	subprocess.call(cmd,shell=True)
	print('gff files combined')
	# sort combined gff by chromosome & ascending order within chromosome
	cmd = 'sortBed -i ' + combinedgff + ' > ' + combinedgff.replace('.gff','_sorted.gff')
	subprocess.call(cmd,shell=True)
	print('Combined gff file sorted')
	# remove temp/intermediate files
	os.remove(combinedgff)
	os.rename(combinedgff.replace('.gff','_sorted.gff'),fastafile.replace('.fas','_TE.gff'))
	return(fastafile.replace('.fas','_TE.gff'))

# function to extract sequences from fasta file according to annotations in gff file
def Extractor(fastafile,gff):
	# extract sequences corresponding to gff annotations to a fasta file (NB: this forces strandedness i.e. reverse-complements annotations on the antisense strand)
	cmd = 'bedtools getfasta -s -fo ' + fastafile.replace('.fas','_TE.fas') + ' -fi ' + fastafile + ' -bed ' + gff
	subprocess.call(cmd,shell=True)
	print('Fasta file written: ' + fastafile.replace('.fas','_TE.fas'))
	# remove temp/intermediate files
	os.remove(fastafile + '.fai')

# function to apply a length screen to a gff file
def GffLengthScreener(gff,minlength=0,maxlength=float('inf')):
	ScreenedGff = ''
	counter = 0
	included = 0
	for line in open(gff,'r'):
		counter += 1
		start = int(line.split('\t')[3])
		end = int(line.split('\t')[4])
		# NB +1 here to include all bases of annotation
		length = end-start+1
		if minlength < length < maxlength:
			ScreenedGff += line
			included += 1
	outfile = open(gff.replace('.gff','_screened.gff'),'wt')
	outfile.write(ScreenedGff)
	outfile.close()
	print('New gff written: ' + str(included) + '/' + str(counter) + ' annotations passed length filter')
	return(gff.replace('.gff','_screened.gff'))

# function to annotate TEs in a fasta file (usually a genome) using RepeatMasker & RepeatModeler - produces gff by default and additional (optional) fasta
def TEAnnotator(fastafile=GenomeFile,outputfasta=True):
	# remove TEAnnotator dir if already exists and replace it with new (blank) dir
	if os.path.isdir('TEAnnotator') is True:
		shutil.rmtree('TEAnnotator')
		print('Old TEAnnotator tempfiles removed')
	os.makedirs('TEAnnotator')
	# run RepeatMasker using species set at start
	RepeatMaskerGFF = RepeatMaskerRunner(infile=GenomeFile,cores=Cores,species=Species,cutoff=CutOff,nolow=NoLow)
	# run RepeatModeler on genome
	RepeatModelerOutput = RepeatModelerRunner(infile=GenomeFile,cores=Cores)
	if outputfasta is True:
		# if RepeatModeler identified no TEs, extract sequences based on the length-screened RepeatMasker gff file
		if RepeatModelerOutput == 'No TEs found':
			print('No TEs found by RepeatModeler - extracting TE sequences based on length-screened species-based gff')
			# screen out TE annotations by length
			FinalGff = GffLengthScreener(gff=RepeatMaskerGFF,minlength=MinLength)
			os.remove(RepeatMaskerGFF)
			Extractor(fastafile=GenomeFile,gff=FinalGff)
		# if RepeatModeler identified some TEs, run RepeatMasker based on the RepeatModeler de novo database, combine the RepeatMasker and RepeatModeler gff files, and extract sequences based on the combined gff file
		else:
			print('Running RepeatMasker based on RepeatModeler models')
			RepeatModelerGFF = RepeatMaskerRunner(infile=GenomeFile,cores=Cores,customlib=RepeatModelerOutput,cutoff=CutOff,nolow=NoLow)
			print('Combining species-based and RepeatModeler-based gff files')
			GffCombinerGFF = GffCombiner(fastafile=GenomeFile,gff1=RepeatMaskerGFF,gff2=RepeatModelerGFF)
			os.remove(RepeatMaskerGFF)
			os.remove(RepeatModelerGFF)
			# screen out TE annotations by length
			FinalGff = GffLengthScreener(gff=GffCombinerGFF,minlength=MinLength)
			print('Extracting TE sequences based on length-screened combined gff')
			os.remove(GffCombinerGFF)
			Extractor(fastafile=GenomeFile,gff=FinalGff)

# find TEs in genome
TEAnnotator(fastafile=GenomeFile)

