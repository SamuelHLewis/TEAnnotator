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

#############################################################
##                      USER OPTIONS                       ##
## (Change nothing else unless you know what you're doing) ##
#############################################################

# input genome fasta file
GenomeFile = 'test.fas'
if GenomeFile.endswith('.fasta'):
	os.rename(GenomeFile,GenomeFile.replace('.fasta','.fas'))
# number of cores to use
Cores = 8
# species to use in RepeatMasker
Species = 'Metazoa'
# default cutoff value for RepeatMasker
CutOff = 250
# low complexity repeats included (True) or excluded (False)
NoLow = True
# minimum length of sequence
MinLength = 500

#####################################################
##                    CORE CODE                    ##
## (Don't touch unless you know what you're doing) ##
#####################################################

# function to run RepeatMasker on genome
def RepeatMaskerRunner(infile,cores=1,species='Metazoa',customlib='',cutoff=225,nolow=False):
	# run RepeatMasker (NB: E.coli insert check turned off)
	# if customlib not specified, find TEs according to species
	if customlib=='':
		# set outfile name early to reflect its origin from species
		outfile = infile.replace('.fas','_' + Species + '_TE.gff')
		if nolow is False:
			cmd = 'RepeatMasker -pa ' + str(cores) + ' -norna -cutoff ' + str(cutoff) + ' -no_is -species ' + species + ' -x -gff -a ' + infile + ' >> report.log'
		elif nolow is True:
			cmd = 'RepeatMasker -pa ' + str(cores) + ' -norna -cutoff ' + str(cutoff) + ' -nolow -no_is -species ' + species + ' -x -gff -a ' + infile + ' >> report.log'
		print('Running RepeatMasker with species library ' + species)
		subprocess.call(cmd,shell=True)
	# if a custom repeat library is specified, find TEs according to the custom library
	else:
		# set outfile name early to reflect its origin from customlib
		outfile = infile.replace('.fas','_customlib_TE.gff')
		if nolow is False:
			cmd = 'RepeatMasker -pa ' + str(cores) + ' -norna -cutoff ' + str(cutoff) + ' -no_is -lib ' + customlib + ' -x -gff -a ' + infile + ' >> report.log'
		elif nolow is True:
			cmd = 'RepeatMasker -pa ' + str(cores) + ' -norna -cutoff ' + str(cutoff) + ' -nolow -no_is -lib ' + customlib + ' -x -gff -a ' + infile + ' >> report.log'
		print('Running RepeatMasker with custom library ' + customlib)
		subprocess.call(cmd,shell=True)
		os.remove(customlib)
	# remove all output files apart from the gff
	os.remove(infile.replace('.fas','.fas.align'))
	os.remove(infile.replace('.fas','.fas.masked'))
	os.remove(infile.replace('.fas','.fas.out'))
	os.remove(infile.replace('.fas','.fas.tbl'))
	# check whether the .cat file has been compressed or not before removing
	if os.path.isfile(infile.replace('.fas','.fas.cat.gz')) is True:
		os.remove(infile.replace('.fas','.fas.cat.gz'))
	if os.path.isfile(infile.replace('.fas','.fas.cat')) is True:
		os.remove(infile.replace('.fas','.fas.cat'))	
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
	return(outfile)

# function to run RepeatModeler on genome
def RepeatModelerRunner(infile,cores=1):
	# build database
	cmd = 'BuildDatabase -name ' + infile.strip('.fas') + ' -engine ncbi ' + infile + ' >> report.log'
	subprocess.call(cmd,shell=True)
	# run RepeatModeler
	cmd = 'RepeatModeler -engine ncbi -pa ' + str(cores) + ' -database ' + infile.strip('.fas') + ' >> report.log'
	subprocess.call(cmd,shell=True)
	# move RepeatModeler output "consensi.fa" (will be in randomly-named subdir) to current directory
	dirs = os.listdir(path='.')
	for i in dirs:
		if i.startswith('RM_'):
			RepModDir = i
	# remove database files	
	os.remove(infile.replace('.fas','.nhr'))
	os.remove(infile.replace('.fas','.nin'))
	os.remove(infile.replace('.fas','.nnd'))
	os.remove(infile.replace('.fas','.nni'))
	os.remove(infile.replace('.fas','.nog'))
	os.remove(infile.replace('.fas','.nsq'))
	os.remove(infile.replace('.fas','.translation'))
	# if RepeatModeler has identified some TEs, move the RepBase-formatted file into the working directory and rename it, remove the other intermediate files & temp dir, and return 'TEs found'
	if os.path.isfile('./' + RepModDir + '/consensi.fa.classified') is True:
		cmd = "mv ./" + RepModDir + '/consensi.fa.classified ' + infile.replace('.fas','_consensi.fa.classified')
		subprocess.call(cmd,shell=True)
		os.remove('unaligned.fa')
		if os.path.isfile('tmpConsensi.fa.cat.all') is True:
			os.remove('tmpConsensi.fa.cat.all')
		shutil.rmtree(RepModDir)
		return(infile.replace('.fas','_consensi.fa.classified'))
	# if RepeatModeler hasn't identified any TEs, remove the temp dir and return 'No TEs found'
	else:
		shutil.rmtree(RepModDir)
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
	# extract sequences corresponding to gff annotations to a fasta file (NB: this will merge overlapping annotations)
	cmd = 'gffread ' + gff + ' -g ' + fastafile + ' -w ' + fastafile.replace('.fas','_TE.fas')
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
	# remove report logfile if already exists
	if os.path.isfile('report.log') is True:
		os.remove('report.log')
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
	

# the important call (returns gff and optional fasta for TEs in fastafile)
TEAnnotator(fastafile=GenomeFile)











