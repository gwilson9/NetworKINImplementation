#!/usr/bin/env python

"""
To be used for generating the blast_out.txt file for NetworKIN

"""




import argparse
import fileinput
import os
import csv
import sys

parser = argparse.ArgumentParser(description ='Generate NetworKIN Blast Input File')

parser.add_argument('-f', "--fasta", type=str,required = True, help="Declare the fasta file")
parser.add_argument('-o', "--organism", type =str, required = True, help="Declare the organism. Mouse (10090) or human (9606)")

args = parser.parse_args()




# Run system binary
def myPopen(cmd):
	try:
		pipe = subprocess.Popen(cmd, shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout = pipe.stdout.readlines()
	except:
		sys.stderr.write('ERROR executing: '+`cmd`+'\n')
		sys.exit()
	else:
		return stdout


command = "%s -a %s -p blastp -e 1e-10 -m 8 -d %s -i %s | sort -k12nr -> blast_out_Ensembl74.txt"%(blastDir, number_of_processes, blastDB, blast_tmpfile.name)

sys.stderr.write("Performing blast search\n")
sys.stderr.write("fn_blast_output: %s"%fn_blast_output)
blast_out = myPopen(command)
