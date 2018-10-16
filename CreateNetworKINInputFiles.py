#!/usr/bin/env python

"""
To be used to convert UniprotIDs to Ensembl74 ENSPs for use with NetworKIN.py

usage Mouse python UniprotIDsToEnsembl74.py -o 10090 ResultsFile.txt

usage Human python UniprotIDsToEnsembl74.py -o 9606 ResultsFile.txt

"""

import argparse
import fileinput
import os
import csv
import re
import sys
import tempfile
import subprocess, fpformat, dircache, random, operator, glob, thread, threading
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description ='Generate NetworKIN Phosphosites Input File')

    parser.add_argument('-o', "--organism", type=str,required = True, help="Declare whether data is Mouse (10090) or Human (9609)")
    parser.add_argument('-r', "--results", type =str, required = True, help="Provide relative path to results file")
    parser.add_argument('-f', "--fasta", type=str, required=True, help="Provide relative path to results file")

    args = parser.parse_args()
    organism = args.organism
    resultsFile = args.results
    fasta = args.fasta
    dir = os.path.dirname(os.path.abspath(__file__))

    sys.stderr.write("Path: %s\n"%dir)
    #sys.stderr.write("Organism: %s\n"%args.organism)
    sys.stderr.write("results: %s\n"%resultsFile)

    if args.organism == "10090" or args.organism == "9606":
        sys.stderr.write("Organism: %s\n"%args.organism)
    else:
        raise Exception("Organism not recognized.\nMust be 9606 for human data or 10090 for mouse data.")

    #uniprotDict[UniprotID] = Ensembl74 ENSP
    #geneNameDict[Gene Name] = Ensembl74 ENSP
    #uniprotDict is prioritized as there are more redundant gene names
    uniprotDict, geneNameDict = makeDicts(organism, dir)

    # ProteinID, Position within Protein, Residue (S,T,Y)
    id_pos_res = writePhosphoFile(dir, uniprotDict, geneNameDict,resultsFile)

    id_seq = readFasta(fasta, uniprotDict, geneNameDict, dir, id_pos_res)

    writeBlastOut(dir, fasta, uniprotDict, geneNameDict, id_pos_res, id_seq)

def myPopen(cmd):
    try:
        pipe = subprocess.Popen(cmd, shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = pipe.stdout.readlines()
    except:
        sys.stderr.write("Error executing: %s"%cmd)
        sys.exit()
    else:
        return stdout

# Will perform blast search between input fasta and the
def writeBlastOut(dir, fasta, uniprotDict, geneNameDict, id_pos_res, id_seq):

    tempfile.tempdir = dir + "/tmp"

    blast_tmpfile = tempfile.NamedTemporaryFile()
    blastDir = dir + "/blast-2.2.17/bin/blastall"
    blastDB = dir + "/" + fasta +"_Ensembl74.fasta"
    number_of_processes = 1

    if not os.path.isfile(blastDB+'.pin'):
        command = "%s/formatdb -i %s"%(blastDir.rsplit("/",1)[0], blastDB)
        sys.stderr.write("Initializing Blast Database\nCommand: %s\n"%command)
        myPopen(command)

    missedIDs = {}
    for id in id_pos_res:
        #sys.stderr.write("id_pos_res id: %s\n"%id)
        try:
            blast_tmpfile.write('>'+id+'\n'+id_seq[id]+'\n')
            #sys.stderr.write("Found Sequence\n")
        except:
            missedIDs[id] = id
            sys.stderr.write("No sequence available for '%s'\n"%(id))

    blast_tmpfile.flush()

    sys.stderr.write("NetworKIN scores not available for %s of %s proteins\n"%(len(missedIDs),len(id_pos_res)))

    #import pdb; pdb.set_trace()

    command = "%s -a %s -p blastp -e 1e-10 -m 8 -d %s -i %s | sort -k12nr -> blast_out_Ensembl74.txt"%(blastDir, number_of_processes, blastDB, blast_tmpfile.name)

    sys.stderr.write("Performing Blast Search\nCommand: %s"%command)
    myPopen(command)


# Read sequences from fasta file
# Function will replace fasta header with corresponding EnsemblID
def readFasta(fastaFile, uniprotDict, geneNameDict,dir, id_pos_res):

    id_seq = {}
    newRecords = []

    with open(dir + "/" + fastaFile, "rU") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            id = record.id
            desc = record.description
            uniprotID = desc.split('|')[1].split('-')[0]
            if "GN=" in desc:
                newID = str.replace(desc, "GN=","*")
                geneName = newID.split('*')[1].split()[0]

            if uniprotID in uniprotDict:
                record.id = uniprotDict[uniprotID]
                record.description = uniprotDict[uniprotID]
                if record.id in id_pos_res:
                    newRecords.append(record)
                    id_seq[record.id] = str(record.seq)
            else:
                try:
                    if geneName in geneNameDict:
                        #sys.stderr.write("Using gene name for: %s\n"%geneName)
                        record.id = geneNameDict[geneName]
                        record.description = geneNameDict[geneName]
                        if record.id in id_pos_res:
                            rewRecords.append(record)
                            id_seq[record.id] = str(record.seq)
                except:
                    pass
                    #sys.stderr.write("No match for %s\n"%desc)

    sys.stderr.write("Length newRecords: %s\n"%len(newRecords))

    SeqIO.write(newRecords, fastaFile + "_Ensembl74.fasta", "fasta")

    return id_seq

def makeDicts(organism, dir):
    uniprotDict = {}
    geneNameDict = {}

    with open(dir +"/Data/20181002_Biomart_IdentifierConversion_Human_Mouse.csv") as infile:
        reader = csv.DictReader(infile)
        for line in reader:
            if organism == "9606":
                if not line["Uniprot_Human"] in (None, ""):
                    uniprotDict[line["Uniprot_Human"]] = line["Ensembl74_Human"]
                if not reader["HGNC"] in (None, ""):
                    geneNameDict[line["HGNC"]] = line["Ensembl74_Human"]
            else:
                if not line["Uniprot_Mouse"] in (None, ""):
                    uniprotDict[line["Uniprot_Mouse"]] = line["Ensembl74_Human"]
                if not line["MGI"] in (None, ""):
                    geneNameDict[line["MGI"]] = line["Ensembl74_Human"]

    return uniprotDict, geneNameDict

#sys.stderr.write("Length uniprotDict: %s\n"%len(uniprotDict))
#sys.stderr.write("Length geneNameDict: %s\n"%len(geneNameDict))

def writePhosphoFile(dir, uniprotDict, geneNameDict, results):

    id_pos_res = {}

    with open(dir + "/PhosphoSites.txt", 'w+') as outfile:
        writer = csv.writer(outfile , delimiter = "\t")
        with open(dir + "/" + results) as infile:
            reader = csv.DictReader(infile)
            for line in reader:
                if line["Majority Protein ID"].split('-')[0] in uniprotDict:
                    id = uniprotDict[line["Majority Protein ID"].split('-')[0]]
                    pos = line["Position within protein"]
                    res = line["Amino Acid"]

                    if id in id_pos_res:
                        id_pos_res[id][pos] = results
                    else:
                        id_pos_res[id] = {pos: res}

                    writer.writerow([uniprotDict[line["Majority Protein ID"].split('-')[0]]]+[line["Position within protein"]]+[line["Amino Acid"]])
                else:
                    if line["Gene name"] in geneNameDict:
                        id = geneNameDict[line["Gene name"]]
                        pos = line["Position within protein"]
                        res = line["Amino Acid"]

                        if id in id_pos_res:
                            id_pos_res[id][pos] = results
                        else:
                            id_pos_res[id] = {pos: res}
                        writer.writerow([geneNameDict[line["Gene name"]]]+[line["Position within protein"]]+[line["Amino Acid"]])

    return id_pos_res

if __name__ == '__main__':
    main()
