'''Accepts fasta files containing amino acid sequence, extracts characters as
defined in an array (.csv file) into a new string and creates new protein Sequences
with respective extract strings. Adds these to a list of extracted Sequence objects and
outputs the list (representing flu antigenic sites) into an amino acid fasta file.'''

'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory,Sept 2017'''

import sys,string,os, time, Bio
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

seqList = [] #list of aa sequence objects to parse for oligo sequences
indexArray = [] # .csv list of aa's corresponding to antigenic site positions
extractedSeqList = [] #list of extracted antigenic sites extracted from seqList
inFileHandle1 = sys.argv[1] #batch fasta file with sequences to be parsed
inFileHandle2 = sys.argv[2] # .csv file containing positions of aa's to extract
outFileHandle = sys.argv[3] #user-specified name for output file of extracted aa seq's
outFile= open(outFileHandle,'w') #open a writable, appendable output file
localtime = time.asctime(time.localtime(time.time())) #date and time of analysis

def extract_aa_from_sequence(record):
    """Extract specific amino acids from SeqRecord, create new SeqRecord and append to list."""
    original_sequence = str(record.seq) #pull out the SeqRecord's Seq object and ToString it
    new_sequence = "" #set variable to empty
    new_id = record.id #store the same sequence id as the original sequence
    #iterate over each position in index array, extract corresponding aa and add to string
    for pos in indexArray:
        char = original_sequence[pos-1] #aa positions must be zero indexed
        new_sequence = new_sequence + char
    rec = SeqRecord(Seq(new_sequence,IUPAC.protein), id = record.id, name = "", description = "")
    extractedSeqList.append(rec) #add new SeqRecord object to the list

with open (inFileHandle2,'r') as inFile2:
    '''Open csv file containing amino acid positions to extract and add to list.'''
    #read items separated by comma's to position list
    positionList = ""   
    for line in inFile2:
        #remove whitespace from the end of each line
        strippedLine = line.rstrip()
        #split the line at commas and assigned the returned list as indexArray
        positionList = strippedLine.split(',')
    #Convert string items in positionList from strings to int and add to indexArray
    for item in positionList:
        indexArray.append(int(item))
    #print number of amino acids to extract and array to console as user check
    print "Amino Acid positions to extract: %i " %(len(indexArray))
    print indexArray
        
with open(inFileHandle1,'r') as inFile:
    '''Open fasta of amino acid sequences to parse, uppercase and add to protein Sequence list.'''
    for record in SeqIO.parse(inFile, "fasta", alphabet=IUPAC.protein):
        record = record.upper()
        seqList.append(record) #add Seq to list of Sequences
    #print number of sequences to be process as user check
    print "\n%i flu sequences will be extracted for antigenic sites..." % len(seqList)
    #parse each target sequence object
    for record in seqList:
        extract_aa_from_sequence(record)

#print original and extracted sequence
for x in range(0, len(seqList)):
    print "Original %s: %i amino acids,\tExtracted: %i" % (seqList[x].id,len(seqList[x]),len(extractedSeqList[x]))

#write fasta file of extracted antigenic sites
SeqIO.write(extractedSeqList,outFile,"fasta")

print("\n%i Sequences Extracted to Output file: %s"  % ((len(extractedSeqList),outFileHandle)))
inFile.close()
inFile2.close()
outFile.close()

