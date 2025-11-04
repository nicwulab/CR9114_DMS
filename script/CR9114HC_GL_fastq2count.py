#!/usr/bin/python
import os
import sys
import glob
import operator
from Bio import SeqIO
from collections import Counter

def ProcessMultilib(count_dict, sampleID, Rfile):
  print ("Reading %s" % Rfile)
  records = SeqIO.parse(Rfile,"fastq")
  variants = [] 
  record_count = 0
  for record in records:
    record_count += 1
    Rseq  = record.seq
    Rroi = Rseq
    if ((Rroi[0:23] == "GTGTACTTGCCGCCGCTCAACCA") and (Rroi[-23:] == "GCTTCTACTAAGGGACCTTCCGT") and (len(Rroi[26:-26]) == 360)): # Only include those that have the correct forward primer sequence, correct reverse primer sequence, and the correct number of nucleotides between the flanks
      Rroi = Rroi[26:-26] # Trim forward and reverse flanks
      variants.append(Rroi)
    #if record_count == 1000: break
  count_dict[sampleID] = Counter(variants)
  return count_dict

def Output(count_dict, outfile):
  print ("Compiling results into %s" % outfile)
  outfile = open(outfile,'w')
  sampleIDs = [sampleID for sampleID in sorted(count_dict.keys(),key=lambda x:int(x))]
  muts = list(set(list([mut for sampleID in sampleIDs for mut in count_dict[sampleID].keys()])))
  outfile.write("\t".join(['mut']+sampleIDs)+"\n")
  for mut in muts:
    out = [mut]
    for sampleID in sampleIDs:
      out.append(count_dict[sampleID][mut])
    outfile.write("\t".join(map(str,out))+"\n")
  outfile.close()

def main():
  outfile = 'result/CR9114HC_GL_count_nuc.tsv'
  merged_fastq_files = glob.glob('fastq_merged/*.assembled.fastq')
  count_dict = {}
  for merged_fastq_file in sorted(merged_fastq_files, key= lambda x:int(x.rsplit('/')[1].rsplit('_')[0])):
    sampleID = merged_fastq_file.rsplit('/')[1].rsplit('_')[0]
    count_dict  = ProcessMultilib(count_dict, sampleID, merged_fastq_file)
    print ('processing sampleID: %s' % sampleID)
  Output(count_dict, outfile)

if __name__ == "__main__":
  main()
