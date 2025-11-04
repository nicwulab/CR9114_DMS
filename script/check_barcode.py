#!/usr/bin/python
import os
import sys
import glob
import itertools
from Bio import SeqIO

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "NNK":"X"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def barcode_permutation(nucleotide_3rd):
  barcode_perms = []
  for i in nucleotide_3rd:
    if i == 'K':
      barcode_perms.append(['G','T'])
    else: 
      assert i in ['A','C','T','G']
      barcode_perms.append([i])
  barcode_list = []
  for barcode_perm in itertools.product(*barcode_perms):
    barcode_list.append(''.join(barcode_perm))
  return barcode_list
  
def WT_barcodes(template, roi_length, t_offset):
  WT_seq = template[t_offset:t_offset+roi_length]
  WT_3rd = WT_seq[2::3].upper()
  return WT_3rd

def extract_barcodes(infile, outfile, ref_dict, refID, roi_length, t_offset):
  print ("writing: %s" % outfile)
  records = SeqIO.parse(infile,"fasta")
  outfile = open(outfile,'w')
  barcode_dict = {}
  outfile.write("barcode"+"\t"+"position"+"\n")
  for record in records:
    ID  = record.id
    seq = record.seq
    cas = int(str(ID).rsplit('_')[0].replace('Cassette',''))
    template = ref_dict[refID]
    flank5 = template[t_offset:cas*24-3]
    flank3 = template[t_offset+cas*24:t_offset+roi_length]
    seq = seq[21:21+24]
    seq = flank5+seq+flank3
    nucleotide_3rd = seq[2::3].upper()
    barcode_list = barcode_permutation(nucleotide_3rd)
    barcode_dict[ID] = barcode_list
    #print ("Barcode:"+"\t"+ID+"\t"+','.join(barcode_list))
    for barcode in barcode_list:
      outfile.write(barcode+"\t"+ID+"\n")
  WT_3rd = WT_barcodes(ref_dict[refID], roi_length, t_offset)
  barcode_dict['WT'] = [WT_3rd]
  outfile.write(WT_3rd+"\t"+'WT'+"\n")
  outfile.close()
  return (barcode_dict)

def compare_barcodes(barcode_dict):
  ID_list = list(barcode_dict.keys())
  overlap_count = 0
  for n in range(len(ID_list)):
    for m in range(len(ID_list)):
      if n < m:
        ID1 = ID_list[n]
        ID2 = ID_list[m]
        ID1_barcodes = set(barcode_dict[ID1])
        ID2_barcodes = set(barcode_dict[ID2])
        overlap_barcodes = ID1_barcodes.intersection(ID2_barcodes)
        if len(overlap_barcodes) != 0:
          overlap_count += 1
          print ("Overlapping barcodes between %s and %s:" % (ID1, ID2))
          print (overlap_barcodes)
  if overlap_count == 0: 
    print ("Good job! No overlapping barcodes!")

def read_ref(reffile):
  records = SeqIO.parse(reffile,"fasta")
  ref_dict = {}
  for record in records:
    ID  = str(record.id) 
    seq = str(record.seq)
    ref_dict[ID] = seq
  return ref_dict

def main():
  ref_dict_GL = read_ref('Fasta/CR9114_HC_GL.fa')
  infile_GL  = "primer/CR9114_germline_HC_lib_Fprimer_bc.fa"
  outfile_GL = "data/barcodes_GL.tsv"
  ref_dict_WT = read_ref('Fasta/CR9114_HC_WT.fa')
  infile_WT  = "primer/CR9114_WT_HC_lib_Fprimer_bc.fa"
  outfile_WT = "data/barcodes_WT.tsv"
  roi_length = 360
  t_offset   = 21
  print ('barcode analysis of GL')
  barcode_dict = extract_barcodes(infile_GL, outfile_GL, ref_dict_GL, 'CR9114_HC_GL', roi_length, t_offset)
  compare_barcodes(barcode_dict)
  print ('barcode analysis of WT')
  barcode_dict = extract_barcodes(infile_WT, outfile_WT, ref_dict_WT, 'CR9114_HC_WT', roi_length, t_offset)
  compare_barcodes(barcode_dict)

if __name__ == "__main__":
  main()
