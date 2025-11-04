#!/usr/bin/python
import os
import sys
import operator
import itertools
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

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

def seq_to_kabat(kabat_file):
  kabat_nums = open(kabat_file,'r').readlines()[0].rstrip().rsplit('_identity,')[-1].rsplit(',')
  seq2kabat_dict = {'WT':'WT'}
  for seq_num, kabat_num in zip(range(len(kabat_nums)),kabat_nums):
    seq2kabat_dict[str(seq_num+1)] = str(kabat_num)
  return seq2kabat_dict

def get_kabat(mut, seq2kabat_dict):
  if mut == 'WT':
    return 'WT'
  else:
    return mut[0] + seq2kabat_dict.get(mut[1:-1], 'NA') + mut[-1]

def get_mut_class(mut):
  if mut == 'WT':
    return 'WT'
  elif mut[-1] == '_':
    return 'nonsense'
  elif mut[0] == mut[-1]:
    return 'silent'
  else:
    return 'missense'

def add_kabat(df, seq2kabat_dict):
  df['mut_kabat'] = df['mut'].apply(lambda mut: get_kabat(mut, seq2kabat_dict))
  return (df)
 
def compute_delta_Kd(df):
  mean_log_Kd_silent = df[df['mut_class'] == 'silent']['minus_log_Kd'].mean()
  mean_log_Kd_silent_rep1 = df[df['mut_class'] == 'silent']['minus_log_Kd_rep1'].mean()
  mean_log_Kd_silent_rep2 = df[df['mut_class'] == 'silent']['minus_log_Kd_rep2'].mean()
  df['minus_delta_log_Kd'] = (df['minus_log_Kd']) - mean_log_Kd_silent
  df['minus_delta_log_Kd_rep1'] = (df['minus_log_Kd_rep1']) - mean_log_Kd_silent_rep1
  df['minus_delta_log_Kd_rep2'] = (df['minus_log_Kd_rep2']) - mean_log_Kd_silent_rep2
  print ("mean of silent: %f" % mean_log_Kd_silent)
  print ('range of log10_Kd (before cutoff): %f to %f' % (df['minus_delta_log_Kd'].min(), df['minus_delta_log_Kd'].max()))
  return (df)

def compute_expression_score(countfile, exp_bins, presort, selection, MFI_dict):
  df = pd.read_csv(countfile, sep="\t")
  for rep in ['rep1', 'rep2']:
    df['ipt_freq_'+rep]  = df[presort[selection][rep]]/df[presort[selection][rep]].sum()
    df['exp_score_'+rep] = (df[exp_bins[selection][rep+'_bin1']]/df[exp_bins[selection][rep+'_bin1']].sum()*MFI_dict[selection][rep+'_bin1'] +
			    df[exp_bins[selection][rep+'_bin2']]/df[exp_bins[selection][rep+'_bin2']].sum()*MFI_dict[selection][rep+'_bin2'] +
			    df[exp_bins[selection][rep+'_bin3']]/df[exp_bins[selection][rep+'_bin3']].sum()*MFI_dict[selection][rep+'_bin3'] +
			    df[exp_bins[selection][rep+'_bin4']]/df[exp_bins[selection][rep+'_bin4']].sum()*MFI_dict[selection][rep+'_bin4'])/ \
			   (df[exp_bins[selection][rep+'_bin1']]/df[exp_bins[selection][rep+'_bin1']].sum() +
			    df[exp_bins[selection][rep+'_bin2']]/df[exp_bins[selection][rep+'_bin2']].sum() +
			    df[exp_bins[selection][rep+'_bin3']]/df[exp_bins[selection][rep+'_bin3']].sum() +
			    df[exp_bins[selection][rep+'_bin4']]/df[exp_bins[selection][rep+'_bin4']].sum())
    df = normalization(df, 'exp_score_'+rep, 'exp_score_'+rep)
  df['exp_score']  = (df['exp_score_rep1'] + df['exp_score_rep2']) / 2
  df['ipt_freq']  = (df['ipt_freq_rep1'] + df['ipt_freq_rep2']) / 2
  return (df)

def get_region(mut_kabat):
  if mut_kabat == 'WT': return 'WT'
  else:
    pos = mut_kabat[1:-1]
    try: 
      pos = int(pos)
    except:
      pos = int(pos[0:-1])
    if pos >= 1 and pos <= 25: return 'FR1'
    elif pos >= 26 and pos <= 32: return 'CDR1'
    elif pos >= 33 and pos <= 51: return 'FR2'
    elif pos >= 52 and pos <= 56: return 'CDR2'
    elif pos >= 57 and pos <= 70: return 'FR3'
    elif pos >= 71 and pos <= 78: return 'DE Loop'
    elif pos >= 79 and pos <= 94: return 'FR3'
    elif pos >= 95 and pos <= 102: return 'CDR3'
    elif pos >= 103 and pos <= 112: return 'FR4'
    else: return 'Other'

def normalization(df, param, new_param):
  w_summary = df.groupby('mut_class')[param].mean()
  w_summary = w_summary.reset_index()
  w_silent   = float(w_summary.loc[w_summary['mut_class']=='silent'][param].iloc[0])
  w_nonsense = float(w_summary.loc[w_summary['mut_class']=='nonsense'][param].iloc[0])
  df[new_param] = (df[param]-w_nonsense)/(w_silent-w_nonsense)
  return (df)
 
def wrapper(dict_mut_dist, infile_rep1, infile_rep2, outfile, countfile, kabat_file, selection):
  exp_bins = {'GL':  {'rep1_bin1':'33','rep1_bin2':'34','rep1_bin3':'35','rep1_bin4':'36',
                      'rep2_bin1':'73','rep2_bin2':'74','rep2_bin3':'75','rep2_bin4':'76'},
              'H1':  {'rep1_bin1':'17-A','rep1_bin2':'17-B','rep1_bin3':'17-C','rep1_bin4':'17-D',
                      'rep2_bin1':'22-A','rep2_bin2':'22-B','rep2_bin3':'22-C','rep2_bin4':'22-D'},
              'H3':  {'rep1_bin1':'15-A','rep1_bin2':'15-B','rep1_bin3':'15-C','rep1_bin4':'15-D',
                      'rep2_bin1':'16-A','rep2_bin2':'16-B','rep2_bin3':'16-C','rep2_bin4':'16-D'},
              'fluB':{'rep1_bin1':'65','rep1_bin2':'66','rep1_bin3':'67','rep1_bin4':'68',
                      'rep2_bin1':'69','rep2_bin2':'70','rep2_bin3':'71','rep2_bin4':'72'}}
  presort  = {'GL':    {'rep1':'81', 'rep2':'81'},
              'H1':    {'rep1':'17-pre', 'rep2':'22-pre'},
              'H3':    {'rep1':'15-pre', 'rep2':'16-pre'},
              'fluB':  {'rep1':'73', 'rep2':'74'}}
  MFI_dict = {'GL':   {'rep1_bin1':29, 'rep1_bin2':468, 'rep1_bin3':1843, 'rep1_bin4':5708,
                       'rep2_bin1':29, 'rep2_bin2':487, 'rep2_bin3':1858, 'rep2_bin4':7130},
              'H1':   {'rep1_bin1':7, 'rep1_bin2':146, 'rep1_bin3':532, 'rep1_bin4':2049,
                       'rep2_bin1':6, 'rep2_bin2':146, 'rep2_bin3':541, 'rep2_bin4':2310},
              'H3':   {'rep1_bin1':10, 'rep1_bin2':146, 'rep1_bin3':570, 'rep1_bin4':2234,
                       'rep2_bin1':13, 'rep2_bin2':137, 'rep2_bin3':525, 'rep2_bin4':1944},
              'fluB': {'rep1_bin1':7, 'rep1_bin2':143, 'rep1_bin3':530, 'rep1_bin4':2158,
                       'rep2_bin1':7, 'rep2_bin2':141, 'rep2_bin3':571, 'rep2_bin4':2807}} 
  seq2kabat_dict = seq_to_kabat(kabat_file)
  df_Kd_rep1  = pd.read_csv(infile_rep1)
  df_Kd_rep2  = pd.read_csv(infile_rep2)[['mut','Kd']]
  df_Kd_rep1.rename(columns={'Kd': 'Kd_rep1'}, inplace=True)
  df_Kd_rep2.rename(columns={'Kd': 'Kd_rep2'}, inplace=True)
  df_exp      = compute_expression_score(countfile, exp_bins, presort, selection, MFI_dict)
  df = df_Kd_rep1.merge(df_Kd_rep2, on='mut', how='outer') \
                 .merge(df_exp, on='mut', how='outer')
  df['minus_log_Kd_rep1'] = -np.log10(df['Kd_rep1'])
  df['minus_log_Kd_rep2'] = -np.log10(df['Kd_rep2'])
  df['minus_log_Kd'] = (df['minus_log_Kd_rep1'] + df['minus_log_Kd_rep2']) / 2 
  df = normalization(df, 'minus_log_Kd_rep1', 'bind_score_rep1')
  df = normalization(df, 'minus_log_Kd_rep2', 'bind_score_rep2')
  df = normalization(df, 'minus_log_Kd', 'bind_score')
  df = add_kabat(df, seq2kabat_dict)
  df = compute_delta_Kd(df)
  df.rename(columns={'p.value': 'p_value'}, inplace=True)
  df['mut_class'] = df['mut'].apply(get_mut_class)
  df['region'] = df['mut_kabat'].apply(get_region)
  if selection in ['H1','H3','fluB']:
    df['codon_dist'] = df['mut'].map(dict_mut_dist['CR9114_WT'])
    df_out = df[['mut', 'mut_kabat', 'codon_dist', 'region', 'mut_class', 'bind_score', 'bind_score_rep1', 'bind_score_rep2',
                 'minus_delta_log_Kd', 'minus_delta_log_Kd_rep1', 'minus_delta_log_Kd_rep2',
                 'minus_log_Kd', 'minus_log_Kd_rep1', 'minus_log_Kd_rep2', 'ipt_freq_rep1','ipt_freq_rep2','ipt_freq']]
  if selection in ['GL']:
    df['codon_dist'] = df['mut'].map(dict_mut_dist['CR9114_GL'])
    df_out = df[['mut', 'mut_kabat', 'codon_dist', 'region', 'mut_class', 'bind_score', 'bind_score_rep1', 'bind_score_rep2',
                 'minus_delta_log_Kd', 'minus_delta_log_Kd_rep1', 'minus_delta_log_Kd_rep2', 
                 'minus_log_Kd', 'minus_log_Kd_rep1', 'minus_log_Kd_rep2','ipt_freq', 'exp_score','exp_score_rep1','exp_score_rep2']]
  print ("writing: %s" % outfile)
  df_out.to_csv(outfile, index=False)
  return (df)

def write_expression_score_file(outfile, df_WT_H1, df_WT_H3, df_WT_fluB):
  df_WT_H1.rename(columns={'exp_score_rep1': 'exp_score_rep1'}, inplace=True)
  df_WT_H1.rename(columns={'exp_score_rep2': 'exp_score_rep2'}, inplace=True)
  df_WT_H1 = df_WT_H1[['mut','mut_class','exp_score_rep1', 'exp_score_rep2']]
  df_WT_H3.rename(columns={'exp_score_rep1': 'exp_score_rep3'}, inplace=True)
  df_WT_H3.rename(columns={'exp_score_rep2': 'exp_score_rep4'}, inplace=True)
  df_WT_H3 = df_WT_H3[['mut','exp_score_rep3', 'exp_score_rep4']]
  df_WT_fluB.rename(columns={'exp_score_rep1': 'exp_score_rep5'}, inplace=True)
  df_WT_fluB.rename(columns={'exp_score_rep2': 'exp_score_rep6'}, inplace=True)
  df_WT_fluB = df_WT_fluB[['mut','exp_score_rep5', 'exp_score_rep6']]
  df = df_WT_H1.merge(df_WT_H3,   on='mut', how='outer') \
               .merge(df_WT_fluB, on='mut', how='outer')
  #df['exp_score'] = df[['exp_score_rep1','exp_score_rep2','exp_score_rep3','exp_score_rep4','exp_score_rep5','exp_score_rep6']].median(axis=1)
  df['exp_score'] = df['exp_score_rep1']
  df['mut_class'] = df['mut'].apply(get_mut_class)
  df = normalization(df, 'exp_score', 'exp_score')
  print ('writing: %s' % outfile)
  df.to_csv(outfile, index=False)

def seq_to_codon_dist(infile):
  AA_codon = {
     'C': ['TGT', 'TGC'],
     'D': ['GAT', 'GAC'],
     'S': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
     'Q': ['CAA', 'CAG'],
     'M': ['ATG'],
     'N': ['AAC', 'AAT'],
     'P': ['CCT', 'CCG', 'CCA', 'CCC'],
     'K': ['AAG', 'AAA'],
     '_': ['TAG', 'TGA', 'TAA'],
     'T': ['ACC', 'ACA', 'ACG', 'ACT'],
     'F': ['TTT', 'TTC'],
     'A': ['GCA', 'GCC', 'GCG', 'GCT'],
     'G': ['GGT', 'GGG', 'GGA', 'GGC'],
     'I': ['ATC', 'ATA', 'ATT'],
     'L': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
     'H': ['CAT', 'CAC'],
     'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
     'W': ['TGG'],
     'V': ['GTA', 'GTC', 'GTG', 'GTT'],
     'E': ['GAG', 'GAA'],
     'Y': ['TAT', 'TAC']}

  dict_mut_dist = {}
  for record in SeqIO.parse(infile, 'fasta'):
    ID = str(record.id)  
    seq = str(record.seq)  
    dict_mut_dist[ID] = defaultdict(int)
    i = 0
    while i < len(seq):
      WT_codon = seq[i:i+3].upper()
      WT_aa = translation(WT_codon)
      aa_pos =  str(int(i/3+1))
      for AA in AA_codon.keys():
        dict_mut_dist[ID][WT_aa+aa_pos+AA] = min([hamming(WT_codon, codon) for codon in AA_codon[AA]])
      i = i + 3
  return (dict_mut_dist)

def main():
  dict_mut_dist = seq_to_codon_dist('Fasta/CR9114_original_seq.fa')
  df_GL_H1   = wrapper(dict_mut_dist, 'result/GL_H1/result/HC_KD_table_rep1.csv','result/GL_H1/result/HC_KD_table_rep2.csv',
                       'result/HC_GL_KD_table_summary_kabat.csv','result/CR9114HC_GL_count_aa.tsv',"result/CR9114HC_GL_kabat_H.csv", 'GL')
  df_WT_H1   = wrapper(dict_mut_dist, 'result/WT_H1/result/HC_KD_table_rep1.csv','result/WT_H1/result/HC_KD_table_rep2.csv',
                       'result/HC_WT_H1_KD_table_summary_kabat.csv','result/CR9114HC_WT_count_aa_H1.tsv',"result/CR9114HC_WT_kabat_H.csv", 'H1')
  df_WT_H3   = wrapper(dict_mut_dist, 'result/WT_H3/result/HC_KD_table_rep1.csv','result/WT_H3/result/HC_KD_table_rep2.csv',
                       'result/HC_WT_H3_KD_table_summary_kabat.csv','result/CR9114HC_WT_count_aa_H3.tsv',"result/CR9114HC_WT_kabat_H.csv", 'H3')
  df_WT_fluB = wrapper(dict_mut_dist, 'result/WT_fluB/result/HC_KD_table_rep1.csv','result/WT_fluB/result/HC_KD_table_rep2.csv',
                       'result/HC_WT_fluB_KD_table_summary_kabat.csv','result/CR9114HC_WT_count_aa_fluB.tsv',"result/CR9114HC_WT_kabat_H.csv", 'fluB')
  write_expression_score_file('result/HC_WT_expression_score_kabat.csv', df_WT_H1, df_WT_H3, df_WT_fluB)

if __name__ == "__main__":
  main()
