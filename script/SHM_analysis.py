#/usr/bin/python
import sys
import pandas as pd
from collections import defaultdict

def KD_to_dict(file_KD):
  df = pd.read_csv(file_KD)
  KD_dict = df.set_index('mut_kabat')['minus_delta_log_Kd'].to_dict()
  return KD_dict

def get_num_pos(SHM):
  pos = SHM[1:-1]
  try:
    return int(pos)
  except:
    return int(pos[0:-1])

def assign_KD_to_ab(outfile, file_ab, KD_dict, CR9114_SHM_list):
  df = pd.read_csv(file_ab)
  print ('writing: %s' % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['ab_name', 'SHM', 'in_CR9114', 'minus_delta_log_KD (GL)'])+"\n")
  for index, info_dict in df.iterrows():
    ab_name  = info_dict['Name']
    if info_dict['Binds to'] != 'HA:Stem': continue
    if str(info_dict['Heavy_V_gene']) == 'nan': continue
    if info_dict['Heavy_V_gene'].rsplit('*')[0] != 'IGHV1-69': continue
    if str(info_dict['Heavy_shm']) == 'nan': continue
    ab_name  = info_dict['Name']
    SHM_list = info_dict['Heavy_shm'].rsplit(',')
    for SHM in SHM_list:
      SHM = SHM[0]+SHM[1:-1].lower()+SHM[-1]
      if get_num_pos(SHM) <= 91:
        if SHM in KD_dict.keys():
          KD = KD_dict[SHM]
          if str(KD) != 'nan':
            print (SHM[1::])
            in_CR9114 = 'yes' if SHM[1::] in CR9114_SHM_list else 'no'
            outfile.write("\t".join(map(str,[ab_name, SHM, in_CR9114, KD]))+"\n")
  outfile.close()

def main():
  outfile = 'result/IGHV1-69_HAstem_Ab_SHM_KD.tsv'
  file_KD = 'result/HC_GL_KD_table_summary_kabat.csv'
  file_ab = 'data/IGHV1-69_HAstem_Abs_SHM.csv'
  CR9114_SHM_list = ['24S','29S','30N','31N','46D','52S','56S','57T','58A','70S','73I','74F','75S','76N','82aN','83T','91F','100bS']
  KD_dict = KD_to_dict(file_KD)
  assign_KD_to_ab(outfile, file_ab, KD_dict, CR9114_SHM_list)

if __name__ == "__main__":
  main()
