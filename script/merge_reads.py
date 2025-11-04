#!/usr/bin/python
import os
import glob

def main():
    R1_files  = glob.glob('fastq/*_R1_*')
    outfolder = 'fastq_merged/'
    for R1_file in R1_files:
        R2_file = R1_file.replace('_R1_','_R2_')
        outfile_prefix = outfolder+'/'+R1_file.rsplit('_R1_')[0].rsplit('/')[1].replace('.fastq.gz','')
        assert(R1_file != R2_file)
        os.system('/Users/wchnicholas/Softwares/PEAR/src/pear -j 7 -f '+R1_file+' -r '+R2_file+' -o '+outfile_prefix)

if __name__ == "__main__":
    main()
