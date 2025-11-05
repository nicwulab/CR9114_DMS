# Deep mutational scanning of germline and somatic CR9114 against different influenza hemagglutinin (HA) proteins

## Introduction
This study performs four deep mutational scanning (DMS) experiments to compare the binding affinity landscapes of the germline and somatic CR9114 against H1 stem, H3 HA, and influenza B HA. Raw read files in fastq format from NIH SRA database: [BioProject PRJNA1284397](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1284397)

## Tite-Seq analysis
[Tite-Seq](https://pubmed.ncbi.nlm.nih.gov/28035901/) is used to quantify the binding affinity of individual mutants. Here is the data analysis workflow:
1. Merging forward and reverse reads by [PEAR](https://cme.h-its.org/exelixis/web/software/pear/doc.html) using [./script/merge_reads.py](./script/merge_reads.py).
2. Count variants at the nucleotide level using [./script/CR9114HC_GL_fastq2count.py](./script/CR9114HC_GL_fastq2count.py) for the germline CR9114 DMS experiments and [./script/CR9114HC_WT_fastq2count.py](./script/CR9114HC_WT_fastq2count.py) for the somatic CR9114 DMS experiments.
3. Convert the count data at the nucleotide level to amino acid level using [./script/CR9114HC_GL_count_nuc2aa.py](./script/CR9114HC_GL_count_nuc2aa.py) for the germline CR9114 DMS experiments and [./script/CR9114HC_WT_count_nuc2aa.py](./script/CR9114HC_WT_count_nuc2aa.py) for the somatic CR9114 DMS experiments.
4. Compute the KD:
  * For germline CR9114 vs H1: see [./result/GL_H1/](./result/GL_H1/)
  * For somatic CR9114 vs H1: see [./result/WT_H1/](./result/WT_H1/)
  * For somatic CR9114 vs H3: see [./result/WT_H3/](./result/WT_H3/)
  * For somatic CR9114 vs BHA: see [./result/WT_fluB/](./result/WT_fluB/)
5. Compile the KD results for downstream analysis and visualization using [./script/compile_TiteSeq.py](./script/compile_TiteSeq.py)

## Plotting
Figures in [./graph/](./graph/) are generated as follows:
1. Heatmap showing the effects of each mutation on binding affinity: [./script/plot_KD_heatmap.R](./script/plot_KD_heatmap.R)
2. Compare the binding affinity effects of individual mutations across different DMS experiments: [./script/plot_compare_selection.R](./script/plot_compare_selection.R)
3. Violin plots showing the distribution of binding affinity effects: [./script/script/plot_KD_dist.R](./script/script/plot_KD_dist.R)
4. Scatterplots comparing the binding affinity results with expression level: [./script/plot_QC.R](./script/plot_QC.R)
5. Scatterplots comparing the results from replicates: [./script/plot_bind_vs_exp.R](./script/plot_bind_vs_exp.R)
6. Heatmap showing the frequency of each mutation in the pre-selection input: [./script/plot_heatmap_freq.R](./script/plot_heatmap_freq.R)
7. Plot the effects of somatic hypermutations in known IGHV1-69 HA stem antibodies on CR9114 germline: [./script/plot_SHM_KD.R](./script/plot_SHM_KD.R)
