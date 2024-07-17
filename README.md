# AFE_quan

## Description

DATTSS is developed for Dynamic analyses of Alternative Tandem TSS events from standard RNA-seq data. A key aspect of DATTSS is that it takes advantage of previously reported tandem TSSs. DATTSS incorporates these known TSSs as prior knowledge and then identifies the used tandem TSSs based on ‘change point’ model, which forms a linear regression model to infer the location of proximal tandem TSSs that can best explain the localized changes of read coverage profiles in the first exon regions of transcripts.


## Diagram illuminates the DATTSS algorithm


## Installation

DATTSS is built on Python, which requires packages ```HTSeq```, ```collections```, ```multiprocessing``` and ```argparse```.

Clone the lastest development version of DATTSS and change directory:

```
  git clone https://github.com/ZhaozzReal/AFE_quan.git
  cd AFE_quan
```

## Usage

```
python AFE_quan.py -b /path/to/allbamfiles.txt -anno /path/to/hg38_AFE_annotation.txt -p 10 -o /path/to/AFE_quan_output.txt
```

 ```/path/to/allbamfiles.txt``` contains all input filenames (BAM format) of samples of interest. 
The expected format is `,` to separate different files:
```
/path/to/sample1.bam,/path/to/sample2.bam,/path/to/sample3.bam,/path/to/sampleN.bam
```

 ```/path/to/hg38_AFE_annotation.txt``` contains non-redundant first exon regions from protein-coding genes that are supported by CAGE peaks, which is built on hg38 genome version.



## Output

In ```/path/to/AFE_quan_output.txt```, each row corresponds to one first exon region and its related features.

The explanation of each column is as follows:
 
 * genename: the HUGO gene symbol
 * first_exon_region：the genomic region and strand information of first exon region ranging from transcription start site to first exon end
 * sample1.bam：the expression level (junction count) of given FE in sample1
 * sample2.bam：the expression level (junction count) of given FE in sample2
 * sampleN.bam：the expression level (junction count) of given FE in sampleN


<br/>
<br/>

## Compare alternative first exon usage between conditions

***Step1: Detect and quantify alternative first exon events***


```
python AFE_compare.py -b /path/to/allbamfiles.txt -anno /path/to/hg38_AFE_annotation.txt -p 10 -o /path/to/AFE_quan_output.txt
```

allbamfiles.txt contains all filename of bamfile between two conditions, as shown below:

```
condition1=/path/to/ctrl1.bam,/path/to/ctrl2.bam 
condition2=/path/to/case1.bam,/path/to/case2.bam
```


***Step2: Infer significantly dysregulated alternative first exons between conditions using DEXSeq model***

We utilizes DEXSeq, the model for differential exon usage analysis based on standard RNA-seq data, to detect differential usage of alternative first exon. This statistical framework could account for biological variability between replicates and is robust to changes in isoform abundance between conditions.


```
Rscript Infer_DU_AFE.R -b /path/to/allbamfiles.txt -c /path/to/AFE_quan_output.txt -d /path/to/DEXSeq_count/ -o /path/to/AFE_DU.txt
```
 ```/path/to/DEXSeq_count/``` is a directory to be created, which is used to store count files for DEXSeq model. 

Final results will be saved in the file ```AFE_DU.txt```.




 
