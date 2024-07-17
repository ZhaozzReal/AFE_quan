# AFE_quan

## Description

DATTSS is developed for Dynamic analyses of Alternative Tandem TSS events from standard RNA-seq data. A key aspect of DATTSS is that it takes advantage of previously reported tandem TSSs. DATTSS incorporates these known TSSs as prior knowledge and then identifies the used tandem TSSs based on ‘change point’ model, which forms a linear regression model to infer the location of proximal tandem TSSs that can best explain the localized changes of read coverage profiles in the first exon regions of transcripts.


## Diagram illuminates the DATTSS algorithm

![image](https://github.com/ZhaozzReal/DATTSS/blob/main/diagram.png)


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

* ```/path/to/allbamfiles.txt``` contains all input filenames (BAM format) of samples of interest. 
The expected format is `,` to separate different files:
```
/path/to/sample1.bam,/path/to/sample2.bam,/path/to/sample3.bam,/path/to/sampleN.bam
```

* ```/path/to/hg38_AFE_annotation.txt``` contains first exon regions from protein-coding genes that are supported by CAGE peaks, which is built on hg38 genome version.*



 
