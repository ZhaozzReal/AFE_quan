from argparse import ArgumentParser,ArgumentTypeError
import HTSeq,collections
from multiprocessing import Pool
parser = ArgumentParser(description = "Quantify first exon expression from standard RNA-seq")
parser.add_argument("-b",dest = 'bamfiles',action = "store",type = str,help = "Input text file with all bamfiles")
parser.add_argument('-anno',dest = 'anno_txt',action = "store",type = str,help = "Input annotation file contains CAGE-supported first exons from protein-coding genes")
parser.add_argument("-p",dest = "processors",action = "store",default = 10,type = int,help = "<INT> Number of processors used [default: 10]")
parser.add_argument("-o",dest = "outfile",action = "store",type = str,help = "Output all first exon expression (count matrix)")
args = parser.parse_args()


def parse_cfgfile(bamfile_txt):
    for line in open(bamfile_txt,"r"):
        lines = line.strip().split("=")
        if lines[0] == "condition1":
            bamfiles_condition1 = lines[1].split(",")
        if lines[0] == "condition2":
            bamfiles_condition2 = lines[1].split(",")
    return bamfiles_condition1,bamfiles_condition2


def Get_junction_num(input_tuple):
    anno,bamfile = input_tuple
    genename,first_exon,strand = anno.split("\t")
    if strand == "+":
        pos = first_exon.split(":")[1].split("-")[1]
    else:
        pos = first_exon.split(":")[1].split("-")[0]
    skip_list = []
    bam_reader = HTSeq.BAM_Reader(bamfile)
    chrom = first_exon.split(':')[0]
    region_fetch = chrom + ":" + str(int(pos) - 50) + "-" + str(int(pos) + 50)
    read_seq = bam_reader.fetch(region = region_fetch)
    for a in read_seq:
        if strand == "+":
            skip_list.extend([int(cigop.ref_iv.start) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
        else:
            skip_list.extend([int(cigop.ref_iv.end) for cigop in a.cigar if cigop.type == "N" and cigop.size >0])
    skip_dict = dict(collections.Counter(skip_list))
    num_list = [ value for key,value in skip_dict.items() if int(pos) - 3 < key < int(pos) + 3]
    if num_list == []:
        num = 0
    else:
        num = max(num_list)
    return num


def AFE_main(anno_txt,all_bamfiles,processors,outfile):
    out = open(outfile,"w")
    out.write("{}\t{}\t{}\n".format("genename","first_exon_region","\t".join(all_bamfiles)))
    for line in open(anno_txt,"r"):
        SYMBOL,first_exon,strand = line.strip().split("\t")
        exon_bamfiles_list = list(zip([line.strip()]*len(all_bamfiles),all_bamfiles))
        pool = Pool(processors)
        num_list = pool.map(Get_junction_num,exon_bamfiles_list)
        out.write("{}\t{}\t{}\n".format(SYMBOL,first_exon + "_" + strand,"\t".join(list(map(str,num_list)))))
    out.close()


bamfiles_condition1,bamfiles_condition2 = parse_cfgfile(args.bamfiles)
all_bamfiles = bamfiles_condition1 + bamfiles_condition2
AFE_main(args.anno_txt,all_bamfiles,args.processors,args.outfile)