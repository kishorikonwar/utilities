from Bio import SeqIO
from Bio.Seq import Seq
import sys
import argparse
import re
import random
import string


class GTFLine:
    def __init__(self, line):
       try:
           fields = [x.strip() for x in line.strip().split('\t')] 
           self.chr = fields[0] 
           self.type = fields[2] 
           self.start = int(fields[3])
           self.end = int(fields[4])
           self.strand = fields[6]
           subfields = [re.sub("\"","", x.strip()).split() for x in fields[8].strip().split(';')] 

           self.gene_id = None
           for pair in subfields:
              if pair[0]=='gene_id':
                self.gene_id = pair[1]
           #self.gene_name = re.sub("\"","",  subfields[3].split()[1])
           #print(fields[8].strip())
       except:
           print("ERROR: in line - " + line)


class GTFReader:
    commPATT=re.compile(r'#')
    def __init__(self, gtf_file):
        try:
           self.fpin = open(gtf_file, 'r')
        except FileNotFoundError: 
           raise

    def __iter__(self):
        return self

    def __next__(self):
        line = self.fpin.readline()
        while self.commPATT.search(line):
           line = self.fpin.readline()

        if not line:
            raise StopIteration
        else:
            gtf_line = GTFLine(line) 
            return gtf_line

        

def read_gtf_file(annotation_file):
   gtfs = []
   gtfreader = GTFReader(annotation_file)
   for gtf in gtfreader:
       gtfs.append(gtf)

   return gtfs


def read_and_replace(args):
   gtfs = read_gtf_file(args.annot)
   numgtfs  = len(gtfs)

   refrecords = list(SeqIO.parse(args.ref, "fasta"))

   seq = str(refrecords[0].seq)
   patt = re.compile(r'^[ATCG]')

   j =0
   i = 10

   output_handle= open(args.output, "w") 
   newrecords = []
   for record in  SeqIO.parse(args.fastq,"fastq"):
      length = len(record.seq)
      
      newseq = seq[gtfs[j%numgtfs].start: gtfs[j%numgtfs].start+length]
      record.seq = Seq(newseq)
      j = j + 1
      SeqIO.write(record, output_handle, "fastq")

INST_NAME="EAS139"
RUN_ID='1'
FLOW_CELL_ID='FC706VJ'
FLOW_CELL_LANE='1'
FLOW_CELL_TILENO='1'
X_COOR='1'
Y_COOR='100'
MATE='1'
FILTERED='Y'
BITS='0'
ISEQ='0'

SCORE='FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'

INDSEQ='GGCAATGG'

def create_read_name(y, m):
    name ='@{}:{}:{}:{}:{}:{}:{}  {}:{}:{}:{}'.format(INST_NAME, RUN_ID, FLOW_CELL_ID, 
                                              FLOW_CELL_LANE, FLOW_CELL_TILENO,
                                              X_COOR, y, m, FILTERED, BITS, ISEQ)

    return name


def randnuc(n):
    return ''.join( [random.choice('ATCG') for x in range(n) ])


UMISEQ='ATGC'
def UMI(i, n):
    
    umi_ind = []
    v = i
    while v!=0 or len(umi_ind)<n:
      r = v%4
      v = (v -r)/4
      umi_ind.append(int(r))

    umi = ''.join( [ UMISEQ[x] for x in umi_ind] )
    return str(umi)
    
    
translator = str.maketrans("ATCG", "TAGC")


def rev_complement(seqstr):
    i = 0
    j = len(seqstr) -1
   
    seq = list(seqstr)  
    while i < j:
       u = seq[i] 
       seq[i] = seq[j] 
       seq[j] = u 
       i += 1
       j -= 1

    revstr = ''.join(seq)
    revcomp = revstr.translate(translator)
    return revcomp


def get_read2(reference, start, end, strand): 
    if strand =="+":
       return reference[start:end]
    return rev_complement(reference[start:end])

def main(args):
    refrecords = list(SeqIO.parse(args.refgenome, "fasta"))
    refseq = str(refrecords[0].seq)
    barcodes = read_barcodes(args.wlist)


    gtfs = read_gtf_file(args.refgtf)
    gene_locs ={}
    for gtf in gtfs:
     if gtf.type=='CDS' and abs(int(gtf.start)- int(gtf.end))>=91 :
        if not gtf.gene_id in gene_locs: 
          gene_locs[gtf.gene_id] = [ int(gtf.start), int(gtf.end), gtf.strand ]

    # read the reference fasta
    refrecords = list(SeqIO.parse(args.refgenome, "fasta"))
    seq = str(refrecords[0].seq)
    
    # read the gene distribution file 
    genes, weights = read_genes_count_dist(args.gene_count_dist)


    # sample the reads
    gene_samples = random.choices(population=genes, weights=weights, k=args.nreads)


    with open(args.output_prefix + "_R1_001.fastq", 'w') as fR1, open(args.output_prefix + "_R2_001.fastq", 'w')  as fR2, open(args.output_prefix + "_I1_001.fastq", 'w') as fI1 :
        for i, gene_id in enumerate(gene_samples): 
          fR1.write(create_read_name(i,'1') +'\n')
          fR1.write(barcodes[i]+randnuc(10) +'\n')
          fR1.write('+' +'\n')
          fR1.write(SCORE[0:26] +'\n')
    
          fR2.write(create_read_name(i,'3') +'\n')
          fR2.write(get_read2(refseq, gene_locs[gene_id][0], gene_locs[gene_id][0] + 91, gene_locs[gene_id][2]) + '\n') 
          fR2.write('+' + '\n')
          fR2.write(SCORE[0:91] + '\n')
    
          fI1.write(create_read_name(i,'2') + '\n')
          fI1.write('ACATTACT' + '\n')
          fI1.write('+' + '\n')
          fI1.write(SCORE[0:8] + '\n')

   
def read_genes_count_dist(gene_count_dist):
    genes =[]
    weights =[]
    with open(gene_count_dist, 'r') as fin:
       for x in fin.readlines():
         x = x.split('\t')
         genes.append(x[0])
         weights.append(float(x[1]))
    return genes, weights


def read_barcodes(barcodefile):
    barcodes = []
    with open(barcodefile, 'r') as fin:
       barcodes = [ x.strip() for x in fin.readlines() ]

    return barcodes


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Replaces the sequence with one from the reference")
    parser.add_argument('--ref-genome', 
                      dest='refgenome', 
                      help='reference sequence file') 

    parser.add_argument('--ref-gtf', 
                      dest='refgtf', 
                      help='annotation file .gtf') 


    parser.add_argument('--white-list', 
                      dest='wlist', 
                      help='white list') 

    parser.add_argument('--nreads', 
                      dest='nreads', type=int, default=0,
                      help='number of reads') 

    parser.add_argument('--gene-count-dist', 
                      dest='gene_count_dist', default=None,
                      help='gene count distribution') 

    parser.add_argument('--op', 
                      dest='output_prefix', 
                      help='output prefix') 

    args = parser.parse_args()
    #read_and_replace(args)
    main(args)

