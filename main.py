#!pip install pysam
#!pip install Bio
import pysam
from Bio import SeqIO
# Инициализация файлов
record = SeqIO.read("/content/drive/MyDrive/Colab Notebooks/referese.fasta", "fasta")
ref_seq = str(record.seq)
samfile = pysam.AlignmentFile("/content/drive/MyDrive/Colab Notebooks/input.sam", "r")
# Инициализация словаря, значений
start_pos = 25100000
end_pos = 25200000
min_cnt = 2
min_frequency = 0.2

counts = {}
for i in range(start_pos, end_pos + 1):
    counts[i] = {}


ref_seq_dict = {}
for i in range(start_pos, end_pos + 1):
  ref_seq_dict[i] = ref_seq[i-start_pos]

# Обработка, заполнение словаря
for read in samfile.fetch():
    if not read.is_unmapped:
        ref_pos = read.reference_start + 1 
        read_seq = read.query_sequence
        for i in range(len(read_seq)):
            if ref_pos+i < start_pos or ref_pos+i > end_pos:
                continue
            if read_seq[i] != ref_seq_dict[ref_pos + i]:
                if  read_seq[i] in counts[ref_pos + i]:
                    counts[ref_pos + i][read_seq[i]] += 1
                else:
                    counts[ref_pos + i][read_seq[i]] = 1

# Поиск генетических вариантов
variants = []
for pos in counts:
    total_reads = sum(counts[pos].values())
    if total_reads > 0:
        for alt in counts[pos]:
            alt_reads = counts[pos][alt]
            if alt_reads >= min_cnt and alt_reads/total_reads > min_frequency:
                variants.append((pos, ref_seq_dict[pos], alt, alt_reads))
              
with open('variants.vcf', 'w') as vcf_file:
    chrom = 'chr1'
    vcf_file.write('##fileformat=VCFv4.3\n')
    vcf_file.write('##INFO=<ID=COUNT,Number=1,Type=Integer,Description="Count">\n')
    vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for variant in variants:
      pos = variant[0] 
      ref = variant[1]
      alt = variant[2]
      count = variant[3]
      vcf_file.write('{}\t{}\t.\t{}\t{}\t.\t.\tCOUNT={}\n'.format(chrom, pos, ref, alt, count))

vcf_file.close()

import pandas as pd
row_n = 1
dict_pd = {}

colNames = ["CHROM", "POS", "REF", "ALT"]
for variant in variants:
  pos = variant[0]
  ref = variant[1]
  alt = variant[2]
  dict_pd[row_n] = ['chr1'] + [pos, ref, alt]
  row_n += 1
df = pd.DataFrame.from_dict(dict_pd, orient='index', columns=colNames)
df.to_csv("output.csv", sep='\t')
