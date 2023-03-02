import pysam
from Bio import SeqIO

# Инициализация файлов
record = SeqIO.read("referese.fasta", "fasta")
ref_seq = record.seq
samfile = pysam.AlignmentFile("input.sam", "r")

# Инициализация словаря, значений
counts = {}
for i in range(len(ref_seq)):
    counts[i] = {}

start_pos = 25100000
end_pos = 25200000
min_cnt = 2
min_frequency = 0.2

# Обработка, заполнение словаря

for read in samfile.fetch():
    if not read.is_unmapped:
        ref_pos = read.reference_start
        read_seq = read.query_sequence

        for i in range(len(read_seq)):
            if ref_pos+i < start_pos or ref_pos+i > end_pos:
                continue
            if read_seq[i] != ref_seq[ref_pos+i-start_pos]:
                if  read_seq[i] in counts[ref_pos+i-start_pos]:
                    counts[ref_pos+i-start_pos][read_seq[i]] += 1
                else:
                    counts[ref_pos+i-start_pos][read_seq[i]] = 1

# Поиск генетических вариантов
variants = []
for pos in counts:
    total_reads = sum(counts[pos].values())
    if total_reads > 0:
        for alt in counts[pos]:
            alt_reads = counts[pos][alt]
            if alt_reads >= min_cnt and alt_reads/total_reads > min_frequency:
                variants.append((pos, ref_seq[pos], alt, alt_reads))

# Сохранение результатов
with open("variants.vcf", "w") as f:
    f.write("Pos in ref\t|\tRef value\t|\tAlt value\t|\tCount\n")
    for variant in variants:
        f.write(f"{variant[0]+start_pos}\t\t\t\t\t{variant[1]}\t\t\t\t\t\t{variant[2]}\t\t\t\t{variant[3]}\n")
