import sys
from Bio import SeqIO

reference = sys.argv[1]
print(reference)
def calculate_genome_length(fasta_file):
  total_length = 0
  with open(fasta_file, 'r') as file:
    for record in SeqIO.parse(file, 'fasta'):
      total_length += len(record.seq)
  with open("genofreq/size_genome.txt", 'w') as file:
    #ifile.write(total_length)  
    print(file.write(str(total_length)))
    #with open("C:/Users/me/Downloads/Documents/lala",mode="w")as f:
    #    print(f.write(a))


#print(total_length, file=open("../genofreq/size_genome.txt", "a"))
calculate_genome_length(reference)



