import argparse
import sys
import os

START_CODON = "ATG"
STOP_CODONS = ["TAG", "TAA", "TGA"]

DNA_MAP = {
  "A": "T",
  "T": "A",
  "G": "C",
  "C": "G"
}
### importing csv as 3 lists
codon = []
letter3 = []
letter1 = []
infile = open("codon_table.csv", encoding="utf8")
infile.readline()
i=1
for line in infile:
    line = line.strip()
    data = line.split(',')
    codon.append(data[0])
    #print(data)
    letter3.append(data[1])
    letter1.append(data[2])

### Making dictionary mapping codon to amino acid single letter code
prot_dict = {}
for i in range(len(codon)):
  prot_dict[codon[i]] = letter1[i]
###

def check_existance(file):
  if not os.path.isfile(file):
    print("File doesn't exist or the path is incorrect: {}.".format(file))
    sys.exit("Check {} exists.".format(file))

# returns complementary dna strand
def find_complement(dna):
  complement = ""
  for bp in dna:
    complement += DNA_MAP[bp]
  
  return complement[::-1] # need to reverse this

# finds all genes in a sequence of dna (starts with "ATG" and ends with "TAG", "TAA", or "TGA")
def find_genes(_dna):
  gene_start = None
  gene_end = None
  genes = []

  readframe_1 = 0
  readframe_2 = 1
  readframe_3 = 2

  # read each group of 3 bases starting from beginning of a reading frame
  # start reading frame 1
  for right in range(readframe_1 + 3, len(_dna) + 1, 3):
    codon = _dna[right-3:right]
    
    if codon == START_CODON:
      gene_start = right-3
    
    if codon in STOP_CODONS and gene_start is not None:
      gene_stop = right
      genes.append(_dna[gene_start:gene_stop])
      gene_start = None
      gene_stop = None
  
  # start reading frame 2
  for right in range(readframe_2 + 3, len(_dna) + 1, 3):
    codon = _dna[right-3:right]
    
    if codon == START_CODON:
      gene_start = right-3
    
    if codon in STOP_CODONS and gene_start is not None:
      gene_stop = right
      genes.append(_dna[gene_start:gene_stop])
      gene_start = None
      gene_stop = None
  
  # start reading frame 3
  for right in range(readframe_3 + 3, len(_dna) + 1, 3):
    codon = _dna[right-3:right]
    
    if codon == START_CODON:
      gene_start = right-3
    
    if codon in STOP_CODONS and gene_start is not None:
      gene_stop = right
      genes.append(_dna[gene_start:gene_stop])
      gene_start = None
      gene_stop = None

  return genes

# generates an mRNA strand from the given dna strand
def dna_to_rna(_dna):
  rna = ""
  for base in _dna:
    if base == "T":
      rna += "U"
    else:
      rna += base
  
  return rna

# generates a protein primary sequence (string of amino acids) from mRNA
def rna_to_amino_acid(mrna):
  protein = ""
  for i in range(0, len(mrna)-3, 3):
    codon = mrna[i:i + 3]
    AA = prot_dict[codon]
    protein += AA
  return(protein)

# completes entire process of converting dna to a protein primary sequence
def find_protein(dna):
  if len(dna) < 3:
    raise Exception("DNA sequence is not long enough")

  genes = find_genes(dna)

  # searches the complement
  complement = find_complement(dna)
  for gene in find_genes(complement):
    genes.append(gene)

  # save detected genes in a file
  with open("output/output_genes.txt", "w") as genes_output:
    output_str = ""
    for gene in genes:
      output_str += gene + ","
    genes_output.write(output_str)

  max = 0
  longest_gene = ""
  for gene in genes:
    if len(gene) > max:
      max = len(gene)
      longest_gene = gene
  
  mRNA = dna_to_rna(longest_gene)
  protein = rna_to_amino_acid(mRNA)

  return protein

def main():
  # takes command line arguments
  parser = argparse.ArgumentParser(description="Process bed files")

  # input file is only mandatory argument
  parser.add_argument("input_file")

  # specifies an output file path. by default this will output to "output/output.txt"
  parser.add_argument("-o", "--output")

  # Format: Position:Base
  parser.add_argument("-m", "--mutation", help="Finds effect of changing base X to Y at position Z")
  parser.add_argument("-p", "--probability", help="Probability that a mutation at position Z impact protein product")

  args = parser.parse_args()

  # checks file is a txt file
  file_suffix = os.path.splitext(args.input_file)[-1].lower()

  if file_suffix != '.txt':
    raise Exception("Error: please ensure only a .txt file is inputted")

  # check input file exists
  check_existance(args.input_file)

  # reads dna seq from input file
  input = open(args.input_file)

  if args.output is not None:
    check_existance("output/{}".format(args.output))
    output_file_name = args.output
  else:
    output_file_name = "output.txt"
  
  
  dna = input.read()

  # normalizes dna to uppercase
  dna = dna.upper()

  # check if dna contains only "A", "T", "C", "G"
  for base in dna:
    if base not in ["A", "T", "C", "G"]:
      raise Exception("Error: unknown base {} found".format(base))

  protein = find_protein(dna)

  # write out the output to a file
  with open("output/{}".format(output_file_name), "w") as output_file:
    output_file.write(protein)

  # mutation feature: determines if base Y at position Z changes the protein
  if args.mutation:
    fields = args.mutation.split(':')
    position = int(fields[0])

    if position > len(dna) or position < 0:
      raise Exception("Error: position {} out of bounds".format(position))

    new_base = fields[1]

    if new_base not in ["A", "T", "C", "G"]:
      raise Exception("Error: unknown base {} found".format(new_base))

    mutated_dna = dna[:position] + new_base + dna[position+1:]
    mutated_protein = find_protein(mutated_dna)

    # determines the effect by comparing protein outputs
    if mutated_protein == protein:
      print("A change at position {_position} to {_new_base} won't impact protein structure"
        .format(_position=position, _new_base=new_base))
    else:
      print("A change at position {_position} to {_new_base} will impact protein structure"
        .format(_position=position, _new_base=new_base))
  
  # finds the probability that a mutation at position X will change the protein
  if args.probability:
    position = int(args.probability)

    if position > len(dna) or position < 0:
      raise Exception("Error: position {} out of bounds".format(position))

    original_base = dna[position]

    # counts number of different proteins, and divides by total outcomes (3 possible bases)
    favorable_outcomes = 0
    total_outcomes = 3

    for sample in ['A', 'T', 'C', 'G']:
      if sample != original_base:
        mutated_dna = dna[:position] + sample + dna[position+1:]
        mutated_protein = find_protein(mutated_dna)

        if mutated_protein != protein:
          favorable_outcomes += 1
    
    print("The probability that a change at position", position, " impact protein structure is:", str(favorable_outcomes/total_outcomes))

  

if __name__ == "__main__":
  main()
