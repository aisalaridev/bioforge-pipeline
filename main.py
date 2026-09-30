# Reading section : Needs to be implemented correctly

import re

header_pattern= r"^>\s*(?P<id>\S+)\s*(?P<desc>.*)$"

with open('fasta.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            result = re.findall(header_pattern, line)
            dna_id , description = result[0]
            #print(f"dna id: {dna_id} || description: {description}")
        elif re.fullmatch(r"[ACGT]+", line):
            dna_seq = re.fullmatch(r"[ACGT]+", line)
            dna_seq = dna_seq.group()
            #print(dna_seq)
#------------------------------------------------------------------------------------------------------------------------------------

# Working with dna sequence section

class DNA:

    def __init__(self, sequence):
        self.seq = sequence

    def check_dna_validity(self):
        for char in self.seq:
            if char in "ATCG":
                return True
            raise ValueError("Invalid DNA sequence")

    def complement(self):

        comp_seq = ''
        for char in self.seq:
            if char == 'A':
                comp_seq += 'T'
            elif char == 'T':
                comp_seq += 'A'
            elif char == 'C':
                comp_seq += 'G'
            elif char == 'G':
                comp_seq += 'C'
            else:
                raise Exception("Sequence is not valid! ")
        return comp_seq

    def reverse_complement(self):

        rev_seq = self.complement()
        rev_seq = rev_seq[::-1]
        return rev_seq

    def dna2rna(self):

        rna_seq = ""

        for char in self.seq:
            if char == 'T':
                rna_seq += 'U'
            else:
                rna_seq += char
        return rna_seq
    
    def gc_content(self):

        g_count = self.seq.count('G')
        c_count = self.seq.count('C')
        return round((g_count + c_count) / len(self.seq), 3)
    
dna = DNA(dna_seq)
#print(dna_seq)
#print(dna.complement())
#print(dna.reverse_complement())
#print(dna.dna2rna())
print(dna.gc_content())
