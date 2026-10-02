#input section solved with help on <team-leader>

import re
import os
import logging
from exceptions import FastaFormatError,BioforgeError,DataFileError,InvalidsequenceError 

os.makedirs("output", exist_ok=True)

logging.basicConfig(filename="output/bioforge.log",level=logging.INFO)
header_pattern= r"^>\s*(?P<id>\S+)\s*(?P<desc>.*)$"

try:
    with open('fasta.txt', 'r', encoding="utf-8") as f:
         line = f.readline()
except:
    logging.error("FileNotFoundError")
    raise FileNotFoundError("File not found!")
else:
    if os.path.getsize('fasta.txt') == 0:
        print("File is empty.")
    
with open('fasta.txt', 'r', encoding="utf-8") as f:
        record={}
        id_saved=[]
        line = f.readline()
        line = line.strip()
        if (line != "") and (not line.startswith(">")):
            try:
                raise FastaFormatError("File starts with sequence or corrupt data.")
            except FastaFormatError as e:
                logging.error(e) 
        f.seek(0)
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                result = re.findall(header_pattern, line)
                dna_id , description = result[0]
                logging.info(f"processing {dna_id}")
                if dna_id in id_saved:
                    logging.warning(f"Duplicate id for {dna_id}")
                else:    
                    id_saved.append(dna_id)
            elif re.fullmatch(r"[ACGT]+", line.upper()):
                dna_seq = re.fullmatch(r"[ACGT]+", line.upper())
                dna_seq = dna_seq.group()
                record[dna_id]=dna_seq
            elif line == "":
                 pass
            else:
                try:
                    raise FastaFormatError("invalid sequence or invalid header")
                except FastaFormatError as e:
                    logging.error(e)
        for i in id_saved:
                if not i in record:
                    print(f"{i} doesn't have any sequences")  

# Working with dna sequence section

class DNA:

    def __init__(self, sequence):
        self.__seq = sequence

    def __str__(self):
        return f"Current sequence: {self.__seq}"

    def check_dna_validity(self):
        for char in self.__seq:
            if char in "ATCG":
                return True
            raise InvalidsequenceError("Invalid DNA sequence")

    def complement(self):

        comp_seq = ''
        for char in self.__seq:
            if char == 'A':
                comp_seq += 'T'
            elif char == 'T':
                comp_seq += 'A'
            elif char == 'C':
                comp_seq += 'G'
            elif char == 'G':
                comp_seq += 'C'
            else:
                raise InvalidsequenceError("Sequence is not valid! ")
        return comp_seq

    def reverse_complement(self):

        rev_seq = self.complement()
        rev_seq = rev_seq[::-1]
        return rev_seq

    def dna2rna(self):

        rna_seq = ""

        for char in self.__seq:
            if char == 'T':
                rna_seq += 'U'
            else:
                rna_seq += char
        return rna_seq
    
    def gc_content(self):

        g_count = self.__seq.count('G')
        c_count = self.__seq.count('C')
        return round((g_count + c_count) / len(self.__seq), 2)
