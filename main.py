import re
import os
import logging
from exceptions import FastaFormatError,BioforgeError,DataFileError,InvalidsequenceError 

os.makedirs("output", exist_ok=True)
logging.basicConfig(filename="output/bioforge.log",level=logging.INFO)

header_pattern= r"^>\s*(?P<id>\S+)\s*(?P<desc>.*)$"
sequence_pattern = r"[ACGT]+"

#Check if file has non-DNA related content

def fasta_is_empty(loc):
    f = None
    try:
        f = open(loc, "r", encoding="utf-8")
    except:
        logging.error("FileNotFoundError")
        raise FileNotFoundError("File not found!")
    else:
        for line in f:
            line = line.strip()
            if ((line.startswith(">") and re.fullmatch(header_pattern, line))
                or re.fullmatch(sequence_pattern, line.upper())):
                return "non-empty file"
        raise FileNotFoundError("No DNA contents in file.")
    finally:
        if f:
            f.close()

#Check if file exists, also checking size of file in section else
#                
f = None
try:
    f = open('fasta.txt', 'r', encoding="utf-8")
except:
    logging.error("FileNotFoundError")
    raise FileNotFoundError("File not found!")
else:
    if os.path.getsize('fasta.txt') == 0:
        raise FileNotFoundError("File is empty!")
finally:
    if f:
        f.close()

#Check fasta file status

fasta_is_empty('fasta.txt')

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
        dna_id = None
        fresh = True
        skip_seq = False
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                result = re.findall(header_pattern, line)
                if not result:
                    logging.error("invalid header")
                    dna_id = None
                    skip_seq = True
                    continue
                dna_id , description = result[0]
                logging.info(f"processing {dna_id}")
                if dna_id in id_saved:
                    logging.warning(f"Duplicate id for {dna_id}")
                else:    
                    id_saved.append(dna_id)
                fresh = True
                skip_seq = False
            elif re.fullmatch(r"[ACGT]+", line.upper()):
                if skip_seq or dna_id is None:
                    continue
                dna_seq = re.fullmatch(r"[ACGT]+", line.upper())
                dna_seq = dna_seq.group()
                if fresh or dna_id not in record:
                    record[dna_id] = dna_seq
                    fresh = False
                else:
                    record[dna_id] += dna_seq
            elif line == "":
                 pass
            else:
                try:
                    raise InvalidsequenceError("invalid sequence or invalid header")
                except InvalidsequenceError as e:
                    logging.error(e)
                skip_seq = True
                if dna_id in record:
                    del record[dna_id]
                fresh = True
        for i in id_saved:
                if not i in record:
                    logging.error(f"{i} doesn't have any sequences")  

# Working with dna sequence section

class DNA:

    def __init__(self, sequence):
        self.__seq = sequence.upper()

    def __str__(self):
        return f"Current sequence: {self.__seq}"

    def check_dna_validity(self):
        if self.__seq == "":
            logging.info(f"Invalid DNA sequence.")
            raise InvalidsequenceError("Invalid DNA sequence")
        for char in self.__seq:
            if not char in "ATCG":
                logging.info(f"Invalid DNA sequence.")
                raise InvalidsequenceError("Invalid DNA sequence")
            else: 
                pass
        return f"Valid"

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
                logging.info(f"Unvalid sequence.")
                raise InvalidsequenceError("Sequence is not valid! ")
        return comp_seq

    def reverse_complement(self):

        rev_seq = self.complement()
        rev_seq = rev_seq[::-1]
        return rev_seq

    def dna2rna(self):
        self.check_dna_validity()

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
        return round((g_count + c_count) / len(self.__seq) * 100, 2)
#checking data files

def check_codon_table(path):
    with open(path,"r",encoding="utf-8") as f:
           n=0
           for line in f:
                    n+=1
                    line=line.strip()
                    if line.startswith("#") or not line:
                        pass
                    elif re.fullmatch(r"^[ACGU]{3}\s+[A-Z]$",line):
                        pass
                    elif re.fullmatch(r"^[ACGU]{3}\s+[*]$",line): 
                        pass
                    else:
                        logging.error(f"{path}is invalid:DataFileError  because of {line} in line {n} ")
                        raise DataFileError(f"{path} is invalid because of {line} in line {n}")
           return f"{path} is valid"
print(check_codon_table("data/codon_table.txt"))
def check_amino_weights(path):
       with open(path,"r",encoding="utf-8") as f:
            n=0
            for line in f:
                n+=1
                line=line.strip()
                if line.startswith("#") or not line:
                    pass
                elif re.fullmatch(r"^[A-Z]{1}\s+\d+[.]\d+$",line):
                    pass
                else:
                    logging.error(f"{path}is invalid:DataFileError  because of {line} in line {n} ") 
                    raise DataFileError(f"{path} is invalid because of {line} in line {n}")
            return f"{path} is valid"
print(check_amino_weights("data/amino_weights.txt"))

from ORF import find_orfs_in_six_frames
from orf import find_all_orfs

for dna_id in id_saved:
    if dna_id not in record:
        continue
    try:
        dna = DNA(record[dna_id])
        print(dna_id, dna.check_dna_validity(), dna.gc_content())
    except InvalidsequenceError as e:
        logging.error(e)
        continue
    for found in find_orfs_in_six_frames(record[dna_id]):
        print(dna_id, found)
    for found in find_all_orfs(record[dna_id]):
        print(found.format_report_entry())