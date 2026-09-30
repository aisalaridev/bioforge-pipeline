import re
import os

header_pattern= r"^>\s*(?P<id>\S+)\s*(?P<desc>.*)$"

try:
    with open('fasta.txt', 'r') as f:
         line = f.readline()
except:
    raise FileExistsError("File not found!")
else:
    if os.path.getsize('fasta.txt') == 0:
        print("File is empty.")
    
with open('fasta.txt', 'r') as f:
        line = f.readline()
        line = line.strip()
        if (line != "") and (not line.startswith(">")):
                raise ValueError("File starts with sequence or corrupt data.")
        f.seek(0)
        data_list = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                result = re.findall(header_pattern, line)
                dna_id , description = result[0]
            elif re.fullmatch(r"[ACGT]+", line):
                dna_seq = re.fullmatch(r"[ACGT]+", line)
                dna_seq = dna_seq.group()
            elif line == "":
                 pass
            else:
                 raise ValueError("Blah Blah")
            