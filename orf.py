class ProteinSequence:    
    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)

class ORF:
    def __init__(self, strand, frame, start_pos, protein, is_complete, orf_id):
        self.strand = strand
        self.frame = frame
        self.start_pos = start_pos
        self.protein = protein
        self.length = protein.length
        self.is_complete = is_complete
        self.orf_id = orf_id

    def format_report_entry(self):
        if self.is_complete:
            status = "complete"
        else:
            status = "partial"   

        return (
            f"{self.orf_id}, strand={self.strand}, frame={self.frame}, "
            f"start={self.start_pos}, length={self.length}, {status}"
        )  

stop_codons = ["UGA", "UAA", "UAG"]

bases = "UCAG"
amino = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codon_table = {}
k = 0
for a in bases:
    for b in bases:
        for c in bases:
            codon_table[a + b + c] = amino[k]
            k += 1

def translate(rna, start, end):
    protein = ""
    for i in range(start, end, 3):
        protein = protein + codon_table[rna[i: i + 3]]
    return protein


#Reverse complement

def reverse_complement(dna):
    paris = {"A": "T", "T": "A", "C": "G", "G": "C"}
    result = ""
    for base in dna:
        result = result + paris[base]
    return result[::-1]    # [::-1] reverses the string

# find all ORFs in one RNA string

def find_orfs(rna):
    result = []

    for frame in range(3):   #three frame; 0, 1, 2
        start_found = False
        start_pos = None

        for i in range(frame, len(rna) - 2, 3):
            codon = rna[i : i + 3]

            if start_found == False:
                if codon == "AUG":
                    start_found = True
                    start_pos = i

            else:
                if codon in stop_codons:
                    result.append((frame, start_pos, i, True))
                    start_found = False  #look for the next AUG

       #loop finished but start is still opene -> incomplete ORF
        if start_found:
            end = frame + ((len(rna) - frame) // 3) * 3
            result.append((frame, start_pos, end, False))

    return result

    # take DNA, return a list of ORF objects (both strands)

def find_all_orfs(dna):
        dna = dna.upper()
        total = len(dna)
        orfs = []
        counter = 1


        for strand in ["+", "-"]:
            if strand == "+":
                seq = dna
            else:
                seq = reverse_complement(dna) 

            rna = seq.replace("T", "U") 

            for frame, start, end, complete, in find_orfs(rna):
                protein = ProteinSequence(translate(rna, start, end))

                if strand == "+":
                    pos = start
                else:
                    pos = total -1 - start

                orfs.append(ORF(strand, frame, pos, protein, complete, f"ORF{counter}"))
                counter += 1

        return orfs

    # test (run only when this file is executed directly)  

if __name__ == "__main__":
        tests =[
            ("complete", "GAACTAATGGCTAATAG"),
            ("incomplete (no stop)", "GAACTAATGGCTAAA"),
            ("frame 2", "GGATGGCTTAAA"),
            ("minus strand", "CTATTTAGCCATCATAGC"),
            ("no AUG", "GAACTAAAAGGC"),
        ]
        for title, dna in tests:
            print("---", title, dna) 
            orfs = find_all_orfs(dna)
            if len(orfs) == 0:
                print("no ORF found")
            for orf in orfs:
                print(orf.format_report_entry())  

                                   
                 



