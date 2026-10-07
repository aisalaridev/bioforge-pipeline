class ORF:
    # Create an ORF object
    def __init__(self, Strand, Frame, Protein, start_pos, is_complete):
                  
        self.Strand = Strand              
        self.Frame = Frame
        self.Protein = Protein                
        self.start_pos = start_pos        
        self.is_complete = is_complete    

    def __str__(self):
        return f"Strand={self.Strand}, Frame={self.Frame}, Protein={self.Protein}, start_pos={self.start_pos}, is_complete={self.is_complete}"
        
def scan_frame(sequence, Frame, Strand):
    orfs = []
    orf_started = False
    orf_start = None

    # Check codons three nucleotides at a time
    for codon_start in range(Frame, len(sequence) - 2, 3):
        # Get three nucleotides
        codon = sequence[codon_start:codon_start + 3]

        if not orf_started:
            #start codon
            if codon == "AUG":
                orf_started = True
                orf_start = codon_start

        else:
            #stop codon
            if codon == "UAA" or codon == "UAG" or codon == "UGA":

                orf_sequence = sequence[orf_start:codon_start + 3]

                orfs.append(ORF(Strand=Strand, Frame=Frame, Protein=orf_sequence, start_pos=orf_start, is_complete=True))

                orf_started = False
                orf_start = None
    #orf naghes
    if orf_started:

        orf_sequence = sequence[orf_start:]

        orfs.append(ORF(Strand=Strand, Frame=Frame, Protein=orf_sequence, start_pos=orf_start, is_complete=False))

    return orfs
#three reading frame
def find_orfs_in_three_frames(sequence, Strand):
    all_orfs = []

    for Frame in range(3):

        orfs = scan_frame(sequence, Frame, Strand)

        all_orfs.extend(orfs)

    return all_orfs
#sixs reading frame
def find_orfs_in_six_frames(sequence):
    import sys
    if hasattr(sys.modules.get("__main__"), "DNA"):
        DNA = sys.modules["__main__"].DNA
    else:
        from main import DNA
    all_orfs = []

    dna = DNA(sequence)
    forward = find_orfs_in_three_frames(dna.dna2rna(), "Forward")

    all_orfs.extend(forward)
    #revers complemt in claas_DNA
    reverse_dna = DNA(dna.reverse_complement())
    reverse = find_orfs_in_three_frames(reverse_dna.dna2rna(), "Reverse")
    for orf in reverse:
        orf.start_pos = len(sequence) - 1 - orf.start_pos

    all_orfs.extend(reverse)

    return all_orfs
