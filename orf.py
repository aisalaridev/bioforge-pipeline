from dna import DNA
from exceptions import DataFileError


class ProteinSequence:
    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)

    def __str__(self):
        return self.sequence

    def molecular_weight(self, weights):
        total = 18.015
        for amino in self.sequence:
            if amino not in weights:
                raise DataFileError(f"no weight for amino acid {amino}")
            total += weights[amino]
        return round(total, 3)


class ORF:
    def __init__(self, protein, strand, frame, start_pos, is_complete, source_id, organism=None):
        self.protein = protein
        self.strand = strand
        self.frame = frame
        self.start_pos = start_pos
        self.is_complete = is_complete
        self.source_id = source_id
        self.organism = organism
        self.motifs = []
        self.orf_id = None

    def status(self):
        if self.is_complete:
            return "Complete"
        return "Incomplete"


class ORFFinder:
    def __init__(self, codon_table):
        self.codon_table = codon_table

    def find(self, dna, source_id, organism=None):
        found = []
        length = len(dna.get_sequence())
        forward = dna.dna2rna()
        reverse = DNA(dna.reverse_complement()).dna2rna()
        for frame in (0, 1, 2):
            found.extend(self._scan(forward, frame, "Forward", length, source_id, organism))
        for frame in (0, 1, 2):
            found.extend(self._scan(reverse, frame, "Reverse", length, source_id, organism))
        return found

    def _scan(self, rna, frame, strand, dna_length, source_id, organism):
        found = []
        index = frame
        while index + 3 <= len(rna):
            codon = rna[index:index + 3]
            if codon != "AUG":
                index += 3
                continue
            protein_text, complete, next_index = self._read(rna, index)
            if strand == "Reverse":
                start_pos = dna_length - 1 - index
            else:
                start_pos = index
            protein = ProteinSequence(protein_text)
            found.append(
                ORF(protein, strand, frame, start_pos, complete, source_id, organism)
            )
            if not complete:
                break
            index = next_index
        return found

    def _read(self, rna, start):
        protein = []
        pos = start
        while pos + 3 <= len(rna):
            codon = rna[pos:pos + 3]
            if codon in ("UAA", "UAG", "UGA"):
                return "".join(protein), True, pos + 3
            if codon not in self.codon_table:
                raise DataFileError(f"codon {codon} was not found in the codon table")
            amino = self.codon_table[codon]
            if amino == "*":
                return "".join(protein), True, pos + 3
            protein.append(amino)
            pos += 3
        return "".join(protein), False, len(rna)


def annotate(orfs):
    number = 1
    for orf in orfs:
        orf.orf_id = f"BFG_{number:03d}"
        number += 1
