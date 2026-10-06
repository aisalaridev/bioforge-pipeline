from exceptions import InvalidSequenceError


class DNA:

    def __init__(self, sequence):
        self.__seq = sequence.upper()

    def __str__(self):
        return f"Current sequence: {self.__seq}"

    def get_sequence(self):
        return self.__seq

    def check_dna_validity(self):
        if self.__seq == "":
            raise InvalidSequenceError("Invalid DNA sequence")
        for char in self.__seq:
            if char not in "ATCG":
                raise InvalidSequenceError("Invalid DNA sequence")
        return "Valid"

    def complement(self):
        self.check_dna_validity()
        comp_seq = []
        for char in self.__seq:
            if char == "A":
                comp_seq.append("T")
            elif char == "T":
                comp_seq.append("A")
            elif char == "C":
                comp_seq.append("G")
            elif char == "G":
                comp_seq.append("C")
        return "".join(comp_seq)

    def reverse_complement(self):
        rev_seq = self.complement()
        rev_seq = rev_seq[::-1]
        return rev_seq

    def dna2rna(self):
        self.check_dna_validity()
        rna_seq = []
        for char in self.__seq:
            if char == "T":
                rna_seq.append("U")
            else:
                rna_seq.append(char)
        return "".join(rna_seq)

    def gc_content(self):
        self.check_dna_validity()
        g_count = self.__seq.count("G")
        c_count = self.__seq.count("C")
        return (g_count + c_count) / len(self.__seq) * 100
