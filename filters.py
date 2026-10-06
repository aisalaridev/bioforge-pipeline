class Motif:
    def __init__(self, sequence, position):
        self.sequence = sequence
        self.position = position


class Filter:
    def apply(self, orfs):
        raise NotImplementedError


class LengthFilter(Filter):
    def __init__(self, min_length):
        self.min_length = min_length

    def apply(self, orfs):
        kept = []
        for orf in orfs:
            if orf.protein.length >= self.min_length:
                kept.append(orf)
        return kept


class WeightFilter(Filter):
    def __init__(self, amino_weights, min_weight=None, max_weight=None):
        self.amino_weights = amino_weights
        self.min_weight = min_weight
        self.max_weight = max_weight

    def apply(self, orfs):
        kept = []
        for orf in orfs:
            weight = orf.protein.molecular_weight(self.amino_weights)
            if self.min_weight is not None and weight < self.min_weight:
                continue
            if self.max_weight is not None and weight > self.max_weight:
                continue
            kept.append(orf)
        return kept


class MotifFilter(Filter):
    def __init__(self, motif):
        if isinstance(motif, str):
            motif = [motif]
        self.motifs = []
        for item in motif:
            text = item.strip().upper()
            if text and text not in self.motifs:
                self.motifs.append(text)

    def identify(self, orfs):
        for orf in orfs:
            hits = []
            sequence = orf.protein.sequence
            for motif in self.motifs:
                start = 0
                while True:
                    found = sequence.find(motif, start)
                    if found < 0:
                        break
                    hits.append(Motif(motif, found))
                    start = found + 1
            hits.sort(key=lambda item: (item.position, item.sequence))
            orf.motifs = hits

    def apply(self, orfs):
        kept = []
        for orf in orfs:
            if orf.motifs:
                kept.append(orf)
        return kept
