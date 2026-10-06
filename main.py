import argparse
import logging
import os
import re
import sys

from dna import DNA
from exceptions import BioForgeError, DataFileError, FastaFormatError, InvalidSequenceError
from filters import LengthFilter, MotifFilter, WeightFilter
from orf import ORFFinder, annotate


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CODON_FILE = os.path.join(BASE_DIR, "data", "codon_table.txt")
WEIGHT_FILE = os.path.join(BASE_DIR, "data", "amino_weights.txt")
HEADER_PATTERN = r"^>\s*(?P<id>\S+)\s*(?P<desc>.*)$"


def setup_logging(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    logger = logging.getLogger("bioforge")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()
    log_path = os.path.join(out_dir, "bioforge.log")
    handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    return logger


def load_codon_table(path):
    logger = logging.getLogger("bioforge")
    table = {}
    try:
        handle = open(path, "r", encoding="utf-8")
    except FileNotFoundError:
        logger.error("data file not found: %s", path)
        raise
    with handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if line == "" or line.startswith("#"):
                continue
            bits = line.split()
            codon_ok = len(bits) == 2 and re.fullmatch(r"[ACGU]{3}", bits[0])
            amino_ok = codon_ok and re.fullmatch(r"[A-Z*]", bits[1])
            if not amino_ok:
                message = f"DataFileError: {path} is invalid on line {line_number}"
                logger.error(message)
                raise DataFileError(message, line_number)
            codon, amino = bits
            if codon in table:
                message = f"DataFileError: {path} repeats {codon} on line {line_number}"
                logger.error(message)
                raise DataFileError(message, line_number)
            table[codon] = amino
    if not table:
        message = f"DataFileError: {path} has no codon rows"
        logger.error(message)
        raise DataFileError(message)
    return table


def load_amino_weights(path):
    logger = logging.getLogger("bioforge")
    weights = {}
    try:
        handle = open(path, "r", encoding="utf-8")
    except FileNotFoundError:
        logger.error("data file not found: %s", path)
        raise
    with handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if line == "" or line.startswith("#"):
                continue
            bits = line.split()
            if len(bits) != 2 or re.fullmatch(r"[A-Z]", bits[0]) is None:
                message = f"DataFileError: {path} is invalid on line {line_number}"
                logger.error(message)
                raise DataFileError(message, line_number)
            amino, raw_weight = bits
            try:
                weight = float(raw_weight)
            except ValueError:
                message = f"DataFileError: {path} is invalid on line {line_number}"
                logger.error(message)
                raise DataFileError(message, line_number)
            if weight <= 0 or amino in weights:
                message = f"DataFileError: {path} is invalid on line {line_number}"
                logger.error(message)
                raise DataFileError(message, line_number)
            weights[amino] = weight
    if not weights:
        message = f"DataFileError: {path} has no weight rows"
        logger.error(message)
        raise DataFileError(message)
    return weights


def read_fasta(path):
    logger = logging.getLogger("bioforge")
    try:
        handle = open(path, "r", encoding="utf-8")
    except FileNotFoundError:
        logger.error("input file not found: %s", path)
        raise
    except OSError:
        message = f"FastaFormatError: cannot read {path}"
        logger.error(message)
        raise FastaFormatError(message)
    with handle:
        raw = handle.read()
    if raw.strip() == "":
        message = "FastaFormatError: FASTA file is empty"
        logger.error(message)
        raise FastaFormatError(message)

    records = []
    seen = set()
    current = None
    saw_header = False
    leading_logged = False

    def flush():
        if current is None:
            return
        sequence = "".join(current["parts"])
        if sequence == "":
            logger.error("FastaFormatError: header %s has no sequence", current["id"])
            return
        current["sequence"] = sequence
        del current["parts"]
        records.append(current)

    for line_no, raw_line in enumerate(raw.splitlines(), start=1):
        line = raw_line.strip()
        if line == "":
            continue
        if line.startswith(">"):
            saw_header = True
            flush()
            match = re.fullmatch(HEADER_PATTERN, line)
            if not match:
                logger.error("FastaFormatError: invalid header on line %s", line_no)
                current = None
                continue
            dna_id = match.group("id")
            description = match.group("desc").strip()
            if dna_id in seen:
                logger.warning("duplicate id %s", dna_id)
            seen.add(dna_id)
            organism = None
            org_match = re.search(r"organism=(\S+)", description)
            if org_match:
                organism = org_match.group(1)
            if description:
                logger.info("processing %s %s", dna_id, description)
            else:
                logger.info("processing %s", dna_id)
            current = {
                "id": dna_id,
                "description": description,
                "organism": organism,
                "parts": [],
            }
            continue
        if current is None:
            if not saw_header:
                if not leading_logged:
                    logger.error("FastaFormatError: sequence appears before the first header")
                    leading_logged = True
            else:
                logger.error("FastaFormatError: sequence on line %s has no valid header", line_no)
            continue
        current["parts"].append(line)

    flush()
    return records


class Pipeline:
    def __init__(self, codon_table, filters, motif_filter=None):
        self.finder = ORFFinder(codon_table)
        self.filters = filters
        self.motif_filter = motif_filter

    def run(self, records):
        logger = logging.getLogger("bioforge")
        orfs = []
        for record in records:
            dna = DNA(record["sequence"])
            try:
                dna.check_dna_validity()
            except InvalidSequenceError:
                logger.error(
                    "InvalidSequenceError: record %s has an invalid sequence",
                    record["id"],
                )
                continue
            logger.info(
                "record %s length %d gc %.2f",
                record["id"],
                len(dna.get_sequence()),
                dna.gc_content(),
            )
            found = self.finder.find(dna, record["id"], record["organism"])
            orfs.extend(found)
        if self.motif_filter is not None:
            self.motif_filter.identify(orfs)
        for current_filter in self.filters:
            orfs = current_filter.apply(orfs)
        annotate(orfs)
        return orfs


def write_report(path, orfs):
    with open(path, "w", encoding="utf-8") as handle:
        if not orfs:
            handle.write("No ORF passed the filters.\n")
            return
        blocks = []
        for orf in orfs:
            lines = [
                f"ID: {orf.orf_id}",
                f"Record: {orf.source_id}",
            ]
            if orf.organism:
                lines.append(f"Organism: {orf.organism}")
            lines.append(f"Strand: {orf.strand}")
            lines.append(f"Frame: {orf.frame}")
            lines.append(f"Start Position: {orf.start_pos}")
            lines.append(f"Protein: {orf.protein.sequence}")
            lines.append(f"Status: {orf.status()}")
            if orf.motifs:
                shown = []
                for hit in orf.motifs:
                    shown.append(f"{hit.sequence}:{hit.position}")
                lines.append("Motifs: " + ", ".join(shown))
            blocks.append("\n".join(lines))
        handle.write("\n\n".join(blocks))
        handle.write("\n")


def parse_args(argv):
    parser = argparse.ArgumentParser(description="BioForge DNA analysis pipeline")
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-length", required=True, type=int, dest="min_length")
    parser.add_argument("--min-weight", type=float, dest="min_weight")
    parser.add_argument("--max-weight", type=float, dest="max_weight")
    parser.add_argument("--motif", action="append", dest="motifs")
    args = parser.parse_args(argv)
    if args.min_length < 0:
        parser.error("min-length cannot be negative")
    if (
        args.min_weight is not None
        and args.max_weight is not None
        and args.min_weight > args.max_weight
    ):
        parser.error("min-weight is greater than max-weight")
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        logger = setup_logging(args.out)
    except OSError:
        print(f"cannot use output directory: {args.out}", file=sys.stderr)
        return 1
    try:
        codon_table = load_codon_table(CODON_FILE)
        amino_weights = load_amino_weights(WEIGHT_FILE)
        records = read_fasta(args.input)
    except FileNotFoundError as err:
        name = err.filename if err.filename else "file"
        print(f"file not found: {name}", file=sys.stderr)
        return 1
    except BioForgeError as err:
        print(err, file=sys.stderr)
        return 1

    filters = [LengthFilter(args.min_length)]
    if args.min_weight is not None or args.max_weight is not None:
        filters.append(WeightFilter(amino_weights, args.min_weight, args.max_weight))
    motif_filter = None
    cleaned = []
    for motif in args.motifs or []:
        text = motif.strip().upper()
        if text and text not in cleaned:
            cleaned.append(text)
    if cleaned:
        motif_filter = MotifFilter(cleaned)
        filters.append(motif_filter)

    pipeline = Pipeline(codon_table, filters, motif_filter)
    try:
        orfs = pipeline.run(records)
    except BioForgeError as err:
        logger.error("%s", err)
        print(err, file=sys.stderr)
        return 1

    report_path = os.path.join(args.out, "report.txt")
    write_report(report_path, orfs)
    logger.info("wrote %s orf(s) to %s", len(orfs), report_path)
    print(f"wrote {len(orfs)} orf(s) to {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
