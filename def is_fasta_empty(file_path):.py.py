def is_fasta_empty(file_path):
    with open(file_path, "r") as f:
        content = f.read()

    return content.strip() == ""


# Example
if is_fasta_empty("fasta.txt"):
    print("FASTA file is empty")
else:
    print("FASTA file is not empty")