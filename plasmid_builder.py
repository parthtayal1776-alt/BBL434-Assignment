import os

MARKER_DICT = {

    # Restriction enzymes
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SalI": "GTCGAC",
    "XbaI": "TCTAGA",
    "SmaI": "CCCGGG",
    "NotI": "GCGGCCGC",
    "SphI": "GCATGC",

    # Antibiotic resistance fragments 
    "Ampicillin": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTG",
    "KanR": "ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAG",
    "CmR": "ATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCCCAATGG",
    "TetR": "ATGTCTAGATTAGATAAAAGTAAAGTGATTTTAAAG",
    "SpecR": "ATGAAAAACATTGTTTATGTTTACGTTGAG",

    # Blue-white selection
    "lacZ_alpha": (
        "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTG"
        "ACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCAC"
        "ATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCG"
    ),

    "Blue_White_Selection": (
        "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTG"
        "ACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCAC"
        "ATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCG"
    )
}

# Read FASTA
def read_fasta(file_path):
    sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence.upper()

# GC Skew ORI Prediction
def find_ori_gc_skew(sequence, window_size=100):
    max_skew = float("-inf")
    ori_position = 0

    for i in range(len(sequence) - window_size):
        window = sequence[i:i + window_size]
        g = window.count("G")
        c = window.count("C")
        if (g + c) > 0:
            skew = (g - c) / (g + c)
            if skew > max_skew:
                max_skew = skew
                ori_position = i

    return ori_position

# Find dominant k-mer near ORI
def dominant_kmer(sequence, start, k=9, region=200):
    region_seq = sequence[start:start + region]
    kmer_count = {}

    for i in range(len(region_seq) - k):
        kmer = region_seq[i:i + k]
        kmer_count[kmer] = kmer_count.get(kmer, 0) + 1

    dominant = max(kmer_count, key=kmer_count.get)
    return dominant

# Parse Design File
def parse_design_file(design_file):

    marker_sequences = []

    with open(design_file, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split(",")

            if len(parts) < 2:
                continue

            marker_name = parts[1].strip()

            # Handle conceptual replication marker
            if marker_name == "High_Copy_Replication":
                print("High copy replication assumed via ORI.")
                continue

            if marker_name not in MARKER_DICT:
                raise ValueError(f"Marker '{marker_name}' not recognized")

            marker_sequences.append(MARKER_DICT[marker_name])

    return marker_sequences

# Write FASTA Output
def write_fasta(output_file, sequence):

    with open(output_file, "w") as f:
        f.write(">Synthetic_Plasmid\n")

        for i in range(0, len(sequence), 70):
            f.write(sequence[i:i+70] + "\n")

# Main Function
def main():

    base_dir = os.path.dirname(os.path.abspath(__file__))

    fasta_file = os.path.join(base_dir, "pUC19.fa")
    design_file = os.path.join(base_dir, "Design_pUC19.txt")
    output_file = os.path.join(base_dir, "Output.Fa")

    # Read backbone
    backbone = read_fasta(fasta_file)

    # Predict ORI
    ori_position = find_ori_gc_skew(backbone)
    print("Predicted ORI position:", ori_position)

    dominant = dominant_kmer(backbone, ori_position)

    # Extract ORI region (500 bp)
    ori_sequence = backbone[ori_position:ori_position + 500]

    # Parse markers
    markers = parse_design_file(design_file)

    # Build final plasmid
    final_sequence = ori_sequence

    for marker_seq in markers:
        final_sequence += marker_seq

    print("Final plasmid length:", len(final_sequence))

    # Write output
    write_fasta(output_file, final_sequence)

    print("Output written to Output.Fa")

if __name__ == "__main__":
    main()
