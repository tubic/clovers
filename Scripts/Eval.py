import sys
import os
import argparse
from typing import List

class Orf:
    """
    Represents an Open Reading Frame (ORF) with host, start, end, and strand information.
    """
    def __init__(self, host: str, start: int, end: int, strand: str):
        self.host = host
        sep = host.rfind('.')
        if sep >= 0: self.host = host[:sep]
        self.start = start
        self.end = end
        self.strand = strand

    def eqstrand(self, query: 'Orf') -> bool:
        """Check if two ORFs are on the same host and strand."""
        return self.host == query.host and self.strand == query.strand

    def match5end(self, query: 'Orf') -> bool:
        """
        Check if the 5' ends match.
        For positive strand: compare start positions;
        For negative strand: compare end positions.
        """
        if self.strand == '+':
            return self.start == query.start
        else:  # strand == '-'
            return self.end == query.end

    def match3end(self, query: 'Orf') -> bool:
        """
        Check if the 3' ends match.
        For positive strand: compare end positions;
        For negative strand: compare start positions.
        """
        if self.strand == '+':
            return self.end == query.end
        else:  # strand == '-'
            return self.start == query.start

    def __eq__(self, query: 'Orf') -> bool:
        """
        Two ORFs are equal if they share the same host/strand and have a matching 5' or 3' end.
        """
        return self.eqstrand(query) and (self.match5end(query) or self.match3end(query))


def parse_gff(file_path: str) -> List[Orf]:
    """
    Parse a GFF file and extract all CDS features as Orf objects.

    Parameters
    ----------
        file_path: str
            Path to the GFF file.

    Returns
    -------
        A list of Orf objects for each CDS entry.
    """
    orfs = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip empty lines and comment lines
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 8:
                    continue  # Malformed line, skip

                # Check that the feature type is "CDS" (column index 2 in 0-based)
                if fields[2] != 'CDS':
                    continue

                host = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]

                orfs.append(Orf(host, start, end, strand))

    except IOError:
        print(f"Error: cannot open file {file_path}.", file=sys.stderr)

    return orfs

def inter_count(pred_orfs: List[Orf], true_orfs: List[Orf]) -> int:
    """
    Compute the number of predicted ORFs that match true ORFs using the __eq__ method.
    Each true ORF can be matched at most once to avoid overcounting.
    """
    matched = [False] * len(true_orfs)
    intersection = 0
    for p in pred_orfs:
        for i, t in enumerate(true_orfs):
            if not matched[i] and p == t:
                matched[i] = True
                intersection += 1
                break
    return intersection

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate predicted ORFs against ground truth ORFs from GFF files."
    )
    parser.add_argument("-p", "--pred", required=True, help="Directory containing predicted GFF files.")
    parser.add_argument("-t", "--true", required=True, help="Directory containing ground truth GFF files.")
    args = parser.parse_args()

    pred_dir = args.pred
    true_dir = args.true

    # Check that both directories exist
    if not os.path.isdir(pred_dir):
        sys.exit(f"Error: prediction directory '{pred_dir}' does not exist.")
    if not os.path.isdir(true_dir):
        sys.exit(f"Error: truth directory '{true_dir}' does not exist.")

    # Collect all .gff files in prediction directory
    pred_files = [f for f in os.listdir(pred_dir) if f.endswith('.gff')]
    if not pred_files:
        sys.exit("Error: no .gff files found in prediction directory.")

    # Prepare data for table
    print(f"Assembly\tPred\tTrue\tTrue_Pos")
    total_files = len(pred_files)
    for idx, fname in enumerate(pred_files):
        print(f"Processing {fname} ({idx}/{total_files})", file=sys.stderr)
        acc = fname.replace(".gff", "")
        pred_path = os.path.join(pred_dir, fname)
        true_path = os.path.join(true_dir, fname)

        if not os.path.exists(true_path):
            print(f"Truth file '{true_path}' not found for prediction '{pred_path}'", file=sys.stderr)

        pred_orfs = parse_gff(pred_path)
        true_orfs = parse_gff(true_path)

        pred_count = len(pred_orfs)
        true_count = len(true_orfs)
        inter = inter_count(pred_orfs, true_orfs)

        print(f"{acc}\t{pred_count}\t{true_count}\t{inter}")    

if __name__ == "__main__":
    main()