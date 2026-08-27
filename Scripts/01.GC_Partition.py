import argparse
import logging
import os
import sys
from collections import defaultdict

from Bio import SeqIO

# FASTA extensions accepted as input.
FASTA_EXTENSIONS = (".fna", ".fasta", ".fa")

# Sequence classes processed by default (subdirectory names).
SEQUENCE_CLASSES = ("Positive", "Negative")


def setup_logger(output_dir):
    """Configure a logger writing to both the console and a log file.

    Args:
        output_dir: Directory in which ``GC_Partition.log`` is created.

    Returns:
        logging.Logger: The configured logger instance.
    """
    logger = logging.getLogger("GC_Partition")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        file_handler = logging.FileHandler(
            os.path.join(output_dir, "GC_Partition.log"), mode="w",
            encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def calculate_gc_content(sequence):
    """Compute the GC content of a nucleotide sequence.

    Args:
        sequence: A Bio.Seq.Seq object or a plain string.

    Returns:
        float: GC content as a percentage in the range [0, 100].
    """
    sequence = str(sequence).upper()
    total_length = len(sequence)
    if total_length == 0:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / total_length * 100.0


def get_gc_bin(gc_content, min_gc, max_gc, step):
    """Map a GC-content value to its bin label.

    Args:
        gc_content: GC content percentage.
        min_gc: Lower bound of the accepted GC range (inclusive).
        max_gc: Upper bound of the accepted GC range (exclusive).
        step: Bin width in percentage points.

    Returns:
        str or None: Bin label such as "25_26", or None if the value lies
        outside the accepted range.
    """
    if gc_content < min_gc or gc_content >= max_gc:
        return None
    lower = int(gc_content) - (int(gc_content) % step)
    lower = max(lower, min_gc)
    upper = lower + step
    return f"{lower:02d}_{upper:02d}"


def list_fasta_files(input_dir):
    """Return the sorted names of all FASTA files in a directory."""
    return sorted(f for f in os.listdir(input_dir)
                  if f.lower().endswith(FASTA_EXTENSIONS))


def process_class(input_dir, output_dir, seq_class, min_gc, max_gc, step,
                  logger):
    """Partition all sequences of one class into GC-bin FASTA files.

    Args:
        input_dir: Root input directory containing the class subdirectory.
        output_dir: Root output directory.
        seq_class: Class name (subdirectory name, e.g. "Positive").
        min_gc: Lower GC bound (inclusive).
        max_gc: Upper GC bound (exclusive).
        step: Bin width in percentage points.
        logger: Logger instance.

    Returns:
        dict or None: Mapping of bin label to sequence count, or None if the
        class could not be processed.
    """
    class_input_dir = os.path.join(input_dir, seq_class)
    if not os.path.isdir(class_input_dir):
        logger.error("Input directory not found: %s", class_input_dir)
        return None

    class_output_dir = os.path.join(output_dir, seq_class)
    os.makedirs(class_output_dir, exist_ok=True)

    fasta_files = list_fasta_files(class_input_dir)
    if not fasta_files:
        logger.warning("No FASTA files found in %s", class_input_dir)
        return None

    logger.info("Processing %s sequences: %d FASTA files",
                seq_class, len(fasta_files))

    bin_counts = defaultdict(int)
    total_sequences = 0
    written_sequences = 0
    file_handles = {}

    try:
        for file_idx, file_name in enumerate(fasta_files, 1):
            file_path = os.path.join(class_input_dir, file_name)
            try:
                for record in SeqIO.parse(file_path, "fasta"):
                    total_sequences += 1
                    gc_content = calculate_gc_content(record.seq)
                    gc_bin = get_gc_bin(gc_content, min_gc, max_gc, step)
                    if gc_bin is None:
                        continue
                    # Lazily open one output handle per GC bin.
                    if gc_bin not in file_handles:
                        bin_path = os.path.join(
                            class_output_dir,
                            f"{seq_class}_{gc_bin}.fasta")
                        file_handles[gc_bin] = open(bin_path, "w")
                    SeqIO.write(record, file_handles[gc_bin], "fasta")
                    bin_counts[gc_bin] += 1
                    written_sequences += 1
            except Exception as exc:  # noqa: BLE001 - keep batch running
                logger.error("Failed to process file %s: %s", file_name, exc)
                continue

            if file_idx % 100 == 0 or file_idx == len(fasta_files):
                logger.info("  %s progress: %d/%d files",
                            seq_class, file_idx, len(fasta_files))
    finally:
        for handle in file_handles.values():
            handle.close()

    logger.info("%s done: %d sequences read, %d written into %d GC bins",
                seq_class, total_sequences, written_sequences,
                len(bin_counts))
    return dict(bin_counts)


def generate_distribution_plot(bin_counts, output_dir, seq_class, logger):
    """Save a bar chart of the GC-bin distribution for one class.

    Args:
        bin_counts: Mapping of bin label to sequence count.
        output_dir: Directory in which the figure is saved.
        seq_class: Class name used in the title and file name.
        logger: Logger instance.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib is not installed; skipping the plot "
                       "for %s", seq_class)
        return

    bins = sorted(bin_counts)
    counts = [bin_counts[b] for b in bins]
    midpoints = [(int(b.split("_")[0]) + int(b.split("_")[1])) / 2
                 for b in bins]

    fig, ax = plt.subplots(figsize=(12, 5))
    color = "#4292c6" if seq_class == "Positive" else "#fb6a4a"
    ax.bar(midpoints, counts, width=0.8, color=color)
    ax.set_title(f"{seq_class} sequences: GC-content distribution")
    ax.set_xlabel("GC content (%)")
    ax.set_ylabel("Number of sequences")
    ax.set_xticks(midpoints[::5])
    ax.set_xticklabels(bins[::5], rotation=45)
    ax.grid(axis="y", ls=":", alpha=0.5)
    fig.tight_layout()

    figure_path = os.path.join(output_dir,
                               f"{seq_class}_GC_distribution.png")
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)
    logger.info("Distribution plot saved: %s", figure_path)


def write_report(results, output_dir, min_gc, max_gc, step, logger):
    """Write the plain-text summary report.

    Args:
        results: Mapping of class name to its bin-count dictionary.
        output_dir: Directory in which the report is saved.
        min_gc: Lower GC bound used for partitioning.
        max_gc: Upper GC bound used for partitioning.
        step: Bin width used for partitioning.
        logger: Logger instance.
    """
    report_path = os.path.join(output_dir, "GC_partition_report.txt")
    with open(report_path, "w", encoding="utf-8") as report:
        report.write("GC-content partition summary report\n")
        report.write("=" * 50 + "\n\n")
        report.write(f"GC range: [{min_gc}%, {max_gc}%)\n")
        report.write(f"Bin width: {step}%\n\n")

        for seq_class, bin_counts in results.items():
            if not bin_counts:
                continue
            total = sum(bin_counts.values())
            report.write(f"{seq_class} sequences\n")
            report.write("-" * 30 + "\n")
            report.write(f"Total: {total}\n\n")
            for gc_bin in sorted(bin_counts):
                count = bin_counts[gc_bin]
                percentage = count / total * 100.0 if total else 0.0
                report.write(f"  {gc_bin}: {count} ({percentage:.2f}%)\n")
            report.write("\n")

    logger.info("Summary report saved: %s", report_path)


def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Partition Positive/Negative sequence datasets into "
                    "GC-content bins.")
    parser.add_argument(
        "--input", required=True,
        help="Input root directory containing the Positive and Negative "
             "subdirectories.")
    parser.add_argument(
        "--output", required=True,
        help="Output root directory for the GC-binned FASTA files.")
    parser.add_argument(
        "--min-gc", type=float, default=20,
        help="Lower GC-content bound in percent, inclusive (default: 20).")
    parser.add_argument(
        "--max-gc", type=float, default=80,
        help="Upper GC-content bound in percent, exclusive (default: 80).")
    parser.add_argument(
        "--step", type=int, default=1,
        help="GC bin width in percentage points (default: 1).")
    parser.add_argument(
        "--no-plot", action="store_true",
        help="Do not generate GC-distribution plots.")
    return parser.parse_args(argv)


def main(argv=None):
    """Entry point of the GC partitioning pipeline."""
    args = parse_args(argv)
    logger = setup_logger(args.output)

    logger.info("=== GC-content dataset partition ===")
    logger.info("Input directory : %s", args.input)
    logger.info("Output directory: %s", args.output)
    logger.info("GC range: [%.0f%%, %.0f%%), bin width: %d%%",
                args.min_gc, args.max_gc, args.step)

    results = {}
    for seq_class in SEQUENCE_CLASSES:
        bin_counts = process_class(args.input, args.output, seq_class,
                                   args.min_gc, args.max_gc, args.step,
                                   logger)
        results[seq_class] = bin_counts
        if bin_counts and not args.no_plot:
            generate_distribution_plot(bin_counts, args.output, seq_class,
                                       logger)

    write_report(results, args.output, args.min_gc, args.max_gc, args.step,
                 logger)
    logger.info("All done.")


if __name__ == "__main__":
    main()
