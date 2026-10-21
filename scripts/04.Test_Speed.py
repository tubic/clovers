import argparse
import csv
import os
import subprocess
import tempfile
import time

from Bio import SeqIO

DEFAULT_TOOLS = ["clovers", "mga", "prodigal"]


def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Benchmark gene-finding tools (clovers/mga/prodigal) "
                    "on a set of genomes.")
    parser.add_argument("--target", required=True,
                        help="File with one genome accession per line.")
    parser.add_argument("--gbff-dir", action="append", required=True,
                        dest="gbff_dirs", metavar="DIR",
                        help="Directory containing <acc>.gbff files; may be "
                             "repeated to search multiple directories.")
    parser.add_argument("--bin", default="bin",
                        help="Directory containing the tool executables "
                             "(default: bin).")
    parser.add_argument("--output", default=None,
                        help="Output CSV path (default: <cwd>/speed_benchmark.csv).")
    parser.add_argument("--tools", default=",".join(DEFAULT_TOOLS),
                        help="Comma-separated tools to benchmark "
                             "(default: clovers,mga,prodigal).")
    parser.add_argument("--threads", type=int, default=16,
                        help="Threads passed to clovers via -T (default: 16).")
    return parser.parse_args(argv)


def find_gbff(acc, gbff_dirs):
    for d in gbff_dirs:
        p = os.path.join(d, acc + ".gbff")
        if os.path.isfile(p):
            return p
    return None


def run_timed(cmd, **kw):
    t0 = time.perf_counter()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   check=True, **kw)
    return time.perf_counter() - t0


def benchmark(acc, gbff_dirs, bin_dir, tools, threads):
    gbff = find_gbff(acc, gbff_dirs)
    if gbff is None:
        raise FileNotFoundError(acc)
    records = list(SeqIO.parse(gbff, "genbank"))
    size_bp = sum(len(r.seq) for r in records)
    row = {"Genome": acc, "Size_bp": size_bp}
    with tempfile.TemporaryDirectory() as wd:
        fasta = os.path.join(wd, "genome.fa")
        SeqIO.write(records, fasta, "fasta")
        if "clovers" in tools:
            row["clovers_sec"] = run_timed(
                [os.path.join(bin_dir, "clovers"), "-i", fasta,
                 "-o", os.devnull, "-T", str(threads), "-q"], cwd=wd)
        if "mga" in tools:
            row["mga_sec"] = run_timed(
                [os.path.join(bin_dir, "mga"), fasta, "-s"], cwd=wd)
        if "prodigal" in tools:
            row["prodigal_sec"] = run_timed(
                [os.path.join(bin_dir, "prodigal"), "-i", fasta, "-q",
                 "-o", os.devnull], cwd=wd)
    return row


def main(argv=None):
    args = parse_args(argv)
    tools = [t.strip() for t in args.tools.split(",") if t.strip()]
    out_csv = args.output or os.path.join(os.getcwd(), "speed_benchmark.csv")
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)

    with open(args.target) as fh:
        targets = [l.strip() for l in fh if l.strip()]

    done = set()
    if os.path.isfile(out_csv):
        with open(out_csv) as fh:
            done = {r["Genome"] for r in csv.DictReader(fh)}

    fields = ["Genome", "Size_bp"] + [t + "_sec" for t in tools]
    with open(out_csv, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        if not done:
            w.writeheader()
        for i, acc in enumerate(targets, 1):
            if acc in done:
                continue
            try:
                row = benchmark(acc, args.gbff_dirs, args.bin, tools,
                                args.threads)
            except Exception as e:  # noqa: BLE001 - keep the batch running
                row = {"Genome": acc, "Size_bp": "", "Status": repr(e)}
                print(f"[{i}/{len(targets)}] {acc} FAILED: {e}", flush=True)
                continue
            w.writerow(row)
            fh.flush()
            timings = " ".join(
                f"{t}={row.get(t + '_sec', '')}" for t in tools)
            print(f"[{i}/{len(targets)}] {acc} size={row['Size_bp']} "
                  f"{timings}", flush=True)

    print(f"Done. results -> {out_csv}", flush=True)


if __name__ == "__main__":
    main()
