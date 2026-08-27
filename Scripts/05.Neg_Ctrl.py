import argparse
import csv
import gzip
import hashlib
import logging
import os
import random
import re
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio import SeqIO
from Bio.Seq import Seq

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MIN_ORF_LEN = 90            # minimum ORF / decoy length (nt, incl. stop)
MIN_GAP_MARGIN = 20         # flanking margin kept around each decoy (nt)
START_CODONS = ("ATG", "GTG", "TTG")
START_CODON_WEIGHTS = (0.80, 0.15, 0.05)
STOP_CODONS = ("TAA", "TAG", "TGA")

STOPS_BY_TABLE = {
    1: ("TAA", "TAG", "TGA"),
    4: ("TAA", "TAG"),
    11: ("TAA", "TAG", "TGA"),
}

MANIFEST_FIELDS = ("accession", "group", "n_records", "record_lengths",
                   "merged_len", "transl_table", "n_orfs", "n_decoys_target",
                   "n_decoys_inserted")


# ---------------------------------------------------------------------------
# ORF model
# ---------------------------------------------------------------------------

def build_table(table_id):
    """Return (stop_set, start_set) for the requested genetic code."""
    stops = set(STOPS_BY_TABLE.get(table_id, STOPS_BY_TABLE[11]))
    return stops, set(START_CODONS)


def find_orfs(seq, min_len, stops, starts):
    """Find all ORFs >= ``min_len`` nt on both strands.

    An ORF is a stop-to-stop segment (same frame) that contains at least one
    start codon; its coordinates run from the first start codon after the
    previous in-frame stop to the closing stop codon (inclusive).  Partial
    ORFs truncated at sequence ends are ignored.

    Args:
        seq: Upper-case DNA string.
        min_len: Minimum ORF length in nt including the stop codon.
        stops: Set of stop codons for the genetic code in use.
        starts: Set of accepted start codons.

    Returns:
        Tuple ``(orfs, intervals)`` where ``orfs`` is a list of dicts with
        keys ``start``/``end`` (0-based half-open) and ``strand``, and
        ``intervals`` is the sorted list of covered (start, end) pairs.
    """
    stop_re = re.compile("|".join(sorted(stops)))
    orfs = []
    for strand, s in ((1, seq), (-1, str(Seq(seq).reverse_complement()))):
        n = len(s)
        for frame in range(3):
            seg_start = frame
            for m in stop_re.finditer(s):
                c = m.start()
                if (c - frame) % 3 or c + 3 > n:
                    continue
                # first start codon within this stop-to-stop segment
                for i in range(seg_start, c, 3):
                    if s[i:i + 3] in starts:
                        if c + 3 - i >= min_len:
                            if strand == 1:
                                orfs.append({"start": i, "end": c + 3,
                                             "strand": strand})
                            else:
                                # map reverse-complement coordinates back
                                # to the plus-strand coordinate system
                                orfs.append({"start": n - (c + 3),
                                             "end": n - i,
                                             "strand": strand})
                        break
                seg_start = c + 3
    orfs.sort(key=lambda o: (o["start"], o["end"]))
    intervals = sorted((o["start"], o["end"]) for o in orfs)
    return orfs, intervals


def covered_gaps(intervals, seq_len, margin):
    """Compute inter-ORF gaps usable for decoy placement.

    Args:
        intervals: Sorted (start, end) pairs covered by ORFs.
        seq_len: Length of the concatenated sequence.
        margin: Flanking margin (nt) that must remain free around a decoy.

    Returns:
        List of (start, end) free intervals, 0-based half-open, already
        shrunk by ``margin`` on each side (only intervals with
        end - start >= MIN_ORF_LEN are kept).
    """
    merged = []
    for s, e in intervals:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    gaps = []
    prev = 0
    for s, e in merged + [[seq_len, seq_len]]:
        lo, hi = prev + margin, s - margin
        if hi - lo >= MIN_ORF_LEN:
            gaps.append((lo, hi))
        prev = e
    return gaps


# ---------------------------------------------------------------------------
# Decoy synthesis
# ---------------------------------------------------------------------------

def codon_usage(orfs, seq, stops):
    """Compute sense-codon and stop-codon usage from the extracted ORFs.

    Codon frequencies are computed over the coding strand of every ORF
    (the genome's gene pool), excluding the terminal stop codon.

    Returns:
        (sense_counter, stop_counter): ``collections.Counter`` objects.
    """
    sense = Counter()
    stop = Counter()
    for o in orfs:
        s = seq[o["start"]:o["end"]]
        if o["strand"] == -1:
            s = str(Seq(s).reverse_complement())
        stop[s[-3:]] += 1
        for i in range(0, len(s) - 3, 3):
            codon = s[i:i + 3]
            if codon not in stops and "N" not in codon:
                sense[codon] += 1
    return sense, stop


def make_weighted_sampler(counter, fallback_keys, rng):
    """Return a zero-arg function sampling keys proportional to counts."""
    keys = [k for k in counter if counter[k] > 0] or list(fallback_keys)
    weights = [counter.get(k, 0) for k in keys]
    if not any(weights):
        weights = [1] * len(keys)
    cumulative = []
    total = 0
    for w in weights:
        total += w
        cumulative.append(total)

    def sample():
        x = rng.uniform(0, total)
        lo, hi = 0, len(cumulative) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if cumulative[mid] < x:
                lo = mid + 1
            else:
                hi = mid
        return keys[lo]

    return sample


def synthesize_decoy(length, sense_sampler, stop_sampler, rng):
    """Build one decoy ORF of ``length`` nt (coding strand).

    Body codons are drawn i.i.d. from the genome codon usage, the start
    codon from the fixed 0.80/0.15/0.05 mixture and the stop codon from the
    genome stop usage, so the translation contains no internal stop.
    """
    n_body = length // 3 - 2  # codons between start and stop
    start = rng.choices(START_CODONS, weights=START_CODON_WEIGHTS, k=1)[0]
    body = "".join(sense_sampler() for _ in range(n_body))
    return start + body + stop_sampler()


def place_decoys(gaps, decoy_lengths, rng):
    """Randomly place decoys inside inter-ORF gaps without overlap.

    Args:
        gaps: Free (start, end) intervals (margins already subtracted).
        decoy_lengths: List of decoy lengths (nt) to place.
        rng: ``random.Random`` instance.

    Returns:
        List of (start, end) placed decoy intervals (0-based half-open);
        shorter than ``decoy_lengths`` when gap space runs out.
    """
    free = [list(g) for g in gaps]
    placed = []
    for length in decoy_lengths:
        # choose uniformly among the gaps that can hold this decoy
        fitting = [i for i, (lo, hi) in enumerate(free)
                   if hi - lo >= length]
        if not fitting:
            continue  # decoy dropped (caller logs the reduction)
        idx = fitting[rng.randrange(len(fitting))]
        lo, hi = free[idx]
        s = rng.randint(lo, hi - length)
        placed.append((s, s + length))
        # split the gap around the placed decoy
        free.pop(idx)
        if s - lo >= MIN_ORF_LEN:
            free.append([lo, s])
        if hi - (s + length) >= MIN_ORF_LEN:
            free.append([s + length, hi])
    placed.sort()
    return placed


# ---------------------------------------------------------------------------
# Per-genome pipeline
# ---------------------------------------------------------------------------

def read_genbank(path):
    """Read every record of a (optionally gzipped) GenBank file.

    Returns:
        (records, transl_table): list of Bio.SeqRecord.SeqRecord and the
        modal ``/transl_table`` value across CDS features (default 11).
    """
    opener = gzip.open if path.endswith(".gz") else open
    tables = []
    with opener(path, "rt") as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    for rec in records:
        for feat in rec.features:
            if feat.type == "CDS":
                val = feat.qualifiers.get("transl_table")
                if val:
                    try:
                        tables.append(int(val[0]))
                    except ValueError:
                        pass
    transl_table = Counter(tables).most_common(1)[0][0] if tables else 11
    if transl_table not in STOPS_BY_TABLE:
        transl_table = 11  # fall back for exotic codes (e.g. 25)
    return records, transl_table


def genome_seed(global_seed, accession):
    """Derive a per-genome deterministic seed from seed + accession."""
    digest = hashlib.sha256(f"{global_seed}:{accession}".encode()).hexdigest()
    return int(digest[:16], 16)


def process_genome(args):
    """Worker: build the decoy genome for one accession.

    Args:
        args: Tuple ``(gbff_path, accession, group, out_dir, seed)``.

    Returns:
        Dict with manifest fields plus ``status``/``message``; never raises.
    """
    gbff_path, accession, group, out_dir, seed = args
    result = {"accession": accession, "group": group, "status": "ok",
              "message": ""}
    fa_path = os.path.join(out_dir, accession + ".fa")
    gff_path = os.path.join(out_dir, accession + ".decoy.gff")
    try:
        records, transl_table = read_genbank(gbff_path)
        if not records:
            raise ValueError("no sequence records parsed")
        seqs = [str(r.seq).upper() for r in records]
        merged = "".join(seqs)
        seq_len = len(merged)
        rng = random.Random(genome_seed(seed, accession))

        stops, starts = build_table(transl_table)
        orfs, intervals = find_orfs(merged, MIN_ORF_LEN, stops, starts)
        sense, stop = codon_usage(orfs, merged, stops)

        n_target = max(20, min(200, round(0.05 * len(orfs))))
        gaps = covered_gaps(intervals, seq_len, MIN_GAP_MARGIN)
        capacity = sum((hi - lo) // MIN_ORF_LEN for lo, hi in gaps)
        if capacity < n_target:
            result["message"] = (f"gap capacity {capacity} < target "
                                 f"{n_target}; N reduced")
            n_target = max(0, min(n_target, capacity))

        # sample decoy lengths from the empirical ORF length distribution;
        # resample (up to 5 draws) lengths that cannot fit any gap, then
        # place longest-first for better packing (lengths themselves are
        # still i.i.d. draws from the ORF length distribution)
        orf_lengths = [o["end"] - o["start"] for o in orfs]
        max_gap = max((hi - lo for lo, hi in gaps), default=0)
        decoy_lengths = []
        for _ in range(n_target):
            length = rng.choice(orf_lengths)
            for _try in range(4):
                if length <= max_gap:
                    break
                length = rng.choice(orf_lengths)
            decoy_lengths.append(length)
        decoy_lengths.sort(reverse=True)
        placements = place_decoys(gaps, decoy_lengths, rng)
        if len(placements) < n_target:
            result["message"] = (f"placed {len(placements)}/{n_target} "
                                 f"decoys; N reduced (insufficient gaps)")
        n_inserted = len(placements)

        sense_sampler = make_weighted_sampler(sense, ("AAA",), rng)
        stop_sampler = make_weighted_sampler(stop, STOP_CODONS, rng)

        # synthesise decoys (coding strand), then insert into the genome
        decoys = []
        for (s, e) in placements:
            decoys.append({"start": s, "end": e,
                           "strand": 1 if rng.random() < 0.5 else -1,
                           "nt": synthesize_decoy(e - s, sense_sampler,
                                                  stop_sampler, rng)})

        out = []
        prev = 0
        for d in decoys:
            out.append(merged[prev:d["start"]])
            insert = d["nt"] if d["strand"] == 1 else str(
                Seq(d["nt"]).reverse_complement())
            out.append(insert)
            prev = d["end"]
        out.append(merged[prev:])
        final_seq = "".join(out)

        # --- validation: every decoy translates without internal stops ---
        table_stops = stops
        for d in decoys:
            frag = final_seq[d["start"]:d["end"]]
            if d["strand"] == -1:
                frag = str(Seq(frag).reverse_complement())
            if frag[:3] not in starts or frag[-3:] not in table_stops:
                raise ValueError(f"decoy at {d['start']} invalid termini")
            body = set(frag[i:i + 3] for i in range(3, len(frag) - 3, 3))
            if body & table_stops:
                raise ValueError(f"decoy at {d['start']} has internal stop")

        with open(fa_path, "w") as fh:
            fh.write(f">{accession} decoy-augmented genome, "
                     f"{len(records)} record(s) concatenated\n")
            for i in range(0, len(final_seq), 70):
                fh.write(final_seq[i:i + 70] + "\n")

        with open(gff_path, "w") as fh:
            fh.write("##gff-version 3\n")
            fh.write(f"##sequence-region {accession} 1 {len(final_seq)}\n")
            for idx, d in enumerate(decoys, 1):
                fh.write("\t".join([
                    accession, "Neg_Ctrl", "CDS",
                    str(d["start"] + 1), str(d["end"]), ".",
                    "+" if d["strand"] == 1 else "-", "0",
                    f"ID=decoy_{idx};Note=synthetic decoy ORF"]) + "\n")

        result.update({
            "n_records": len(records),
            "record_lengths": ";".join(str(len(s)) for s in seqs),
            "merged_len": seq_len,
            "transl_table": transl_table,
            "n_orfs": len(orfs),
            "n_decoys_target": n_target,
            "n_decoys_inserted": n_inserted,
        })
    except Exception as exc:  # noqa: BLE001 - keep the batch running
        result["status"] = "failed"
        result["message"] = f"{type(exc).__name__}: {exc}"
    return result


# ---------------------------------------------------------------------------
# Batch driver
# ---------------------------------------------------------------------------

def list_group_genomes(gbff_root, group):
    """Return [(gbff_path, accession)] for one group."""
    group_dir = os.path.join(gbff_root, group)
    files = {}
    for name in sorted(os.listdir(group_dir)):
        acc = name
        for suffix in (".gbff.gz", ".gbff", ".fa.gz", ".fa"):
            if acc.endswith(suffix):
                acc = acc[:-len(suffix)]
                break
        else:
            continue
        files.setdefault(acc, os.path.join(group_dir, name))
    return [(p, a) for a, p in sorted(files.items())]


def load_done_accessions(manifest_path):
    """Read accessions already present in a manifest (for resume)."""
    done = set()
    if os.path.exists(manifest_path):
        with open(manifest_path, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("accession"):
                    done.add(row["accession"])
    return done


def process_group(group, genomes, out_root, seed, workers, logger):
    """Run the decoy pipeline for one group and write its manifest."""
    out_dir = os.path.join(out_root, group)
    os.makedirs(out_dir, exist_ok=True)
    manifest_path = os.path.join(out_root, f"manifest_{group}.csv")
    done = load_done_accessions(manifest_path)

    pending = [(p, a) for p, a in genomes
               if a not in done
               or not (os.path.exists(os.path.join(out_dir, a + ".fa"))
                       and os.path.exists(
                           os.path.join(out_dir, a + ".decoy.gff")))]
    logger.info("%s: %d genomes total, %d already done, %d pending",
                group, len(genomes), len(genomes) - len(pending),
                len(pending))
    if not pending:
        return

    tasks = [(p, a, group, out_dir, seed) for p, a in pending]
    rows, failed = [], 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(process_genome, t): t[1] for t in tasks}
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            if res["status"] != "ok":
                failed += 1
                logger.error("%s/%s failed: %s",
                             group, res["accession"], res["message"])
            else:
                if res["message"]:
                    logger.info("%s/%s: %s",
                                group, res["accession"], res["message"])
                rows.append(res)
            if i % 50 == 0 or i == len(tasks):
                logger.info("  %s progress: %d/%d (failed so far: %d)",
                            group, i, len(tasks), failed)
                flush_manifest(manifest_path, rows, append=True)
                rows.clear()
    flush_manifest(manifest_path, rows, append=True)
    logger.info("%s done: %d ok, %d failed", group,
                len(tasks) - failed, failed)


def flush_manifest(manifest_path, rows, append):
    """Append manifest rows (creating the header on first write)."""
    if not rows:
        return
    write_header = not (append and os.path.exists(manifest_path))
    with open(manifest_path, "a" if append else "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS,
                                extrasaction="ignore")
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def setup_logger(output_root):
    """Configure console + file logging (``Neg_Ctrl.log``)."""
    os.makedirs(output_root, exist_ok=True)
    logger = logging.getLogger("Neg_Ctrl")
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    fh = logging.FileHandler(os.path.join(
        output_root, "Neg_Ctrl.log"), mode="a", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(ch)
    logger.addHandler(fh)
    return logger


def parse_args(argv=None):
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description="Generate decoy spiked-in negative-control genomes.")
    p.add_argument("--gbff-root", help="Root with dirs containing GenBank files")
    p.add_argument("--output-root", help="Output root for decoy genomes ")
    p.add_argument("--groups", nargs="+", help="Groups to process.")
    p.add_argument("--seed", type=int, default=42, help="Global random seed (default: 42).")
    p.add_argument("--workers", type=int, default=8, help="Parallel worker processes (default: 8).")
    p.add_argument("--limit", type=int, default=None, help="Process at most N pending genomes per group.")
    return p.parse_args(argv)


def main(argv=None):
    """Entry point."""
    args = parse_args(argv)
    logger = setup_logger(args.output_root)
    logger.info("=== decoy genome generation === seed=%d workers=%d",
                args.seed, args.workers)
    for group in args.groups:
        genomes = list_group_genomes(args.gbff_root, group)
        if args.limit:
            genomes = genomes[:args.limit]
        process_group(group, genomes, args.output_root, args.seed,
                      args.workers, logger)
    logger.info("All done.")


if __name__ == "__main__":
    main()
