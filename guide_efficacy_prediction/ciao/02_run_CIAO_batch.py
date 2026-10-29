#!/usr/bin/env python3
"""Wrapper for CIAO predictions across multiple genes, with optional parallel jobs.

Inputs: a single-record genome FASTA, its matching NCBI GFF and optional locus-tag
list. Upstream: 01_fetch_CIAO_model.py and the environment in README.md.
Outputs: <output>/<locus_tag>/gRNAs.csv and run.log; one run_config.json records
inputs/settings so --resume cannot silently mix different analyses. No scoring
or guide-selection calculations are implemented here; the supplied predictor
performs them. Use -h for settings and examples in README.md.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import importlib.util
from importlib import metadata
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.parse import unquote

BIOTYPES = {"protein_coding", "rRNA", "tRNA", "ncRNA", "pseudogene"}
HERE = Path(__file__).resolve().parent

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def annotations(path):
    result = []
    with path.open() as f:
        for line in f:
            if line.startswith("#") or "Gene;gene" not in line:
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 9:
                continue
            attrs = dict(x.split("=", 1) for x in row[8].split(";") if "=" in x)
            identifier = attrs.get("ID", "")
            # Match the predictor's NCBI gene-ID parsing convention.
            locus = identifier.split("-")[1] if "-" in identifier else ""
            if not locus or "Name" not in attrs or "gene_biotype" not in attrs:
                raise ValueError("GFF lacks the NCBI gene attributes expected by the predictor.")
            result.append({"id": locus, "name": unquote(attrs["Name"]),
                           "biotype": attrs["gene_biotype"], "seqid": row[0]})
    if not result:
        raise ValueError("No compatible NCBI gene records found in GFF.")
    if len({r["id"] for r in result}) != len(result):
        raise ValueError("GFF contains duplicate gene IDs.")
    return result

def valid_result(path):
    if not path.is_file():
        return False
    try:
        with path.open(newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            required = {"gRNA_ID", "SequenceID", "gRNA sequence (20nt)", "Activity score"}
            if not required.issubset(reader.fieldnames or []):
                return False
            n = 0
            for row in reader:
                if not math.isfinite(float(row["Activity score"])):
                    return False
                if not re.fullmatch("[ACGT]{20}", row["gRNA sequence (20nt)"]):
                    return False
                n += 1
            return n > 0
    except (OSError, ValueError, KeyError, TypeError):
        return False

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--gff", type=Path, required=True)
    p.add_argument("--genes", type=Path, help="One locus tag or unique gene name per line; default: all compatible GFF genes.")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=1, help="Concurrent gene processes; default 1.")
    p.add_argument("--resume", action="store_true", help="Reuse complete outputs only when run inputs/settings match.")
    p.add_argument("--check-only", action="store_true", help="Validate reference/targets and report dependencies, without creating outputs.")
    p.add_argument("--biotype", choices=sorted(BIOTYPES) + ["all"], default="all")
    p.add_argument("--min-gc", type=float, default=30)
    p.add_argument("--max-gc", type=float, default=85)
    p.add_argument("--if-promoter", type=int, choices=[0, 1], default=0,
                   help="Uniform first-gene-in-operon flag; default 0, as in the predictor.")
    args = p.parse_args()
    if args.jobs < 1 or not 0 <= args.min_gc <= args.max_gc <= 100:
        p.error("Use jobs >=1 and 0 <= min-gc <= max-gc <=100.")
    fasta, gff, output = args.fasta.resolve(), args.gff.resolve(), args.output.resolve()
    if not fasta.is_file() or not gff.is_file():
        p.error("FASTA or GFF file is missing.")
    with fasta.open() as f:
        records = [x[1:].split()[0] for x in f if x.startswith(">")]
    if len(records) != 1:
        p.error("This genome/GFF workflow requires a single-record reference FASTA.")
    all_genes = annotations(gff)
    if set(r["seqid"] for r in all_genes) != set(records):
        p.error("GFF sequence identifiers do not match the reference FASTA.")
    by_id = {r["id"]:r for r in all_genes}
    targets = list(by_id)
    if args.genes:
        targets = []
        for value in args.genes.read_text().splitlines():
            value = value.strip()
            if not value or value.startswith("#"):
                continue
            matches = [by_id[value]] if value in by_id else [r for r in all_genes if r["name"] == value]
            if len(matches) != 1:
                p.error("Unknown or ambiguous gene: " + value + "; use its locus tag.")
            targets.append(matches[0]["id"])
        targets = list(dict.fromkeys(targets))
    allowed = BIOTYPES if args.biotype == "all" else {args.biotype}
    excluded = [x for x in targets if by_id[x]["biotype"] not in allowed]
    targets = [x for x in targets if by_id[x]["biotype"] in allowed]
    if excluded:
        print("Excluded by predictor biotype selection:", ", ".join(excluded), flush=True)
    if not targets:
        p.error("No genes remain after biotype selection.")
    if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", x) for x in targets):
        p.error("Gene IDs must be suitable for directory names.")
    print("Targets:", len(targets), "of", len(all_genes), "annotated genes", flush=True)
    predictor = HERE / "CRISPRi_design_MERF_progress.py"
    model_files = [HERE / "saved_model" / x for x in ["CRISPRi_headers.sav", "CRISPRi_model.sav"]]
    missing = [str(x) for x in model_files if not x.is_file()]
    for module in ["numpy", "pandas", "matplotlib", "Bio", "regex", "sklearn"]:
        if importlib.util.find_spec(module) is None:
            missing.append("Python module " + module)
    for tool in ["RNAduplex2.1.9h", "RNAduplex2.4.14", "RNAfold"]:
        if not shutil.which(tool):
            missing.append(tool)
    if missing:
        for x in missing:
            print("Missing dependency:", x, file=sys.stderr)
        if args.check_only:
            return 1
        p.error("Install missing dependencies; see README.md and 01_fetch_CIAO_model.py.")
    fold_version = subprocess.run(["RNAfold", "--version"], capture_output=True, text=True)
    fold_version_text = (fold_version.stdout + fold_version.stderr).strip()
    if fold_version.returncode or not re.search(r"(?<![0-9.])2\.1\.9h(?![A-Za-z0-9.])", fold_version_text):
        p.error("RNAfold must report version 2.1.9h for this adaptation; observed: " + fold_version_text)
    # Check serialization once before launching a genome-wide run.
    probe = subprocess.run([sys.executable, "-c",
        "import pickle,sys; "
        "h=pickle.load(open(sys.argv[1],'rb')); "
        "m=pickle.load(open(sys.argv[2],'rb')); "
        "assert len(h)==m.n_features_in_; print('Model loading passed:',len(h),'features')",
        str(model_files[0]), str(model_files[1])], capture_output=True, text=True)
    if probe.returncode:
        print(probe.stdout + probe.stderr, file=sys.stderr)
        p.error("Model loading failed in this Python environment; see environment.yml.")
    print(probe.stdout.strip(), flush=True)
    if args.check_only:
        print("References, targets, dependencies and model loading passed checks. RNA feature calculation was not run.")
        return 0
    versions = {}
    for distribution in ["numpy", "pandas", "scikit-learn", "biopython", "regex", "matplotlib", "tqdm", "shap"]:
        try:
            versions[distribution] = metadata.version(distribution)
        except metadata.PackageNotFoundError:
            versions[distribution] = "not installed"
    config = {"fasta_sha256":digest(fasta), "gff_sha256":digest(gff),
              "predictor_sha256":digest(predictor),
              "models":{x.name:digest(x) for x in model_files},
              "python":sys.version, "packages":versions, "rnafold_version":fold_version_text,
              "biotype":args.biotype, "min_gc":args.min_gc,
              "max_gc":args.max_gc, "if_promoter":args.if_promoter,
              "rna_tools":{x:{"path":shutil.which(x), "sha256":digest(Path(shutil.which(x)))}
                           for x in ["RNAduplex2.1.9h", "RNAduplex2.4.14", "RNAfold"]}}
    if output.exists():
        cfg = output / "run_config.json"
        if not args.resume or not cfg.is_file() or json.loads(cfg.read_text()) != config:
            p.error("Output folder exists with different/unknown settings, or --resume was not given; choose a new folder.")
    else:
        output.mkdir(parents=True)
        (output / "run_config.json").write_text(json.dumps(config, indent=2) + "\n")

    def run_gene(locus):
        dest = output / locus
        if dest.exists():
            if args.resume and valid_result(dest / "gRNAs.csv"):
                return locus, "reused"
            if args.resume and (dest / "no_guides.txt").is_file():
                return locus, "no_guides"
            return locus, "failed: existing incomplete folder; inspect it or choose a new output folder"
        work = output / ("." + locus + ".running")
        if work.exists():
            return locus, "failed: unfinished work folder " + str(work)
        command = [sys.executable, str(predictor), str(fasta), "-gff", str(gff),
                   "-targeting_genes", locus, "-if_promoter", str(args.if_promoter),
                   "-mingc", str(args.min_gc), "-maxgc", str(args.max_gc),
                   "-b", args.biotype, "-shap", "no", "-v", "-o", str(work)]
        log = output / ("." + locus + ".run.log")
        try:
            with log.open("w") as handle:
                result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
            if not work.is_dir():
                return locus, "failed: predictor did not create output; see " + str(log)
            log.replace(work / "run.log")
            if result.returncode == 0 and valid_result(work / "gRNAs.csv"):
                work.replace(dest)
                return locus, "completed"
            if result.returncode == 0 and "Error: No gRNA was found" in (work / "run.log").read_text():
                (work / "no_guides.txt").write_text("No guides passed the predictor's design filters.\n")
                work.replace(dest)
                return locus, "no_guides"
            return locus, "failed: see " + str(work / "run.log")
        except Exception as exc:
            return locus, "failed: " + str(exc)

    statuses = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_gene, locus) for locus in targets]
        for i, future in enumerate(as_completed(futures), 1):
            locus, status = future.result()
            statuses.append(status)
            print("[{}/{}] {}: {}".format(i, len(targets), locus, status), flush=True)
    failures = sum(x.startswith("failed:") for x in statuses)
    print("Finished:", len(statuses), "genes;", failures, "failures;",
          sum(x == "no_guides" for x in statuses), "without eligible guides.")
    return 1 if failures else 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as exc:
        print("ERROR:", exc, file=sys.stderr)
        raise SystemExit(2)
