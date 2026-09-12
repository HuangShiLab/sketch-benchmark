#!/usr/bin/env python
"""Append one timing row (shared schema, all tasks) from a `/usr/bin/time -v` file.
task,tool,tool_version,family,stage,dataset,params,threads,rep,wall_s,user_s,sys_s,max_rss_kb,input_bp,sketch_bytes,index_bytes,node,cpu_model,timestamp
input_bp comes from --bp-sidecar (reads.bp written by stage 02) when given, else from the
inputs (seqkit if available, else a Python count). Appends under an exclusive lock."""
import argparse, datetime, fcntl, os, re, socket, subprocess, sys, gzip
from pathlib import Path

FAMILY = {"mash": "F1", "bindash": "F1", "mash_screen": "F1", "sourmash": "F2", "sourmash_gather": "F2", "sourmash_compare": "F2",
          "sourmash_read": "F2", "sylph": "F2", "sylph_profile": "F2", "sylph_profile_bc": "F2", "sylph_read": "F2", "yacht": "F2",
          "cmash": "F2", "maxgeomhash": "F2", "i2brad": "F2", "syn2bani": "F2", "fast2brad_m": "F2", "hulk": "F3", "gsearch": "F3", "dashing2": "F4", "hypergen": "F7",
          "deacon": "F6", "deacon_panhuman": "F6", "deacon_syncmer": "F6", "kraken2": "F6",
          "bowtie2": "aln", "hostile": "aln", "minimap2": "aln", "minimap2_cov": "aln", "skani": "aln", "fastani": "aln", "simka": "count"}
HEADER = "task,tool,tool_version,family,stage,dataset,params,threads,rep,wall_s,user_s,sys_s,max_rss_kb,input_bp,sketch_bytes,index_bytes,node,cpu_model,timestamp\n"

def parse_time(path):
    d = {}
    for l in open(path):
        if "Elapsed (wall clock)" in l:
            t = l.split("):")[-1].strip().split(":")           # "h:mm:ss" or "m:ss" after the label's closing "):"
            t = [float(x) for x in t]; d["wall"] = sum(v * 60 ** i for i, v in enumerate(reversed(t)))
        elif "User time" in l: d["user"] = float(l.split(":")[-1])
        elif "System time" in l: d["sys"] = float(l.split(":")[-1])
        elif "Maximum resident set size" in l: d["rss"] = int(l.split(":")[-1])
    return d

def du(path):
    if not path or not os.path.exists(path): return 0
    p = Path(path)
    return p.stat().st_size if p.is_file() else sum(f.stat().st_size for f in p.rglob("*") if f.is_file())

def count_bp(files):
    total = 0
    for f in files:
        try:
            r = subprocess.run(["seqkit", "stats", "-T", f], capture_output=True, text=True, check=True)
            total += int(r.stdout.splitlines()[1].split("\t")[4]); continue
        except Exception:
            pass
        op = gzip.open if f.endswith(".gz") else open
        with op(f, "rt") as fh:
            for i, l in enumerate(fh):
                if i % 4 == 1: total += len(l.strip())
    return total

def tool_version(tool):
    exe = {"sourmash_gather": "sourmash", "sourmash_compare": "sourmash", "sourmash_read": "sourmash", "sylph_profile": "sylph",
           "sylph_profile_bc": "sylph", "sylph_read": "sylph", "mash_screen": "mash", "minimap2_cov": "minimap2",
           "deacon_panhuman": "deacon", "deacon_syncmer": "deacon-syncmer", "fastani": "fastANI"}.get(tool, tool)
    try:
        r = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=20)
        return (r.stdout or r.stderr).strip().splitlines()[0].replace(",", " ")
    except Exception:
        return "NA"

def cpu_model():
    try:
        for l in open("/proc/cpuinfo"):
            if l.startswith("model name"): return l.split(":", 1)[1].strip().replace(",", " ")
    except Exception: pass
    return "NA"

ap = argparse.ArgumentParser()
for k in ("task", "tool", "stage", "dataset", "params", "threads", "rep", "time-file", "out"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--input", nargs="*", default=[]); ap.add_argument("--bp-sidecar"); ap.add_argument("--index"); ap.add_argument("--sketch")
a = ap.parse_args()
t = parse_time(a.time_file)
bp = int(open(a.bp_sidecar).read().split()[0]) if a.bp_sidecar and os.path.exists(a.bp_sidecar) else (count_bp(a.input) if a.input else 0)
row = [a.task, a.tool, tool_version(a.tool), FAMILY.get(a.tool, "NA"), a.stage, a.dataset, a.params.replace(",", ";"), a.threads, a.rep,
       f"{t.get('wall', 0):.2f}", f"{t.get('user', 0):.2f}", f"{t.get('sys', 0):.2f}", t.get("rss", 0), bp, du(a.sketch), du(a.index),
       socket.gethostname(), cpu_model(), datetime.datetime.now().isoformat(timespec="seconds")]
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
with open(a.out, "a") as fh:
    fcntl.flock(fh, fcntl.LOCK_EX)
    if fh.tell() == 0: fh.write(HEADER)
    fh.write(",".join(map(str, row)) + "\n")
