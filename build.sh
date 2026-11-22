#!/bin/bash
# Build recipe of this solver, as used for every build of the paper's runs.
# Run it from this directory:  bash build.sh <output directory>   ->  <output directory>/solver
set -e
set -o pipefail
CXX=${CXX:-g++}
command -v "$CXX" > /dev/null || { echo "no C++ compiler ($CXX) on PATH"; exit 1; }
for f in code/main.cpp code/methods/MCTS.cpp ext/cxxopts.hpp ext/fmt/core.h params.txt; do
  [ -f "$f" ] || { echo "$f is missing"; exit 1; }
done
python3 - <<'PY' || exit 1
import re, sys
seen = set()
for l in open("params.txt").read().splitlines():
    l = l.strip()
    if not l or l.startswith("#"):
        continue
    if not re.match(r"^[A-Za-z_]+\s+\S+$", l):
        sys.exit("params.txt: a line must read  option value  (found: %s)" % l)
    k = l.split()[0]
    if k in ('problem', 'instance', 'rand_seed', 'time_limit', 'target', 'output_directory', 'help'):
        sys.exit("params.txt: %s is set by the launcher and may not appear here" % k)
    if k in seen:
        sys.exit("params.txt: %s occurs twice" % k)
    seen.add(k)
PY
INC="-Iext -Icode -DFMT_HEADER_ONLY"
if [ "$2" = check ]; then                   # compile errors only, no binary
  for f in code/main.cpp code/*/*.cpp; do
    "$CXX" -O0 -DNDEBUG -std=c++17 -w -fsyntax-only $INC "$f" 2>&1 | tail -n 60
  done
  exit 0
fi
O="$1/obj"; rm -rf "$O"; mkdir -p "$O"
FLAGS="-O3 -flto=auto -DNDEBUG -std=c++17 -w $INC"
export CXX FLAGS O
ls code/main.cpp code/*/*.cpp | xargs -P 4 -I{} sh -c '"$CXX" $FLAGS -c "{}" -o "$O/$(echo "{}" | tr / _).o" 2>&1 | tail -n 60; test -s "$O/$(echo "{}" | tr / _).o"'
"$CXX" -O3 -flto=auto -s "$O"/*.o -o "$1/prog" 2>&1 | tail -n 60
rm -rf "$O"
cp params.txt "$1/"
cat > "$1/solver" <<'LAUNCHER'
#!/usr/bin/env python3
"""solver <instance.wcol> <solution> -t <seconds> -seed <int>: one run of the program (one thread). The program reads
an instance by NAME from ../instances/reduced_wvcp/<name>.col (edges) and <name>.col.w (weights, one per line)
relative to its working directory, and writes a CSV file whose last line holds its best colouring; this launcher
writes those two files into <solution>.work/, runs the program there, turns the last CSV line into the solution
file (one colour per vertex) and removes the directory. The target option is always 0 (never reached)."""
import os
import shutil
import signal
import subprocess
import sys
here = os.path.dirname(os.path.abspath(__file__))
a = sys.argv[1:]
inst, sol = os.path.abspath(a[0]), os.path.abspath(a[1])
T = float(a[a.index("-t") + 1])
seed = a[a.index("-seed") + 1]
W = sol + ".work"
shutil.rmtree(W, ignore_errors=True)
for d in ("build", "instances/reduced_wvcp", "out/tbt"):
    os.makedirs(os.path.join(W, d))
name = "g"
n = 0
wts, edges = {}, []
with open(inst) as f:
    for line in f:
        t = line.split()
        if not t:
            continue
        if t[0] == "p":
            n = int(t[2])
        elif t[0] == "n":
            wts[int(t[1])] = t[2]
        elif t[0] == "e":
            edges.append("e %s %s\n" % (t[1], t[2]))
with open(os.path.join(W, "instances", "reduced_wvcp", name + ".col"), "w") as f:
    f.write("p edge %d %d\n" % (n, len(edges)))
    f.writelines(edges)
with open(os.path.join(W, "instances", "reduced_wvcp", name + ".col.w"), "w") as f:
    f.write("".join("%s\n" % wts[i] for i in range(1, n + 1)))
cmd = [os.path.join(here, "prog"), "--problem", "wvcp", "--instance", name, "--rand_seed", seed,
       "--time_limit", "%d" % max(1, int(T - 3.0)), "--target", "0", "--output_directory", os.path.join(W, "out")]
for l in open(os.path.join(here, "params.txt")).read().splitlines():
    l = l.strip()
    if l and not l.startswith("#"):
        k, v = l.split()
        cmd += ["--" + k, v]
p = subprocess.Popen(cmd, cwd=os.path.join(W, "build"))


def stop(signum, frame):
    try:
        p.terminate()                 # the program then ends its current step and writes its last line
    except OSError:
        pass


signal.signal(signal.SIGTERM, stop)
signal.signal(signal.SIGINT, stop)
rc = p.wait()
last = None
for fn in ("%s_%s.csv" % (name, seed), "%s_%s.csv.running" % (name, seed)):
    q = os.path.join(W, "out", fn)
    if os.path.exists(q):
        with open(q) as f:
            for line in f:
                if line.startswith("#") or line.startswith("turn,") or not line.strip():
                    continue
                last = line.strip()
                t = last.split(",")
                if len(t) >= 4:
                    print("TRACE %s %s" % (t[1], t[-2]))
        break
if last:
    cols = last.split(",")[-1].split(":")
    if len(cols) == n:
        with open(sol, "w") as f:
            f.write("".join("%s\n" % c for c in cols))
shutil.rmtree(W, ignore_errors=True)
sys.exit(rc if rc >= 0 else 1)
LAUNCHER
chmod +x "$1/solver"
