#!/bin/bash
# Run the small Fe16+ test and compare the results against the reference.
# Usage: bash run.sh [nprocs]
# Set BIN to use executables other than the ones in <repo>/bin.
set -euo pipefail

nprocs=${1:-1}
here=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:-$here/../../bin}
BIN=$(cd "$BIN" && pwd)

# A hung program should fail the test, not stall it.
limit=300

# Stop on any program that fails, and show the end of its output.
step() {
    local name=$1 log=$2
    shift 2
    if ! timeout "$limit" "$@" < "${input:-/dev/null}" > "$log" 2>&1; then
        echo "FAIL: $name failed or timed out, last lines of $log:"
        tail -20 "$log"
        exit 1
    fi
}

# Work in a temp directory, which is kept for inspection.
work=$(mktemp -d "${TMPDIR:-/tmp}/pci_small.XXXXXX")
echo "working in $work"
cp -r "$here/inputs/." "$work"

# Basis: form orbitals from each HFD input in turn.
cd "$work/basis"
for i in 1 2 3; do
    cp "HFD$i.INP" HFD.INP
    step hfd "hfd$i.out" "$BIN/hfd"
    cp HFD.DAT "HFD$i.DAT"
done

# bass builds on HFD2.DAT and reads extra orbitals from HFD3.DAT.
# The first run always warns, so repeat until clean as basis.py does.
cp HFD2.DAT HFD.DAT
clean=0
for try in 1 2 3 4 5; do
    input=bass.in step bass bass.out "$BIN/bass"
    if ! grep -q -E 'failed|error|warning' bass.out; then
        clean=1
        break
    fi
done
if [ "$clean" -ne 1 ]; then
    echo 'FAIL: bass still reports warnings after 5 attempts'
    exit 1
fi

# CI: one run for each parity.
status=0
for dir in even0 odd1; do
    cd "$work/$dir"
    cp ../basis/HFD.DAT .
    step pbasc pbasc.out mpirun --oversubscribe -n "$nprocs" "$BIN/pbasc"
    step pconf pconf.out mpirun --oversubscribe -n "$nprocs" "$BIN/pconf"
    echo "$dir, $nprocs rank(s):"
    python3 "$here/../compare.py" "$here/ref/$dir/pconf.csv" pconf.csv || status=1
done

# Link the CI files of the odd levels into a pdtm directory, as dtm.py does.
link_odd() {
    for file in CONF.INP CONF.DAT CONF.INT CONF.DET CONF.XIJ CONFSTR.RES; do
        ln -s "../odd1/$file" "$file"
    done
}

# Transition matrix elements from odd levels (CONF) to even levels (CONF1).
cd "$work/tm"
link_odd
ln -s ../even0/CONF.DET CONF1.DET
ln -s ../even0/CONF.XIJ CONF1.XIJ
ln -s ../even0/CONFSTR.RES CONFSTR1.RES
step pdtm pdtm.out mpirun --oversubscribe -n "$nprocs" "$BIN/pdtm"
echo "tm, $nprocs rank(s):"
python3 "$here/../compare.py" "$here/ref/tm/tm.csv" tm.csv || status=1

# g-factors and hyperfine constants of the odd levels.
# The even block has Jm=0, which makes its g-factors and A constants zero.
cd "$work/dm"
link_odd
step pdtm pdtm.out mpirun --oversubscribe -n "$nprocs" "$BIN/pdtm"
for file in GF A_hf B_hf; do
    echo "dm $file, $nprocs rank(s):"
    python3 "$here/../compare.py" "$here/ref/dm/$file.RES" "$file.RES" || status=1
done

# Repeat pconf in the Emu regime, which treats only the first Nemu configurations in full.
for dir in even0 odd1; do
    mkdir "$work/${dir}_emu"
    cd "$work/${dir}_emu"
    cp ../$dir/{CONF.DAT,CONF.INT,CONF.GNT,HFD.DAT,ci.in} .
    sed '/^Kbrt=/a Nemu= 100' ../$dir/CONF.INP > CONF.INP
    step pconf pconf.out mpirun --oversubscribe -n "$nprocs" "$BIN/pconf"
    if ! grep -q 'Regime Emu' pconf.out; then
        echo 'FAIL: pconf did not enter the Emu regime'
        exit 1
    fi
    echo "${dir}_emu, $nprocs rank(s):"
    python3 "$here/../compare.py" "$here/ref/${dir}_emu/pconf.csv" pconf.csv || status=1
done

exit $status
