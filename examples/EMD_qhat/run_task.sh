#!/bin/bash
# Run one hydro task of tasks/*.txt and keep only its spectrum file.
#   run_task.sh <task> [threads]
# runJetscape writes the final-state hadrons into the named pipe
# /tmp/fifo/<task>_final_state_hadrons.dat, which jetana reads on the fly, so
# no event file is stored; the result is spec/<task>.spec.
# The work directory is a private copy on /tmp that links to the build, except
# for the files MUSIC writes (music_input, surface*), so concurrent tasks on
# one node do not overwrite each other.
# Environment:
#   XSCAPE       X-SCAPE source tree (default: two levels above this folder)
#   BUILD        X-SCAPE build directory (default: $XSCAPE/build)
#   KEEP_EVENTS  1: also keep the final-state hadrons, gzipped, in events/<task>.dat.gz
#                (about 10-90 MB per task), so observables can be added later
n=$1
W=$(cd "$(dirname "$0")" && pwd)
XSCAPE=${XSCAPE:-$(cd "$W/../.." && pwd)}
B=${BUILD:-$XSCAPE/build}
export OMP_NUM_THREADS=${2:-2}
[ -z "$n" ] && { echo "usage: run_task.sh <task> [threads]"; exit 2; }
[ -f "$W/done/$n" ] && { echo "$n already done"; exit 0; }
mkdir -p "$W/spec" "$W/logs" "$W/done" /tmp/fifo
[ "$KEEP_EVENTS" = 1 ] && mkdir -p "$W/events"
s=$(date +%s)
T=/tmp/emdq_work_$n; rm -rf "$T"; mkdir -p "$T"; cd "$T"
for f in "$B"/*; do b=$(basename "$f")
  case $b in surface*|music_input|*.log|CMake*|Makefile|cmake_install.cmake|src|examples|external_packages) ;;
    *) ln -s "$f" "$b" ;; esac
done
cp "$B/music_input" . 2>/dev/null
f=/tmp/fifo/${n}_final_state_hadrons.dat; g=/tmp/fifo/${n}_ana; rm -f "$f" "$g"; mkfifo "$f"
if [ "$KEEP_EVENTS" = 1 ]; then   # tee the stream into jetana and a gzipped copy
  mkfifo "$g"
  "$W/jetana" "$W/spec/$n.spec" "$g" > "$W/logs/$n.ana.log" 2>&1 & apid=$!
  ( tee "$g" < "$f" | gzip -1 > "$W/events/$n.dat.gz" ) & tpid=$!
else
  "$W/jetana" "$W/spec/$n.spec" "$f" > "$W/logs/$n.ana.log" 2>&1 & apid=$!; tpid=$apid
fi
./runJetscape "$W/configs/$n.xml" "$XSCAPE/config/jetscape_main.xml" > "$T/run.log" 2>&1; rc=$?
kill -0 $tpid 2>/dev/null && timeout 10 sh -c ": > $f"   # release the reader if JETSCAPE died early
wait $tpid; wait $apid; arc=$?
tail -200 "$T/run.log" > "$W/logs/$n.log"
cd /; rm -rf "$T" "$f" "$g"
echo "$n rc=$rc ana=$arc $(( $(date +%s) - s ))s $(hostname -s)" >> "$W/logs/progress.txt"
[ $rc -eq 0 ] && [ $arc -eq 0 ] && touch "$W/done/$n"
