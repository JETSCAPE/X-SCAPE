#!/bin/bash
# Run the brick configurations of tasks/bricks.txt (about 80 s each, 70 runs).
#   run_bricks.sh [parallel jobs, default 4]
# Each run writes bricks/<set>/<name>_final_state_partons.dat and a log.
# Environment: XSCAPE, BUILD as in run_task.sh.
W=$(cd "$(dirname "$0")" && pwd)
XSCAPE=${XSCAPE:-$(cd "$W/../.." && pwd)}
B=${BUILD:-$XSCAPE/build}
cd "$B" || exit 1
xargs -P "${1:-4}" -I{} sh -c "./runJetscape '$W/configs/{}.xml' '$XSCAPE/config/jetscape_main.xml' > '$W/bricks/{}.log' 2>&1 && echo '{} ok' || echo '{} FAILED'" < "$W/tasks/bricks.txt"
