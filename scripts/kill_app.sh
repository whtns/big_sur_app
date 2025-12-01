#!/usr/bin/env bash
# Safe helper to stop BigSuR processes (gunicorn, celery, redis) started by the app.
# Usage: run from project root: `bash scripts/kill_app.sh`

set -o pipefail
# do not exit immediately on errors; try to clean up as much as possible
shopt -s nullglob

RUN_DIR="./run"
LOG_DIR="./logs"

echo "[kill_app] run dir: $RUN_DIR, logs: $LOG_DIR"

if [[ ! -d "$RUN_DIR" ]]; then
  echo "[kill_app] run directory not found: $RUN_DIR (continuing)"
fi

echo "[kill_app] Stopping processes referenced by pidfiles in $RUN_DIR"
for pidfile in "$RUN_DIR"/*.pid; do
  [ -e "$pidfile" ] || continue
  pid=$(cat "$pidfile" 2>/dev/null || true)
  echo "[kill_app] pidfile $pidfile -> pid=$pid"
  if [[ -z "$pid" ]]; then
    echo "[kill_app] empty pidfile, removing $pidfile"
    rm -f "$pidfile"
    continue
  fi
  if kill -0 "$pid" >/dev/null 2>&1; then
    echo "[kill_app] sending TERM to $pid"
    kill "$pid" >/dev/null 2>&1 || true
    # wait up to 5s
    for i in 1 2 3 4 5; do
      if kill -0 "$pid" >/dev/null 2>&1; then
        sleep 1
      else
        break
      fi
    done
    if kill -0 "$pid" >/dev/null 2>&1; then
      echo "[kill_app] pid $pid did not exit, sending KILL"
      kill -9 "$pid" >/dev/null 2>&1 || true
    else
      echo "[kill_app] pid $pid terminated"
    fi
  else
    echo "[kill_app] pid $pid not running"
  fi
  rm -f "$pidfile"
done

echo "[kill_app] Checking common ports (8050 for gunicorn, 6379 for redis)"
kill_port() {
  port="$1"
  proto=${2:-TCP}
  pids=$(lsof -t -i${proto}:${port} -sTCP:LISTEN 2>/dev/null || true)
  if [[ -n "$pids" ]]; then
    echo "[kill_app] processes listening on port $port: $pids"
    for p in $pids; do
      if kill -0 "$p" >/dev/null 2>&1; then
        echo "[kill_app] sending TERM to $p (port $port)"
        kill "$p" >/dev/null 2>&1 || true
      fi
    done
    sleep 1
    for p in $pids; do
      if kill -0 "$p" >/dev/null 2>&1; then
        echo "[kill_app] forcing KILL $p (port $port)"
        kill -9 "$p" >/dev/null 2>&1 || true
      fi
    done
  else
    echo "[kill_app] no listeners on port $port"
  fi
}

kill_port 8050
kill_port 6379

echo "[kill_app] Attempting graceful redis shutdown via redis-cli (if available)"
if command -v redis-cli >/dev/null 2>&1; then
  if redis-cli -p 6379 ping >/dev/null 2>&1; then
    echo "[kill_app] redis-cli reports redis running; issuing shutdown"
    redis-cli -p 6379 shutdown || true
  else
    echo "[kill_app] redis-cli: redis not responding on 6379"
  fi
else
  echo "[kill_app] redis-cli not found; skipped graceful redis shutdown"
fi

echo "[kill_app] Killing common process names (gunicorn, celery, redis-server)"
for name in gunicorn celery redis-server; do
  pids=$(pgrep -f "$name" || true)
  if [[ -n "$pids" ]]; then
    echo "[kill_app] processes matching '$name': $pids"
    # try graceful pkill first
    pkill -f "$name" >/dev/null 2>&1 || true
    sleep 1
    # force remaining
    pids_remain=$(pgrep -f "$name" || true)
    if [[ -n "$pids_remain" ]]; then
      echo "[kill_app] forcing kill of: $pids_remain"
      pkill -9 -f "$name" >/dev/null 2>&1 || true
    fi
  else
    echo "[kill_app] no processes matching '$name'"
  fi
done

echo "[kill_app] Final cleanup: remove any remaining pidfiles in $RUN_DIR"
rm -f "$RUN_DIR"/*.pid || true

echo "[kill_app] Summary:"
echo " - Ports:"
lsof -nP -iTCP:8050 -sTCP:LISTEN 2>/dev/null || echo "  8050 clear"
lsof -nP -iTCP:6379 -sTCP:LISTEN 2>/dev/null || echo "  6379 clear"
echo " - Processes (pgrep checks):"
pgrep -af gunicorn || echo "  no gunicorn"
pgrep -af celery || echo "  no celery"
pgrep -af redis-server || echo "  no redis-server"

echo "[kill_app] Done. If you want this script executable: chmod +x scripts/kill_app.sh"
