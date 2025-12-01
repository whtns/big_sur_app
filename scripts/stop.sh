#!/usr/bin/env zsh
set -euo pipefail

# Minimal stop script to cleanly stop processes started by `start.sh`.
# - Kills processes referenced by pidfiles in `run/`
# - Attempts to shutdown Redis via redis-cli (if available)

RUN_DIR=./run
LOG_DIR=./logs

echo "[stop] stopping services and cleaning pidfiles in $RUN_DIR"

if [[ -d "$RUN_DIR" ]]; then
  for pidfile in $RUN_DIR/*.pid; do
    if [[ -e "$pidfile" ]]; then
      pid=$(cat "$pidfile" 2>/dev/null || true)
      if [[ -n "$pid" ]]; then
        if kill -0 "$pid" >/dev/null 2>&1; then
          echo "[stop] sending TERM to pid $pid (from $pidfile)"
          kill "$pid"
          # wait up to 5s
          for i in {1..5}; do
            if kill -0 "$pid" >/dev/null 2>&1; then
              sleep 1
            else
              break
            fi
          done
          if kill -0 "$pid" >/dev/null 2>&1; then
            echo "[stop] pid $pid did not exit; forcing KILL"
            kill -9 "$pid" || true
          fi
        else
          echo "[stop] pid $pid (from $pidfile) not running; removing stale pidfile"
        fi
      fi
      rm -f "$pidfile"
    fi
  done
else
  echo "[stop] no run directory found at $RUN_DIR"
fi

# Try to shutdown redis cleanly
if command -v redis-cli >/dev/null 2>&1; then
  if redis-cli -p 6379 ping >/dev/null 2>&1; then
    echo "[stop] shutting down redis on port 6379"
    redis-cli -p 6379 shutdown || true
  else
    echo "[stop] redis not responding on 6379"
  fi
else
  echo "[stop] redis-cli not found; you may need to stop redis manually"
fi

echo "[stop] done."
