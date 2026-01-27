#!/usr/bin/env zsh
set -euo pipefail

# Minimal startup script for local development
# - Creates `run/` and `logs/`
# - Starts Redis (if not already listening on port 6379)
# - Starts Celery workers with pidfiles under `run/`
# - Starts Gunicorn with pidfile under `run/`

PYTHONPATH=./src
export PYTHONPATH

RUN_DIR=./run
LOG_DIR=./logs
DEBUG_MODE=0

# If invoked with --debug or -d, run the web app in foreground (dev mode)
if [[ "${1:-}" == "--debug" || "${1:-}" == "-d" ]]; then
  DEBUG_MODE=1
  echo "[start] debug mode enabled: running app in foreground"
fi
REDIS_PORT=6379
CELERY_APP="tasks.celery:task_queue"

mkdir -p "$RUN_DIR" "$LOG_DIR" ./redis_data

is_listening() {
  lsof -nP -iTCP:$1 -sTCP:LISTEN >/dev/null 2>&1
}

echo "[start] ensuring runtime directories exist: $RUN_DIR, $LOG_DIR"

# Generate redis.conf with absolute paths
python3 scripts/generate_redis_conf.py

if is_listening $REDIS_PORT; then
  echo "[start] redis already listening on port $REDIS_PORT"
else
  if command -v redis-server >/dev/null 2>&1; then
    echo "[start] starting redis-server on port $REDIS_PORT with redis.conf"
    redis-server redis.conf
    sleep 0.5
    if is_listening $REDIS_PORT; then
      echo "[start] redis started"
    else
      echo "[start] failed to start redis; please start it manually or check $LOG_DIR/redis.log"
    fi
  else
    echo "[start] redis-server not found. Install via brew or run redis manually."
  fi
fi

echo "[start] checking for existing celery pidfiles in $RUN_DIR"
setopt null_glob
for p in $RUN_DIR/celery-*.pid; do
  if [[ -e "$p" ]]; then
    pid=$(cat "$p" 2>/dev/null || true)
    if [[ -n "$pid" ]] && kill -0 "$pid" >/dev/null 2>&1; then
      echo "[start] celery pidfile $p exists and process $pid is running; skipping new celery start"
      CELERY_SKIP=1
    else
      echo "[start] removing stale pidfile $p"
      rm -f "$p"
    fi
  fi
done

# Also check system-wide celery pidfiles (e.g., /var/run/celery) to avoid
# accidentally starting duplicate workers created by other launchers.
if [ -d /var/run/celery ]; then
  for p in /var/run/celery/*.pid; do
    if [[ -e "$p" ]]; then
      pid=$(cat "$p" 2>/dev/null || true)
      if [[ -n "$pid" ]] && kill -0 "$pid" >/dev/null 2>&1; then
        echo "[start] found system celery pidfile $p with running pid $pid; skipping new celery start"
        CELERY_SKIP=1
      else
        echo "[start] removing stale system pidfile $p"
        rm -f "$p"
      fi
    fi
  done
fi

if [[ "${CELERY_SKIP:-0}" == "1" ]]; then
  echo "[start] celery already running according to pidfiles; not starting new workers"
else
  echo "[start] launching celery multi (workers will write pidfiles to $RUN_DIR and logs to $LOG_DIR)"
  celery multi start worker -A $CELERY_APP --pidfile="$RUN_DIR/celery-%n.pid" --logfile="$LOG_DIR/celery-%n.log" --loglevel=info
fi

# Start the web server. In debug mode, run the app in foreground (so ipdb works).
if [[ "$DEBUG_MODE" -eq 1 ]]; then
  echo "[start] Running web app in foreground for debug"
  echo "[start] Logs: $LOG_DIR"
  PYTHON_CMD=$(command -v python || command -v python3 || echo python)
  exec "$PYTHON_CMD" src/index.py
else
  # Start Gunicorn (use the active Python interpreter so module-installed gunicorn is used)
  GUNICORN_PID="$RUN_DIR/gunicorn.pid"
  if [[ -e "$GUNICORN_PID" ]]; then
    pid=$(cat "$GUNICORN_PID" 2>/dev/null || true)
    if [[ -n "$pid" ]] && kill -0 "$pid" >/dev/null 2>&1; then
      echo "[start] gunicorn already running with pid $pid; skipping"
    else
      echo "[start] removing stale gunicorn pidfile"
      rm -f "$GUNICORN_PID"
    fi
  fi

  echo "[start] launching gunicorn via 'python -m gunicorn'"
  # Prefer the python on PATH (the active venv/conda env). Fallback to 'python3' then 'python'.
  PYTHON_CMD=$(command -v python || command -v python3 || echo python)
  nohup "$PYTHON_CMD" -m gunicorn --pythonpath src -w 4 -b 0.0.0.0:8050 index:server \
    --pid "$GUNICORN_PID" --access-logfile "$LOG_DIR/gunicorn-access.log" --error-logfile "$LOG_DIR/gunicorn-error.log" > "$LOG_DIR/gunicorn-stdout.log" 2>&1 &

  echo "[start] done. Logs: $LOG_DIR"
fi
