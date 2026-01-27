#!/bin/bash

# Usage: ./BigSuR.sh [debug]
# If run with argument 'debug', the web app runs in foreground via `python src/index.py`.

## Ensure directories exist and generate redis.conf
mkdir -p run logs redis_data
python3 scripts/generate_redis_conf.py

## Startup the REDIS server with config file
redis-server redis.conf

## Setup the celery task queue
PYTHONPATH=./src celery multi start worker0 -A tasks.celery:task_queue &

# Start the BigSuR web server (debug mode runs foreground python server)
if [[ "${1:-}" == "debug" || "${1:-}" == "--debug" ]]; then
	echo "[BigSuR] debug mode: running app in foreground"
	python src/index.py
else
	gunicorn -b 0.0.0.0:8050 --workers=4 --threads=2 --worker-class gevent --max-requests=25 --timeout=1440 --pythonpath src index:server
fi