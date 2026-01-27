#!/usr/bin/env python3
"""
Generate redis.conf with absolute paths based on current directory.
This ensures Redis can start from any working directory.
"""
import os
import sys

# Get the directory where this script is located (should be in scripts/)
script_dir = os.path.dirname(os.path.abspath(__file__))
# Go up one level to the app root
app_root = os.path.dirname(script_dir) if os.path.basename(script_dir) == 'scripts' else script_dir

# Create directories if they don't exist
os.makedirs(os.path.join(app_root, 'run'), exist_ok=True)
os.makedirs(os.path.join(app_root, 'logs'), exist_ok=True)
os.makedirs(os.path.join(app_root, 'redis_data'), exist_ok=True)

redis_conf = f"""# BigSuR Redis Configuration
# Auto-generated with absolute paths

# Network
bind 127.0.0.1
port 6379
timeout 0
tcp-keepalive 300

# General
daemonize yes
pidfile {os.path.join(app_root, 'run', 'redis.pid')}
loglevel notice
logfile {os.path.join(app_root, 'logs', 'redis.log')}

# Data persistence
dir {os.path.join(app_root, 'redis_data')}
dbfilename dump.rdb
save 900 1
save 300 10
save 60 10000

# Memory
maxmemory 256mb
maxmemory-policy allkeys-lru

# Append only file
appendonly no

# Slow log
slowlog-log-slower-than 10000
slowlog-max-len 128
"""

config_path = os.path.join(app_root, 'redis.conf')
with open(config_path, 'w') as f:
    f.write(redis_conf)

print(f"Generated {config_path}")
print(f"App root: {app_root}")
