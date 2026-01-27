# Docker Files Summary

## Files Created

### Core Docker Files

1. **Dockerfile** - All-in-one container (optional, not used with docker-compose)
2. **Dockerfile.web** - Optimized web server image (Gunicorn)
3. **Dockerfile.celery** - Optimized Celery worker image
4. **docker-compose.yml** - Main orchestration file (redis, web, celery-worker)
5. **.dockerignore** - Excludes unnecessary files from Docker context

### Compose Overrides

6. **docker-compose.dev.yml** - Development settings (live reload, debug flags)
7. **docker-compose.prod.yml** - Production settings (resource limits, optimizations)

### Configuration

8. **nginx.conf** - Reverse proxy configuration (optional)

### Documentation & Scripts

9. **DOCKER_GUIDE.md** - Comprehensive Docker usage guide
10. **docker-setup.sh** - Quick setup and initialization script

## Quick Reference

### Start Development Environment
```bash
docker-compose up -d
```

### Start Production Environment
```bash
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d
```

### View Logs
```bash
docker-compose logs -f web
docker-compose logs -f celery-worker
docker-compose logs -f redis
```

### Access Container Shell
```bash
docker-compose exec web bash
docker-compose exec celery-worker bash
docker-compose exec redis sh
```

### Stop All Services
```bash
docker-compose down
```

### Clean Up
```bash
# Stop and remove containers, networks
docker-compose down

# Also remove volumes (data)
docker-compose down -v
```

## Architecture

```
┌─────────────────────────────────────────┐
│         Docker Compose Network          │
├─────────────────────────────────────────┤
│                                         │
│  ┌──────────────┐                       │
│  │   Nginx      │ :80, :443             │
│  │  (optional)  │                       │
│  └───────┬──────┘                       │
│          │                              │
│  ┌───────▼──────────────┐               │
│  │  Web Server (port)   │               │
│  │  Gunicorn :8050      │◄─────┐        │
│  │  (4 workers)         │      │        │
│  └───────┬──────────────┘      │        │
│          │                     │        │
│  ┌───────▼───────────────┐  ┌──┴────────┤
│  │  Celery Worker        │  │           │
│  │  (tasks.celery)       │  │  Redis    │
│  │  (4 concurrency)      │  │  :6379    │
│  └───────────────────────┘  │           │
│                             │ (cache,   │
│                             │  broker)  │
│                             └───────────┘
│                                         │
└─────────────────────────────────────────┘

Host Machine:
:8050  ──► Web Server
:6379  ──► Redis (if exposed)
:80    ──► Nginx (if used)
```

## Key Environment Variables

The services communicate via:
- `CELERY_BROKER_URL=redis://redis:6379/0`
- `CELERY_RESULT_BACKEND=redis://redis:6379/0`
- `PYTHONUNBUFFERED=1` (real-time logging)

## Volume Mounts

| Volume | Purpose | Mount Point |
|--------|---------|------------|
| `redis_data` | Redis persistence | `/data` |
| `./src` | Source code | `/app/src` |
| `./data` | Input data files | `/app/data` |
| `./logs` | Application logs | `/app/logs` |
| `./tmp` | Temporary files | `/app/tmp` |

## Health Checks

All services have health checks configured:
- **Redis**: `redis-cli ping`
- **Web**: HTTP request to `http://localhost:8050/`
- **Celery**: Included in container logs

## Networking

All services are on the `bigsur-network` bridge network:
- `redis` ← internal hostname for Redis
- `web` ← internal hostname for web server
- `celery-worker` ← internal hostname for Celery

Services can communicate using these hostnames internally.

## Performance Tuning

### Development (docker-compose.yml)
- 4 Gunicorn workers
- 4 Celery concurrency
- No resource limits

### Production (docker-compose.prod.yml)
- 8 Gunicorn workers
- Resource limits applied
- Restart policies enabled
- Optional Nginx reverse proxy
- Can scale Celery workers: `--scale celery-worker=3`

## Next Steps

1. **Build images**: `docker-compose build`
2. **Start services**: `docker-compose up -d`
3. **Access app**: `http://localhost:8050`
4. **Check logs**: `docker-compose logs -f`
5. **For production**: Add SSL, configure domain, scale workers

## Troubleshooting

See [DOCKER_GUIDE.md](DOCKER_GUIDE.md#troubleshooting) for detailed troubleshooting steps.

Common issues:
- Port already in use: `lsof -i :8050` or `:6379`
- Services won't connect: Check `docker-compose ps` and logs
- Out of memory: Run `docker system prune`
- Database initialization: Check Redis connection and Celery status
