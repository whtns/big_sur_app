# Dockerization Guide for BigSuR

## Overview

This Docker setup containerizes the entire BigSuR application with three separate services:
- **Redis**: Message broker and cache
- **Celery Worker**: Background task processor
- **Web Server**: Gunicorn-based web application

## Quick Start

### Prerequisites
- Docker (version 20.10+)
- Docker Compose (version 1.29+)

### Running the Application

```bash
# Navigate to project directory
cd /path/to/bigsur_app

# Build and start all services
docker-compose up -d

# View logs
docker-compose logs -f

# Check service status
docker-compose ps
```

The application will be available at `http://localhost:8050`

## Services Details

### Redis Service
- **Container**: `bigsur-redis`
- **Image**: `redis:7-alpine`
- **Port**: `6379` (internal: `redis:6379`)
- **Purpose**: Caching and task queue broker
- **Volume**: `redis_data` (persists to disk with AOF)

### Celery Worker Service
- **Container**: `bigsur-celery-worker`
- **Base Image**: Python 3.10 slim
- **Purpose**: Processes background tasks (email, long-running computations)
- **Concurrency**: 4 workers
- **Connected to**: Redis

### Web Service
- **Container**: `bigsur-web`
- **Base Image**: Python 3.10 slim
- **Port**: `8050`
- **Purpose**: Main web application (Dash + Flask)
- **Server**: Gunicorn (4 workers)
- **Connected to**: Redis and Celery

## Common Commands

### View logs for specific service
```bash
docker-compose logs web        # Web server logs
docker-compose logs celery-worker  # Celery logs
docker-compose logs redis      # Redis logs
```

### Access a service's shell
```bash
docker-compose exec web bash
docker-compose exec celery-worker bash
docker-compose exec redis sh
```

### Stop services
```bash
# Stop but keep containers
docker-compose stop

# Stop and remove containers
docker-compose down

# Stop and remove everything including volumes
docker-compose down -v
```

### Rebuild services
```bash
# Rebuild all images
docker-compose build --no-cache

# Rebuild specific service
docker-compose build --no-cache web
```

## Environment Variables

The services use these key environment variables (set in docker-compose.yml):
- `PYTHONUNBUFFERED=1`: Unbuffered Python output for real-time logs
- `CELERY_BROKER_URL=redis://redis:6379/0`: Celery broker connection
- `CELERY_RESULT_BACKEND=redis://redis:6379/0`: Celery result storage

Add more environment variables in `docker-compose.yml` as needed.

## Volumes & Data

### Mounted Volumes
- `redis_data`: Redis persistence data
- `./src`: Application source code (live reloading)
- `./data`: Input data files (h5ad, etc.)
- `./logs`: Application logs
- `./tmp`: Temporary files

### Data Management
```bash
# View volume details
docker volume inspect bigsur_redis_data

# Clear Redis data (resets cache & queued tasks)
docker-compose down -v
```

## Configuration Files

### Dockerfile
- Main image for all-in-one deployment (if needed)
- Not used with docker-compose by default

### Dockerfile.web
- Optimized for web server
- Runs Gunicorn with 4 workers
- Health checks enabled

### Dockerfile.celery
- Optimized for background task processing
- Runs Celery worker with 4 concurrency

### docker-compose.yml
- Orchestrates all services
- Sets up networking and volumes
- Configures health checks and dependencies

## Troubleshooting

### Services won't start
```bash
# Check Docker daemon is running
docker ps

# View detailed logs
docker-compose logs --tail=100 web celery-worker redis

# Check port conflicts (8050, 6379)
lsof -i :8050
lsof -i :6379
```

### Cannot connect to Redis
```bash
# Verify Redis is healthy
docker-compose exec redis redis-cli ping

# Check network connectivity
docker-compose exec web ping redis
```

### Celery tasks not being processed
```bash
# Check Celery worker status
docker-compose logs celery-worker | grep -i error

# Inspect Celery queue
docker-compose exec redis redis-cli
> KEYS celery*
```

### Out of memory/resources
```bash
# Clean up stopped containers and dangling images
docker system prune

# Clean up with volumes included
docker system prune -a --volumes
```

## Production Deployment

For production, consider:

1. **Use environment-specific configs**:
   ```bash
   docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d
   ```

2. **Add reverse proxy** (Nginx/Traefik):
   ```yaml
   # Add to docker-compose.yml
   nginx:
     image: nginx:alpine
     ports:
       - "80:80"
       - "443:443"
     # ... SSL config, reverse proxy setup
   ```

3. **Scale Celery workers**:
   ```bash
   docker-compose up -d --scale celery-worker=3
   ```

4. **Use Docker Secrets** for sensitive data (credentials, API keys)

5. **Enable resource limits** in docker-compose.yml:
   ```yaml
   web:
     deploy:
       resources:
         limits:
           cpus: '1'
           memory: 512M
   ```

6. **Add monitoring** (Prometheus, Grafana, DataDog)

## Network Communication

```
┌─────────────────────────────────────────┐
│         bigsur-network (bridge)         │
├─────────────────────────────────────────┤
│  web ←→ redis  ←→  celery-worker        │
│  web ← exposed to host:8050             │
│  redis ← exposed to host:6379 (opt)     │
└─────────────────────────────────────────┘
```

Services communicate via internal hostnames:
- `redis`: Redis service
- `web`: Web service (only accessible within network)
- `celery-worker`: Celery service

## Additional Notes

- All services restart automatically if they crash (default policy)
- Logs are preserved across restarts
- Source code (`./src`) is mounted for live reloading during development
- Data (`./data`) should be pre-populated on host before running
- Database connections must be configured in `src/configmodule/default_config.py` or via environment variables

## Debugging

Enable debug mode in development:
```bash
# Run web service without Gunicorn (direct Flask dev server)
docker-compose run web python src/index.py
```

## Next Steps

1. **Customize base images** if needed (e.g., use `python:3.10-bullseye` for more tools)
2. **Add SSL/TLS** for production with reverse proxy
3. **Set up CI/CD** to auto-build and push images to registry
4. **Configure persistent databases** (PostgreSQL, MongoDB) if needed
5. **Add Nginx/Traefik** for load balancing and routing
