#!/bin/bash
# Quick Docker setup script for BigSuR

set -e

echo "=== BigSuR Docker Setup ===" 

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "❌ Docker is not installed. Please install Docker first."
    echo "   See: https://docs.docker.com/get-docker/"
    exit 1
fi

# Check if Docker Compose is installed
if ! command -v docker-compose &> /dev/null; then
    echo "❌ Docker Compose is not installed. Please install Docker Compose."
    echo "   See: https://docs.docker.com/compose/install/"
    exit 1
fi

echo "✓ Docker and Docker Compose are installed"
echo ""

# Build images
echo "📦 Building Docker images..."
docker-compose build

echo ""
echo "✓ Build complete"
echo ""

# Optional: Start services
read -p "Start services now? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "🚀 Starting services..."
    docker-compose up -d
    echo ""
    echo "✓ Services started"
    echo ""
    echo "📊 Service Status:"
    docker-compose ps
    echo ""
    echo "🌐 Web Application: http://localhost:8050"
    echo "📊 Redis: localhost:6379"
    echo ""
    echo "📝 View logs: docker-compose logs -f"
    echo "🛑 Stop services: docker-compose down"
else
    echo "ℹ️  To start services later, run: docker-compose up -d"
fi
