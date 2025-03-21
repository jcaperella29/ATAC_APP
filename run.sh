#!/bin/bash

APP_NAME="atac-shiny"
IMAGE_TAG="atac-shiny-app:latest"
PORT=8787

echo "🔧 Building $APP_NAME Docker image..."
docker build -t $IMAGE_TAG .

if [ $? -ne 0 ]; then
    echo "❌ Docker build failed. Check logs above."
    exit 1
fi

echo "🚀 Running $APP_NAME on http://localhost:$PORT"
docker run --rm -p $PORT:8080 \
    -v $(pwd):/app \
    $IMAGE_TAG
