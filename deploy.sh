#!/bin/bash

# Check if podman is installed
if ! command -v podman &> /dev/null; then
    echo "Error: podman is not installed or not in PATH."
    exit 1
fi

# Parse command line arguments
if [ $# -eq 0 ]; then
    echo "Error: DATA_PATH is required."
    echo "Usage: $0 DATA_PATH"
    echo "Use --help for more information."
    exit 1
fi

if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    echo "Usage: $0 DATA_PATH"
    echo "Deploy the STOP gAPP container."
    echo ""
    echo "Arguments:"
    echo "  DATA_PATH  Path to the data directory (required)"
    echo ""
    echo "Options:"
    echo "  --help, -h  Show this help message"
    exit 0
fi

DATA_PATH="$1"

# Check if the data path is a valid directory
if [ ! -d "$DATA_PATH" ]; then
    echo "Error: '$DATA_PATH' does not exist or is not a directory."
    exit 1
fi

# Pull the latest image from GitHub Container Registry
if ! podman pull ghcr.io/Computational-Rare-Disease-Genomics-WHG/the_stop_gapp:latest; then
    echo "Error: Failed to pull the image."
    exit 1
fi

# Run the container
if ! podman run -d \
  --name stopgapp \
  --restart unless-stopped \
  --hostname stopgapp \
  -p 3838:3838 \
  -v "$DATA_PATH":/home/shiny-app/data \
  -e SHINY_PROD=1 \
  ghcr.io/Computational-Rare-Disease-Genomics-WHG/the_stop_gapp:latest \
  Rscript /home/shiny-app/app-prod.R; then
    echo "Error: Failed to run the container."
    exit 1
fi

echo "Container 'stopgapp' started successfully."