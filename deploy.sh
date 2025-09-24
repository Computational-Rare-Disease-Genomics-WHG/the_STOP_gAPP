#!/bin/bash

# Default data path
DATA_PATH="./data"

# Parse command line arguments
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    echo "Usage: $0 [DATA_PATH]"
    echo "Deploy the STOP gAPP container."
    echo ""
    echo "Arguments:"
    echo "  DATA_PATH  Path to the data directory (default: ./data)"
    echo ""
    echo "Options:"
    echo "  --help, -h  Show this help message"
    exit 0
fi

if [ -n "$1" ]; then
    DATA_PATH="$1"
fi

# Pull the latest image from GitHub Container Registry
podman pull ghcr.io/Computational-Rare-Disease-Genomics-WHG/the_STOP_gAPP:latest

# Run the container
podman run -d \
  --name stopgapp \
  --restart unless-stopped \
  --hostname stopgapp \
  -p 3838:3838 \
  -v "$DATA_PATH":/home/shiny-app/data \
  -e SHINY_PROD=1 \
  ghcr.io/Computational-Rare-Disease-Genomics-WHG/the_STOP_gAPP:latest \
  Rscript /home/shiny-app/app-prod.R