#!/bin/bash
set -euo pipefail

# Deploy StopGAPP using rootless Podman + systemd user unit (no quadlets, no compose).

# Configuration
IMAGE="ghcr.io/computational-rare-disease-genomics-whg/the_stop_gapp:latest"
CTR_NAME="stopgapp"
PORT=3838
UNIT_FILE="container-${CTR_NAME}.service"
SYSTEMD_USER_DIR="${HOME}/.config/systemd/user"

# Parse arguments
DATA_PATH=""

usage(){
  echo "Usage: $0 --data-path PATH"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-path)
      DATA_PATH="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      ;;
  esac
done

# Validate
command -v podman >/dev/null || { echo "Error: podman not found"; exit 1; }
if [[ -z "$DATA_PATH" || ! -d "$DATA_PATH" ]]; then
  echo "Error: invalid data-path: $DATA_PATH"
  exit 1
fi

mkdir -p "$SYSTEMD_USER_DIR"

echo ">> Pulling image $IMAGE"
podman pull "$IMAGE"

# Stop old systemd unit if exists
echo ">> Stopping old systemd unit, if any"
systemctl --user stop "$UNIT_FILE" 2>/dev/null || true

# Remove old container
echo ">> Removing existing container (if any)"
podman rm -f "$CTR_NAME" 2>/dev/null || true

# Run new container
echo ">> Running container $CTR_NAME"
podman run -d \
  --replace \
  --name "$CTR_NAME" \
  --restart unless-stopped \
  --hostname "$CTR_NAME" \
  -p "$PORT:$PORT" \
  -v "${DATA_PATH}:/home/shiny-app/data:Z" \
  -e SHINY_PROD=1 \
  "$IMAGE" \
  Rscript /home/shiny-app/app-prod.R

# Generate systemd unit
echo ">> Generating systemd unit"
podman generate systemd \
  --name "$CTR_NAME" \
  --new \
  --files \
  --restart-policy always

# Move the generated unit file to user units directory
echo ">> Installing systemd user unit"
mv -f "$UNIT_FILE" "$SYSTEMD_USER_DIR/$UNIT_FILE"

# Reload systemd user daemon
echo ">> Reloading user systemd"
systemctl --user daemon-reload

# Enable and start
echo ">> Enabling + starting $UNIT_FILE"
systemctl --user enable --now "$UNIT_FILE"

# --- Check lingering status ---
LINGER_STATUS=$(loginctl show-user "$USER" | awk -F= '/Linger/ {print $2}')

if [[ "$LINGER_STATUS" != "yes" ]]; then
  echo ""
  echo "⚠️  Lingering is not enabled for user '$USER'."
  echo "    Without lingering, the container will stop when you log out or on reboot."
  echo "    To fix this, run:"
  echo "      loginctl enable-linger $USER"
  echo ""
else
  echo "✅ Lingering is enabled for user '$USER' (services persist after logout/reboot)."
fi

cat <<EOF

✅ Deployed StopGAPP under systemd --user

Service unit: ${UNIT_FILE}
Container:    ${CTR_NAME}
Image:        ${IMAGE}

Commands:
  systemctl --user status ${UNIT_FILE}
  journalctl --user -u ${UNIT_FILE} -f
  podman logs -f ${CTR_NAME}
  podman pull ${IMAGE} && systemctl --user restart ${UNIT_FILE}

EOF
