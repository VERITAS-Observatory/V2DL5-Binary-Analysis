#!/bin/bash
# Nightly time bins using Eventdisplay container and podman
#

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <anasum_file>"
    echo
    exit 1
fi

nightly_time_bins()
{
    ANASUM_FILE="/data/$(basename $1)"
    podman run --platform linux/amd64 --rm -it -v "$(pwd)/:/workdir" \
        -v "$(dirname $1):/data" \
        ghcr.io/veritas-observatory/eventdisplay_v4:main \
        /bin/bash -c \
        "cd /workdir/; root -l -q -b 'nightly_time_bins.C(\"$ANASUM_FILE\")';"
}

echo "Determining nightly time bins..."
nightly_time_bins "$1" "$2"
