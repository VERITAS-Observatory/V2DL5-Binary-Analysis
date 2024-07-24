#!/bin/bash
# Generate run lists from time bins
#

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <anasum_file> <time_interval_file> <write run list (true/false)>"
    echo "  generate run lists from time bins"
    echo
    exit 1
fi

 runlist_from_time_bins()
 {
    ANASUM_FILE="/data/$(basename $1)"
    TIME_FILE="/time_data/$(basename $2)"
    podman run --platform linux/amd64 --rm -it -v "$(pwd)/:/workdir" \
        -v "$(dirname $1):/data" \
        -v "$(dirname $2):/time_data" \
        ghcr.io/veritas-observatory/eventdisplay_v4:main \
        /bin/bash -c \
        "cd /workdir/; root -l -q -b 'runlist_from_time_bins.C(\"$ANASUM_FILE\", \"$TIME_FILE\", $3)';"
 }

 echo "Generate run lists..."
runlist_from_time_bins "$1" "$2" "$3"
