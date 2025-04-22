#!/bin/bash
# Light-curve analysis using Eventdisplay container and podman
#

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <anasum_file> <time_interval_file>"
    echo "  (use 'RUNWISE' for time interval file for run-wise analysis)"
    echo
    exit 1
fi

light_curves()
{
    ANASUM_FILE="/data/$(basename "$1")"
    TIME_FILE="/time_data/$(basename "$2")"
    podman run --platform linux/amd64 --rm -it -v "$(pwd)/:/workdir" \
        -v "$(dirname "$1"):/data" \
        -v "$(dirname "$2"):/time_data" \
        ghcr.io/veritas-observatory/eventdisplay_v4:main \
        /bin/bash -c \
        "cd /workdir/; root -l -q -b 'light_curve_analysis.C(\"$ANASUM_FILE\", \"$TIME_FILE\", $3, $4, $5)';"
}

echo "Running light-curve analysis..."
ETRESH="0.3"
MINSIG="-5."
MINEVT="-10."
MINSIG="2."
MINEVT="2."
light_curves "$1" "$2" "$ETRESH" "$MINSIG" "$MINEVT"
