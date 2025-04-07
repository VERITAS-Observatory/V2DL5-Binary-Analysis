#!/bin/bash
# spectral plotting using Eventdisplay container and podman
#

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <anasum_file> <output file> <config file>"
    echo
    exit 1
fi

plot_spectra()
{
    ANASUM_FILE="/data/$(basename $1)"
    OUTPUT_FILE="/output/$(basename $2)"
    CONFIG_FILE="/config/$(basename $3)"
    echo "OUTPUT_FILE: $OUTPUT_FILE"
    echo "CONFIG_FILE: $CONFIG_FILE"
    podman run --platform linux/amd64 --rm -it -v "$(pwd)/:/workdir" \
        -v "$(dirname $1):/data" \
        -v "$(dirname $2):/output" \
        -v "$(dirname $3):/config" \
        ghcr.io/veritas-observatory/eventdisplay_v4:main \
        /bin/bash -c \
        "cd /workdir/; root -l -q -b 'spectral_analysis.C(\"$ANASUM_FILE\", \"$CONFIG_FILE\",  \"$OUTPUT_FILE\")';"
}

echo "Running spectral analysis..."
plot_spectra "$1" "$2" "$3"
