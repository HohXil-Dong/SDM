#!/usr/bin/env bash
# Plot an SDM scalar InSAR fit (LOS or azimuth) with GMT 6. Run from the data directory:
#   bash plot_insar_fit.sh -i DT106_los.dat -o DT106_los_fit
# Input: latitude longitude observation residual prediction correction.
# Fit = prediction + correction; residual = observation - fit.

# ---- Configuration: paths are relative to the working directory ----
INPUT="DT106_az.dat"
OUT="DT106_az_fit"              # Produces OUT.pdf.
R="97/100/33.8/35.4"          # GMT region: west/east/south/north.
J="M8c"                      # Mercator projection, width of each panel in cm.
POINT_SIZE="0.06c"           # Circle diameter.
UNIT="cm"                    # Label only; no conversion of values or observation sign.
CMAP="vik"
DATA_RANGE="-150/150/10"      # Observation/fit CPT: min/max/step in input units.
RESIDUAL_RANGE="-150/150/5"     # Independent residual CPT: min/max/step.
FAULT_PATCHES="slip_model_patches.dat"  # 12-column UL/UR/LR/LL corner table; "" disables.
FAULT_TRACE="slip_model_trace.gmt"      # Two-column segment traces; "" disables.

set -euo pipefail
# Keep concurrent plotting scripts in separate GMT modern sessions.
export GMT_SESSION_NAME=$$

usage() {
    cat <<'EOF'
Usage: bash plot_insar_fit.sh [-i INPUT] [-o PREFIX] [-h]
  -i FILE    SDM scalar InSAR fit file (default: DT106_az.dat).
  -o PREFIX  Output prefix without .pdf/.png (default: DT106_az_fit).
  -h         Show this help.
Input: one header line, then latitude longitude observation residual prediction correction.
Output: PREFIX.pdf with Observation, Fit and Residual panels. Paths are relative to cwd.
Fit = prediction + correction; residual = observation - fit.
Units/range/color scales use the configuration above; no displacement conversion is applied.
EOF
}

while getopts ':i:o:h' option; do
    case "$option" in
        i) INPUT="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        h) usage; exit 0 ;;
        :) echo "Error: -$OPTARG requires an argument." >&2; exit 2 ;;
        \?) echo "Error: unknown option -$OPTARG; use -h." >&2; exit 2 ;;
    esac
    [[ -n "$OPTARG" ]] || { echo "Error: -$option requires a nonempty argument." >&2; exit 2; }
done
shift "$((OPTIND - 1))"
(( $# == 0 )) || { echo "Error: unexpected positional arguments; use -h." >&2; exit 2; }
case "$OUT" in
    ""|*.[pP][dD][fF]|*.[pP][nN][gG])
        echo "Error: -o expects a nonempty prefix without .pdf/.png." >&2; exit 2 ;;
esac
command -v gmt >/dev/null || { echo "Error: GMT 6 is required." >&2; exit 1; }
[[ -r "$INPUT" ]] || { echo "Error: cannot read $INPUT" >&2; exit 1; }
for file in "$FAULT_PATCHES" "$FAULT_TRACE"; do
    [[ -z "$file" ]] && continue
    [[ -r "$file" ]] || {
        echo "Error: cannot read $file; run python plot_slip.py slip_model.dat first." >&2
        exit 1
    }
done

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
if [[ -n "$FAULT_PATCHES" ]]; then
    # One closed map polygon per patch; depths remain in the source table.
    awk '/^#/ || NF==0 {next}
        NF!=12 {print "Error: expected 12 corner columns at line " NR > "/dev/stderr"; exit 1}
        {print ">"; print $1,$2; print $4,$5; print $7,$8; print $10,$11; print $1,$2}
    ' "$FAULT_PATCHES" > "$TMP/patches.xy"
fi

# Keep original coordinates, units, point order and observation sign; column 4 in
# the SDM output excludes correction, so recompute the inversion residual.
awk '
function fail(message) { print "Error: " message > "/dev/stderr"; exit 1 }
NR==1 {next}
{
    if (NF!=6) fail("expected six columns at input line " NR)
    for (i=1; i<=NF; i++)
        if ($i !~ /^[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?$/)
            fail("nonnumeric value at input line " NR)
    fit=$5+$6
    printf "%.5f %.5f %.10g %.10g %.10g\n", $2,$1,$3,fit,$3-fit
    points++
}
END { if (!points) fail("no InSAR data") }
' "$INPUT" > "$TMP/fit.dat"
mkdir -p "$(dirname "$OUT")"

columns=(2 3 4)              # Zero-based columns: Observation Fit Residual.
titles=(Observation Fit Residual)
cpts=("$TMP/data.cpt" "$TMP/data.cpt" "$TMP/residual.cpt")
labels=("Displacement ($UNIT)" "Displacement ($UNIT)" "Residual ($UNIT)")

# Render to a temporary prefix: GMT's session metadata alters spaces in prefixes.
gmt begin "$TMP/insar_fit" pdf
    gmt set MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 9p FONT_LABEL 10p \
        FONT_TITLE 11p FORMAT_GEO_MAP ddd:mmF
    gmt makecpt -C"$CMAP" -T"$DATA_RANGE" -Z -D -H > "$TMP/data.cpt"
    gmt makecpt -C"$CMAP" -T"$RESIDUAL_RANGE" -Z -D -H > "$TMP/residual.cpt"

    # Let GMT derive the Mercator panel height from R/J; reserve space for
    # titles and colorbars. -D leaves frame drawing to each basemap below.
    gmt subplot begin 1x3 -Fs"${J#M}"/0 -R"$R" -J"$J" \
        -M0.4c/0.2c -Cn0.6c -Cs1.5c -Baf -BWSne -D
    for i in 0 1 2; do
        gmt subplot set 0,"$i"
        # Use the same explicit R/J for every layer and its coordinate frame.
        gmt basemap -R"$R" -J"$J" -Baf -BWSne+t"${titles[$i]}"
        if [[ -n "$FAULT_PATCHES" ]]; then
            gmt plot "$TMP/patches.xy" -R"$R" -J"$J" -W0.25p,gray60
        fi
        gmt plot "$TMP/fit.dat" -R"$R" -J"$J" -i0,1,"${columns[$i]}" \
            -Sc"$POINT_SIZE" -C"${cpts[$i]}"
        if [[ -n "$FAULT_TRACE" ]]; then
            gmt plot "$FAULT_TRACE" -R"$R" -J"$J" -W1p,black
        fi
        gmt colorbar -R"$R" -J"$J" -C"${cpts[$i]}" -DJBC+w5c/0.2c+h+o0/0.7c \
            -Baf -Bx+l"${labels[$i]}"
    done
    gmt subplot end
gmt end
mv "$TMP/insar_fit.pdf" "${OUT}.pdf"

echo "Figure: ${OUT}.pdf"
