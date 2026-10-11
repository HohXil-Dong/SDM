#!/usr/bin/env bash
# Plot SDM horizontal/vertical GPS fits with GMT 6. Run from the data directory:
#   bash plot_gps_fit.sh -e all_ew.dat -n all_ns.dat -u up.dat -o gps_fit
#   bash plot_gps_fit.sh -u up.dat -o gps_up
# Input: latitude longitude observation residual prediction correction.
# E/N/[U] must have the same station order and use the same displacement unit.

# ---- Configuration: paths are relative to the working directory ----
EAST="all_ew.dat"
NORTH="all_ns.dat"
UP=""                       # Optional U file.
OUT="gps_fit"               # Prefix: OUT_horizontal.pdf, OUT_vertical.pdf, OUT_merged.dat.
R="97.2/99.8/33.8/35.4"         # GMT region: west/east/south/north.
J="M14c"                    # Mercator projection, map width in cm.
SCALE="0.02c"               # Plot length per input displacement unit.
REF="50"                    # Reference-arrow magnitude in input units.
UNIT="cm"                   # Label only; no conversion (Maduo datunit = 0.01 m).
FAULT_PATCHES="slip_model_patches.dat"  # 12-column UL/UR/LR/LL corner table; "" disables.
FAULT_TRACE="slip_model_trace.gmt"      # Two-column segment traces; "" disables.

set -euo pipefail
# Keep concurrent plotting scripts in separate GMT modern sessions.
export GMT_SESSION_NAME=$$

usage() {
    cat <<'EOF'
Usage: bash plot_gps_fit.sh [-e EAST -n NORTH] [-u UP] [-o PREFIX] [-h]
  -e FILE    East SDM fit file (must be paired with -n).
  -n FILE    North SDM fit file (must be paired with -e).
  -u FILE    Up SDM fit file; may be used alone.
  -o PREFIX  Output prefix without .pdf/.png (default: gps_fit).
  -h         Show this help.
With no component options, use the configured E/N/[U] files.
With any component option, use only the explicitly supplied components.
Input: one header line, then latitude longitude observation residual prediction correction.
Files share station coordinates/order and displacement units. Paths are relative to cwd.
Outputs: PREFIX_horizontal.pdf and/or PREFIX_vertical.pdf, plus PREFIX_merged.dat.
Vertical arrows: upward = uplift; downward = subsidence. Units/range/scale use configuration above.
EOF
}

read_component() {
    # Output: lat lon observation fit residual. Keep coordinate tokens for pairing.
    # SDM column 4 excludes correction: fit = column 5 + column 6.
    awk '
    function fail(message) { print "Error: " FILENAME ": " message > "/dev/stderr"; exit 1 }
    NR==1 {next}
    {
        if (NF != 6) fail("expected six columns at input line " NR)
        for (i=1; i<=NF; i++)
            if ($i !~ /^[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?$/)
                fail("nonnumeric value at input line " NR)
        fit=$5+$6
        printf "%s %s %.17g %.17g %.17g\n", $1,$2,$3,fit,$3-fit
    }
    END { if (NR <= 1) fail("no station data") }
    ' "$1"
}

vertical_vectors() {
    # GMT screen angles: 90 = up, 270 = down; zero displacement has no arrow.
    # The explicit c suffix gives the same cm/unit scale as horizontal arrows.
    awk -v column="$2" -v scale="${SCALE%c}" 'NR>1 && $column!=0 {
        printf "%.5f %.5f %d %.10gc\n", $1,$2,($column>0 ? 90 : 270), \
            sqrt($column*$column)*scale
    }' "$1"
}

EXPLICIT_COMPONENTS=0
EAST_OPT=""; NORTH_OPT=""; UP_OPT=""
while getopts ':e:n:u:o:h' option; do
    case "$option" in
        e) EAST_OPT="$OPTARG"; EXPLICIT_COMPONENTS=1 ;;
        n) NORTH_OPT="$OPTARG"; EXPLICIT_COMPONENTS=1 ;;
        u) UP_OPT="$OPTARG"; EXPLICIT_COMPONENTS=1 ;;
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
if (( EXPLICIT_COMPONENTS )); then
    EAST="$EAST_OPT"; NORTH="$NORTH_OPT"; UP="$UP_OPT"
fi
if [[ -n "$EAST" && -z "$NORTH" ]] || [[ -z "$EAST" && -n "$NORTH" ]]; then
    echo "Error: -e and -n must be supplied together." >&2
    exit 2
fi
files=()
HAS_HORIZONTAL=0
HAS_UP=0
if [[ -n "$EAST" ]]; then
    files+=("$EAST" "$NORTH")
    HAS_HORIZONTAL=1
fi
if [[ -n "$UP" ]]; then
    files+=("$UP")
    HAS_UP=1
fi
(( ${#files[@]} > 0 )) || { echo "Error: supply an E/N pair or a U file." >&2; exit 2; }
command -v gmt >/dev/null || { echo "Error: GMT 6 is required." >&2; exit 1; }
for file in "${files[@]}"; do
    [[ -r "$file" ]] || { echo "Error: cannot read $file" >&2; exit 1; }
done
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
component_files=()
for i in "${!files[@]}"; do
    component_files+=("$TMP/component_$i.dat")
    read_component "${files[$i]}" > "$TMP/component_$i.dat"
done

# Pair rows explicitly and fail if component coordinates/order do not match.
# Prepared component rows each have five columns. All units remain unchanged.
paste "${component_files[@]}" | awk -v has_horizontal="$HAS_HORIZONTAL" -v has_up="$HAS_UP" '
function fail(message) { print "Error: " message > "/dev/stderr"; exit 1 }
BEGIN { print "# lon lat E_obs N_obs U_obs E_fit N_fit U_fit E_res N_res U_res" }
{
    if (NF != (has_horizontal ? 10 : 0) + (has_up ? 5 : 0))
        fail("component station counts differ at data row " NR)
    for (i=6; i<=NF; i+=5)
        if ($1 != $i || $2 != $(i+1))
            fail("component station coordinates/order differ at data row " NR)
    if (seen[$1 "," $2]++)
        fail("duplicate station coordinate at data row " NR)
    e_obs=n_obs=e_fit=n_fit=e_res=n_res="NaN"
    u_obs=u_fit=u_res="NaN"
    if (has_horizontal) {
        e_obs=sprintf("%.8g", $3); n_obs=sprintf("%.8g", $8)
        e_fit=sprintf("%.8g", $4); n_fit=sprintf("%.8g", $9)
        e_res=sprintf("%.8g", $5); n_res=sprintf("%.8g", $10)
    }
    if (has_up) {
        offset=has_horizontal ? 10 : 0
        u_obs=sprintf("%.8g", $(offset+3))
        u_fit=sprintf("%.8g", $(offset+4))
        u_res=sprintf("%.8g", $(offset+5))
    }
    printf "%.5f %.5f %s %s %s %s %s %s %s %s %s\n", \
        $2,$1,e_obs,n_obs,u_obs,e_fit,n_fit,u_fit,e_res,n_res,u_res
}
END { if (NR == 0) fail("no station data") }
' > "$TMP/merged.dat"
mkdir -p "$(dirname "$OUT")"
MERGED="${OUT}_merged.dat"
mv "$TMP/merged.dat" "$MERGED"

# GMT velo reads lon lat E N sigma_E sigma_N correlation.
# SDM output has no uncertainties, so only displacement arrows are drawn.
if (( HAS_HORIZONTAL )); then
    awk 'NR>1 {print $1,$2,$3,$4,0,0,0}' "$MERGED" > "$TMP/observation.xy"
    awk 'NR>1 {print $1,$2,$6,$7,0,0,0}' "$MERGED" > "$TMP/fit.xy"
fi
if (( HAS_UP )); then
    vertical_vectors "$MERGED" 5 > "$TMP/up_observation.xy"
    vertical_vectors "$MERGED" 8 > "$TMP/up_fit.xy"
fi
read -r XREF YREF < <(awk -v region="$R" 'BEGIN {
    split(region, r, "/"); print r[1]+0.08*(r[2]-r[1]), r[3]+0.08*(r[4]-r[3])
}')

plot_gps_map() {
    local component="$1"
    # Render to a temporary prefix: GMT's session metadata alters spaces in prefixes.
    gmt begin "$TMP/gps_${component}" pdf
    gmt set MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 9p FONT_LABEL 10p \
        FORMAT_GEO_MAP ddd:mmF
    gmt basemap -R"$R" -J"$J" -Baf -BWSne
    if [[ -n "$FAULT_PATCHES" ]]; then
        gmt plot "$TMP/patches.xy" -W0.25p,gray60
    fi
    awk 'NR>1 {print $1,$2}' "$MERGED" | \
        gmt plot -Sc0.06c -Gwhite -W0.4p,gray50
    if [[ "$component" == horizontal ]]; then
        gmt velo "$TMP/observation.xy" -Se"$SCALE/0+f0" \
            -A0.24c+e+p0.9p -W1.2p,gray20 -Ggray20
        gmt velo "$TMP/fit.xy" -Se"$SCALE/0+f0" \
            -A0.18c+e+p0.7p -W0.8p,#D55E00 -G#D55E00
    else
        if [[ -s "$TMP/up_observation.xy" ]]; then
            gmt plot "$TMP/up_observation.xy" -Sv0.24c+e+n0.24c/0 -W1.2p,gray20 -Ggray20
        fi
        if [[ -s "$TMP/up_fit.xy" ]]; then
            gmt plot "$TMP/up_fit.xy" -Sv0.18c+e+n0.18c/0 -W0.8p,#D55E00 -G#D55E00
        fi
    fi
    if [[ -n "$FAULT_TRACE" ]]; then
        gmt plot "$FAULT_TRACE" -W1p,black
    fi

    if [[ "$component" == horizontal ]]; then
        printf '%s %s %s 0 0 0 0\n' "$XREF" "$YREF" "$REF" | \
            gmt velo -Se"$SCALE/0+f0" -A0.22c+e+p0.8p -W1p,gray20 -Ggray20
        printf '%s %s %s %s\n' "$XREF" "$YREF" "$REF" "$UNIT" | \
            gmt text -F+f8p,Helvetica,gray20+jBL -D0/0.24c
        gmt legend -DjTR+w2.8c+o0.2c -F+gwhite <<'EOF'
S 0.3c - 0.5c - 1.2p,gray20 0.7c Observation
S 0.3c - 0.5c - 0.8p,#D55E00 0.7c Fit
EOF
    else
        awk -v x="$XREF" -v y="$YREF" -v ref="$REF" -v scale="${SCALE%c}" \
            'BEGIN {printf "%.5f %.5f 90 %.10gc\n", x,y,ref*scale}' | \
            gmt plot -Sv0.22c+e -W1p,gray20 -Ggray20
        printf '%s %s %s %s\n' "$XREF" "$YREF" "$REF" "$UNIT" | \
            gmt text -F+f8p,Helvetica,gray20+jBL -D0.15c/0.12c
        gmt legend -DjTR+w5.8c+o0.2c -F+gwhite <<'EOF'
S 0.3c - 0.5c - 1.2p,gray20 0.7c Observation
S 0.3c - 0.5c - 0.8p,#D55E00 0.7c Fit
L 8p Helvetica L Vertical displacement: upward = uplift
L 8p Helvetica L Downward = subsidence
EOF
    fi
    gmt end
    mv "$TMP/gps_${component}.pdf" "${OUT}_${component}.pdf"
    echo "Figure: ${OUT}_${component}.pdf"
}

if (( HAS_HORIZONTAL )); then
    plot_gps_map horizontal
fi
if (( HAS_UP )); then
    plot_gps_map vertical
fi
echo "Merged data: $MERGED"
