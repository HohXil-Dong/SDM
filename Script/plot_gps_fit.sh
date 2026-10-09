#!/usr/bin/env bash
# Plot SDM horizontal GPS fits with GMT 6. Run from the data directory:
#   bash plot_sdm_gps_fit.sh
# Input: latitude longitude observation residual prediction correction.
# E/N/[U] must have the same station order and use the same displacement unit.

# ---- Configuration: paths are relative to the working directory ----
EAST="tri_ew.dat"
NORTH="tri_ns.dat"
UP=""                       # Optional U file; retained in the table, not mapped.
OUT="gps_fit"               # Produces OUT.pdf and OUT_merged.dat.
R="97.5/99.5/34.0/35.0"         # GMT region: west/east/south/north.
J="M14c"                    # Mercator projection, map width in cm.
SCALE="0.02c"               # Plot length per input displacement unit.
REF="50"                    # Reference-arrow magnitude in input units.
UNIT="cm"                   # Label only; no conversion (Maduo datunit = 0.01 m).

set -euo pipefail

if (( $# )); then
    echo "Edit the configuration at the top, then run without arguments." >&2
    exit 2
fi
command -v gmt >/dev/null || { echo "Error: GMT 6 is required." >&2; exit 1; }
files=("$EAST" "$NORTH")
HAS_UP=0
if [[ -n "$UP" ]]; then
    files+=("$UP")
    HAS_UP=1
fi
for file in "${files[@]}"; do
    [[ -r "$file" ]] || { echo "Error: cannot read $file" >&2; exit 1; }
done

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
for i in "${!files[@]}"; do
    tail -n +2 "${files[$i]}" > "$TMP/component_$i.dat"
done

# Pair rows explicitly and fail if component coordinates/order do not match.
# SDM's fourth column excludes correction; fit = column 5 + column 6,
# and the residual below is observation - fit. All units remain unchanged.
paste "$TMP"/component_*.dat | awk -v has_up="$HAS_UP" '
function fail(message) { print "Error: " message > "/dev/stderr"; exit 1 }
BEGIN { print "# lon lat E_obs N_obs U_obs E_fit N_fit U_fit E_res N_res U_res" }
{
    if (NF != (has_up ? 18 : 12))
        fail("expected six columns per component at data row " NR)
    for (i=1; i<=NF; i++)
        if ($i !~ /^[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?$/)
            fail("nonnumeric value at data row " NR)
    if ($1 != $7 || $2 != $8 || (has_up && ($1 != $13 || $2 != $14)))
        fail("component station coordinates/order differ at data row " NR)
    if (seen[$1 "," $2]++)
        fail("duplicate station coordinate at data row " NR)
    e_fit=$5+$6; n_fit=$11+$12
    u_obs=u_fit=u_res="NaN"
    if (has_up) {
        u_obs=sprintf("%.8g", $15)
        u_fit=sprintf("%.8g", $17+$18)
        u_res=sprintf("%.8g", $15-$17-$18)
    }
    printf "%.5f %.5f %.8g %.8g %s %.8g %.8g %s %.8g %.8g %s\n", \
        $2,$1,$3,$9,u_obs,e_fit,n_fit,u_fit,$3-e_fit,$9-n_fit,u_res
}
END { if (NR == 0) fail("no station data") }
' > "$TMP/merged.dat"
mkdir -p "$(dirname "$OUT")"
MERGED="${OUT}_merged.dat"
mv "$TMP/merged.dat" "$MERGED"

# GMT velo reads lon lat E N sigma_E sigma_N correlation.
# SDM output has no uncertainties, so only displacement arrows are drawn.
awk 'NR>1 {print $1,$2,$3,$4,0,0,0}' "$MERGED" > "$TMP/observed.xy"
awk 'NR>1 {print $1,$2,$6,$7,0,0,0}' "$MERGED" > "$TMP/fitted.xy"
read -r XREF YREF < <(awk -v region="$R" 'BEGIN {
    split(region, r, "/"); print r[1]+0.08*(r[2]-r[1]), r[3]+0.08*(r[4]-r[3])
}')

gmt begin "$OUT" pdf
    gmt set MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 9p FONT_LABEL 10p \
        FORMAT_GEO_MAP ddd:mmF
    gmt basemap -R"$R" -J"$J" -Baf -BWSne
    awk 'NR>1 {print $1,$2}' "$MERGED" | \
        gmt plot -Sc0.06c -Gwhite -W0.4p,gray50
    gmt velo "$TMP/observed.xy" -Se"$SCALE/0+f0" \
        -A0.24c+e+p0.9p -W1.2p,gray20 -Ggray20
    gmt velo "$TMP/fitted.xy" -Se"$SCALE/0+f0" \
        -A0.18c+e+p0.7p -W0.8p,#D55E00 -G#D55E00

    printf '%s %s %s 0 0 0 0\n' "$XREF" "$YREF" "$REF" | \
        gmt velo -Se"$SCALE/0+f0" -A0.22c+e+p0.8p -W1p,gray20 -Ggray20
    printf '%s %s %s %s\n' "$XREF" "$YREF" "$REF" "$UNIT" | \
        gmt text -F+f8p,Helvetica,gray20+jBL -D0/0.24c
    gmt legend -DjTR+w2.8c+o0.2c -F+gwhite <<'EOF'
S 0.3c - 0.5c - 1.2p,gray20 0.7c Observed
S 0.3c - 0.5c - 0.8p,#D55E00 0.7c Fitted
EOF
gmt end

echo "Merged data: $MERGED"
echo "Figure: ${OUT}.pdf"
