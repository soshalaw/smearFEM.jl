#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# convert_pdfs_to_png.sh
# Recursively find PDFs under a directory and convert to PNG using ImageMagick.
# Supports options for density (DPI), resize, and processing all pages.

usage() {
  cat <<EOF
Usage: $(basename "$0") [options] <root_dir>

Options:
  -d DENSITY    Set ImageMagick density (DPI). Default: 300
  -r RESIZE     Pass a resize string to ImageMagick (e.g. 1024x1024 or 50%)
  -a            Convert all pages of each PDF (creates multiple PNGs: name-pageNNN.png)
  -q QUALITY    Set JPEG/PNG quality (default 90)
  -h            Show this help

Example:
  $(basename "$0") -d 200 -r 1200x -a ./my/folder

Notes:
 - Requires ImageMagick (magick or convert). Prefer 'magick' (ImageMagick 7+).
 - Be careful with very high densities: they create large images and use lots of memory.
EOF
}

# defaults
DENSITY=300
RESIZE=""
ALL_PAGES=0
QUALITY=90
JOBS=0

while getopts ":d:r:ahq:j:" opt; do
  case ${opt} in
    d ) DENSITY=$OPTARG ;;
    r ) RESIZE=$OPTARG ;;
    a ) ALL_PAGES=1 ;;
    j ) JOBS=$OPTARG ;;
    q ) QUALITY=$OPTARG ;;
    h ) usage; exit 0 ;;
    \? ) echo "Invalid Option: -$OPTARG" 1>&2; usage; exit 1 ;;
    : ) echo "Invalid Option: -$OPTARG requires an argument" 1>&2; usage; exit 1 ;;
  esac
done
shift $((OPTIND -1))

if [ $# -lt 1 ]; then
  echo "Error: root directory required.\n"
  usage
  exit 2
fi

ROOT_DIR="$1"

# choose ImageMagick command
if command -v magick >/dev/null 2>&1; then
  IM_CMD=(magick)
elif command -v convert >/dev/null 2>&1; then
  IM_CMD=(convert)
else
  echo "ImageMagick not found (magick or convert). Install ImageMagick and retry." >&2
  exit 3
fi

# Build a temporary NUL-separated list of PDFs to handle weird filenames
PDF_LIST=$(mktemp)
find "$ROOT_DIR" -type f \( -iname '*.pdf' \) -print0 > "$PDF_LIST"

# determine number of jobs: default to number of processors if not supplied
if [ "$JOBS" -le 0 ]; then
  if command -v nproc >/dev/null 2>&1; then
    JOBS=$(nproc)
  else
    JOBS=4
  fi
fi

echo "Found PDFs: $(tr -cd '\0' < "$PDF_LIST" | wc -c)  (jobs=$JOBS)"

# Prefer xargs -0 -P for parallel execution (GNU xargs)
if xargs --help 2>&1 | grep -q "-P"; then
  cat "$PDF_LIST" | xargs -0 -n1 -P "$JOBS" -I{} bash -c '
    pdf="{}"; dir=$(dirname -- "$pdf"); base=$(basename -- "$pdf"); name="${base%.*}";
    if [ "$ALL_PAGES" -eq 1 ]; then
      out_pattern="$dir/${name}-page%03d.png"; first_out=$(printf "$dir/${name}-page%03d.png" 0);
      if [ -f "$first_out" ]; then
        echo "Replacing existing files for: $pdf -> ${name}-pageNNN.png"
        # remove any existing page files so ImageMagick will write fresh ones
        rm -f "$dir/${name}-page"*.png || true
      fi
      echo "Converting all pages: $pdf -> ${name}-pageNNN.png (density=$DENSITY, quality=$QUALITY)"
      "${IM_CMD[@]}" -density "$DENSITY" "$pdf" $( [ -n "$RESIZE" ] && printf '%s ' -resize "$RESIZE" ) -quality "$QUALITY" "$out_pattern"
    else
      out="$dir/${name}.png"
      if [ -f "$out" ]; then
        echo "Replacing existing file: $out"
        rm -f "$out" || true
      fi
      echo "Converting: $pdf -> $out (density=$DENSITY, quality=$QUALITY)"
      "${IM_CMD[@]}" -density "$DENSITY" "$pdf[0]" $( [ -n "$RESIZE" ] && printf '%s ' -resize "$RESIZE" ) -quality "$QUALITY" "$out"
    fi
  '
else
  # Fallback: sequential processing
  echo "xargs -P not available; falling back to sequential processing"
  cat "$PDF_LIST" | while IFS= read -r -d '' pdf; do
    dir=$(dirname -- "$pdf")
    base=$(basename -- "$pdf")
    name="${base%.*}"
    if [ "$ALL_PAGES" -eq 1 ]; then
      out_pattern="$dir/${name}-page%03d.png"
      first_out=$(printf "$dir/${name}-page%03d.png" 0)
      if [ -f "$first_out" ]; then
        echo "Replacing existing files for: $pdf -> ${name}-pageNNN.png"
        rm -f "$dir/${name}-page"*.png || true
      fi
      echo "Converting all pages: $pdf -> ${name}-pageNNN.png (density=$DENSITY, quality=$QUALITY)"
      cmd=("${IM_CMD[@]}" -density "$DENSITY" "$pdf")
      if [ -n "$RESIZE" ]; then
        cmd+=( -resize "$RESIZE" )
      fi
      cmd+=( -quality "$QUALITY" "$out_pattern" )
      "${cmd[@]}"
    else
      out="$dir/${name}.png"
      if [ -f "$out" ]; then
        echo "Replacing existing file: $out"
        rm -f "$out" || true
      fi
      echo "Converting: $pdf -> $out (density=$DENSITY, quality=$QUALITY)"
      cmd=("${IM_CMD[@]}" -density "$DENSITY" "$pdf[0]")
      if [ -n "$RESIZE" ]; then
        cmd+=( -resize "$RESIZE" )
      fi
      cmd+=( -quality "$QUALITY" "$out" )
      "${cmd[@]}"
    fi
  done < "$PDF_LIST"
fi

rm -f "$PDF_LIST"

exit 0
