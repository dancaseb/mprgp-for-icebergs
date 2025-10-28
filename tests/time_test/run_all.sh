ElmerGrid 1 2 square.grd

OUT="time_results.txt"
: > "$OUT" 


echo "=== MPRGP_ELMER ===" >> "$OUT"
{ ./mprgp_elmer.sh |  grep "SOLVER TOTAL TIME"; } &>> "$OUT"

echo "=== ELMER ===" >> "$OUT"
{ ./elmer.sh |  grep "SOLVER TOTAL TIME"; } &>> "$OUT"