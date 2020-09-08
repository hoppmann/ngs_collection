columns=$(cat /data/programs/scripts/hoppmann/00-diagnostics/01-filter/01-columns.sh | grep -v "^#" \
| grep -v "^$" | tr "\n" ", " | sed -r 's/(.*),/\1\n/')

echo $columns
