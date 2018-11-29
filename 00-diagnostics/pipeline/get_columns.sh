columns=$(cat /data/programs/scripts/hoppmann/00-diagnostics/pipeline/00-columns.sh | grep -v "^#" \
| grep -v "^$" | tr "\n" ", " | sed -r 's/(.*),/\1\n/')
