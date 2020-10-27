awk -F $'\t' '{ t = $1; $1 = $2; $2 = t; print; }' OFS=$'\t' $1 > $2
