#/bin/bash

loci=("tbx5" "hand2" "tbx20")
chroms=("chr5" "chr8" "chr9")
cstart=(118200000 55500000 24000000)
cstop=(121000000 59000000 26000000)

here=$(realpath ..)
enh_set=$here/enhancer_sets/e105h-enh_set.bed
out=$here/enhancer_sets/loci_with_manual

for i in "${!loci[@]}"; do
        touch $out/"${loci[i]}"_enh.txt
        awk -v chr="${chroms[i]}" -v begin="${cstart[i]}" -v end="${cstop[i]}" '
        {
                if ($1==chr && $2>= begin && $3 <= end) {print $0}
        }' "$enh_set" > $out/"${loci[i]}"_enh.txt
done
