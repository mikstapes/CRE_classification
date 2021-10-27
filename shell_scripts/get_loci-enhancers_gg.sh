#/bin/bash

loci=("tbx5" "hand2" "tbx20")
chroms=("chr15" "chr4" "chr2")
cstart=(11400000 42150000 45800000)
cstop=(12977000 44200000 48000000)

here=$(realpath ..)
enh_set=$here/enhancer_sets/HH22h-enh_set.bed
out=$here/enhancer_sets/loci_with_manual

for i in "${!loci[@]}"; do
        touch $out/"${loci[i]}"_enh.txt
        awk -v chr="${chroms[i]}" -v begin="${cstart[i]}" -v end="${cstop[i]}" '
        {
                if ($1==chr && $2>= begin && $3 <= end) {print $0}
        }' "$enh_set" > $out/"${loci[i]}"_enh_gg.txt
done
