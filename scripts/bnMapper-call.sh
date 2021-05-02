## /bin/bash

#set -eo pipefail

export PATH=/home/phan/miniconda3/envs/chip_atac/bin:$PATH

here=$(realpath .)
enh_dir=$here/enhancer_sets
out=$here/seq_aln

chain_dir=/project/MDL_Ibrahim/MP_all/references/pairwise-aln
gchain=$chain_dir/mm10.galGal6.all.chain
mchain=$chain_dir/galGal6.mm10.all.chain




# ====== Convert bPk into BED4

#cut -f 1,2,3,4 $peaks_dir/*.broadPeak > $peaks_dir/22H-r1_peaks.bed


# ====== Chain mapping

for e in $(cat ./chick-enhs); do
	bnMapper.py \
	$enh_dir/"$e".bed \
	"$mchain" -i BED -f BED12 --threshold 0.05 -o $out/"$e"_gg6_to_mm10.bed -v info
done

for m in $(cat ./mouse-enhs); do
	bnMapper.py \
	"$enh_dir"/"$m".bed \
	"$gchain" -i BED -f BED12 --threshold 0.05  \
	-o $out/"$m"_mm10_to_gg6.bed -v info
done



		
		


