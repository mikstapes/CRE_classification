#! /bin/bash
set -euo pipefail

##Quick heatmap generation with deeptools 
# Plot h27ac, h27me3 around enhancer sets after pks calling, enh calling, and DA analysis

export PATH=/home/phan/miniconda3/envs/chip_atac/bin/:$PATH

here=$(realpath .)
enh_dir=$here/enhancer_sets/IPP
bw_dir=$here/data/bw
outdir=$here/deeptools_matrices
plotdir=$here/plots/deeptools

#score

chick_bws_histone=$(echo "$bw_dir"/{k27ac_22H.bw,k27me3_22H.bw})
chick_bws_atac="$bw_dir"/atac_22H.bw

#regions

mouse_enhs=$(echo $enh_dir/{e105h_ENH_SC.bed,e105h_ENH_PC.bed,e105h_ENH_NC.bed})

#matrix out

mouse_mat_histone=$outdir/e105h_proj_enhancers_histone_5kb.gz 
mouse_mat_atac=$outdir/e105h_proj_enhancers_atac_5kb.gz 


## Compute matrices

computeMatrix reference-point \
  --referencePoint center \
  --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
  --scoreFileName $chick_bws_histone \
  --samplesLabel H3K27ac H3K27me3  \
  --regionsFileName $mouse_enhs \
  --outFileName "$mouse_mat_histone" \
  -p 8

computeMatrix reference-point \
  --referencePoint center \
  --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
  --scoreFileName $chick_bws_atac \
  --samplesLabel ATAC \
  --regionsFileName $mouse_enhs \
  --outFileName "$mouse_mat_atac" \
  -p 8


## Plot heatmaps

Matrix="histone atac"

for M in $Matrix; do

  plotHeatmap \
    --matrixFile $outdir/e105h_proj_enhancers_"$M"_5kb.gz \
    --outFileName $plotdir/e105h_proj_enhancers_"$M"_5kb.pdf \
    --regionsLabel SC PC NC \
    --missingDataColor '#45085c' \
    --colorMap viridis \
    --zMin 0 \
    --refPointLabel ENH \
    --zMax auto \
    --whatToShow 'plot, heatmap and colorbar' \
    --dpi 300 \
    --heatmapHeight 50 \
    --heatmapWidth 12

done





# computeMatrix reference-point \
#   --referencePoint center \
#   --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
#   --scoreFileName $mouse_bws \
#   --samplesLabel E10.5-heart E10.5FL \
#   --regionsFileName $mouse_enhs \
#   --outFileName "$mouse_mat" \
#   -p 8

#mouse_bws_histone=$(echo "$bw_dir"/{k27ac_e105h.bw,k27me3_e105h.bw})
#mouse_bws_atac=$bw_dir/atac_e105h.bw
#chick_enhs=$(echo "$enh_dir"/{HH22_ENH_SC.bed,HH22_ENH_PC.bed,HH22_ENH_NC.bed})
#chick_mat_histone=$outdir/HH22h_proj_enhancers_histone_5kb.gz 
#chick_mat_atac=$outdir/HH22h_proj_enhancers_atac_5kb.gz 
