#! /bin/bash
set -euo pipefail

##Quick heatmap generation with deeptools 
# Plot h27ac, h27me3 around enhancer sets after pks calling, enh calling, and DA analysis

export PATH=/home/phan/miniconda3/envs/chip_atac/bin/:$PATH

here=$(realpath .)
enh_dir=$here/enhancer_sets/all_filtered
bw_dir=$here/data/bw
outdir=$here/deeptools_matrices

chick_bws=$(echo $bw_dir/{k27ac_22H.bw,k27ac_22FL.bw})
mouse_bws=$(echo $bw_dir/{k27ac_e105h.bw,k27ac_E105limb.bw})

chick_enhs=$(echo $enh_dir/{HH22-heart_enh.bed,HH24-FL_enh.bed,HH22_heart-FL-shared_enh.bed})
mouse_enhs=$(echo $enh_dir/{e105-heart_enh.bed,e115-FL_enh.bed,e105_heart-FL-shared_enh.bed})

chick_mat=$outdir/chick_enhancers_histones_5kb.gz 
mouse_mat=$outdir/mouse_enhancers_histones_5kb.gz

chick_maps=$here/plots/deeptools/chick_enhancers_histones_5kb.pdf
mouse_maps=$here/plots/deeptools/mouse_enhancers_histones_5kb.pdf



computeMatrix reference-point \
  --referencePoint center \
  --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
  --scoreFileName $chick_bws \
  --samplesLabel HH22-heart HH22-FL  \
  --regionsFileName $chick_enhs \
  --outFileName "$chick_mat" \
  -p 8


computeMatrix reference-point \
  --referencePoint center \
  --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
  --scoreFileName $mouse_bws \
  --samplesLabel E10.5-heart E10.5FL \
  --regionsFileName $mouse_enhs \
  --outFileName "$mouse_mat" \
  -p 8



MATRICES="chick mouse"

for mat in $MATRICES; do
plotHeatmap \
  --matrixFile "$outdir"/"$mat"_enhancers_histones_5kb.gz \
  --outFileName "$here"/plots/"$mat"_enhancers_histones_5kb.pdf \
  --regionsLabel Heart Limb Shared \
  --missingDataColor '#252841' \
  --colorList '#252841,#b1cdca' \
  --zMin 0 \
  --refPointLabel ENH \
  --zMax auto \
  --whatToShow 'plot, heatmap and colorbar' \
  --dpi 300 \
  --heatmapHeight 50 \
  --heatmapWidth 12
done


