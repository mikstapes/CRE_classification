#! /bin/bash
set -euo pipefail

##Quick heatmap generation with deeptools 
# Plot h27ac, h27me3 around enhancer sets after pks calling, enh calling, and DA analysis

export PATH=/home/phan/miniconda3/envs/chip_atac/bin/:$PATH

here=$(realpath ..)
enh_dir=$here/enhancer_sets
bw_dir=$here/data/bw
outdir=$here/deeptools_matrices

chick_bws="$bw_dir"/{k27ac_22H.bw,k27ac_22FL.bw,k27me3_22H.bw,k27me3_22FL.bw}
mouse_bws="$bw_dir"/{k27ac_e105h.bw,k27ac_E115FL.bw,k27me3_e105h.bw,k27me3_e115FL.bw}

chick_enhs="$enh_dir"/{HH22-enh_heart,HH24-enh_FL,HH22_shared-enh}
mouse_enhs="$enh_dir"/{E105-enh_heart,E115-enh_FL,E105_shared-enh}

chick_mat=$outdir/chick_enhancers_histones_5kb.gz 
mouse_mat=$outdir/mouse_enhancers_histones_5kb.gz

chick_maps=$here/plots/deeptools/chick_enhancers_histones_5kb.pdf
mouse_maps=$here/plots/deeptools/mouse_enhancers_histones_5kb.pdf



computeMatrix reference-point \
  --referencePoint center \
  --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
  --scoreFileName $chick_bws \
  --samplesLabel H3K27ac_Heart H3K27ac_FL H3K27me3_Heart H3K27me3_FL \
  --regionsFileName $chick_enhs \
  --outFileName "$chick_mat" \
  -p 8


computeMatrix reference-point \
  --referencePoint center \
  --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
  --scoreFileName $mouse_bws \
  --samplesLabel H3K27ac_Heart H3K27ac_FL H3K27me3_Heart H3K27me3_FL \
  --regionsFileName $mouse_enhs \
  --outFileName "$mouse_mat" \
  -p 8



MATRICES="$chick_mat $mouse_mat"

for mat in $MATRICES; do
  
plotHeatmap \
  --matrixFile "$mat" \
  --outFileName
  --regionsLabel Heart Limb Shared \
  --missingDataColor '#383d61' \
  --colorMap '#383d61,#9CBFBB' \
  --zMin 0 \
  --refPointLabel Enh \
  --whatToShow 'heatmap and colorbar' \
  --zMax auto \
  --xAxisLabel '' \
  --dpi 100 \
  --heatmapHeight 15 \
  --heatmapWidth 3
done

if [ $? -ne 0 ]
then
  echo "!! plotting failed"
  exit 1
else
  echo "Plotting done"
fi

    


  
    
