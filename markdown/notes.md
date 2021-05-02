# Enhancers analysis notes

## Enhancer calling & distal enhancer filtering

See get_enhancers.Rmd

## Enhancer-target genes

1. Get insulating domains/TAD coords, either by tool for genome-wide or manual for specific loci
- https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html
- 
2. Get expression counts (salmon?) of all tx within domains
3. Get contact counts between enhancers and gene/tx TSS/proms within domains
4. Correlate exp with contacts --> set of active enhancers/gene-tx

# R notes

Base dir in a Rmd is the dir the .Rmd script lives in. 
Base dir in a .R script is the top level of the project file. 
Use here to clear confusions!

See: https://jenrichmond.rbind.io/post/how-to-use-the-here-package/

