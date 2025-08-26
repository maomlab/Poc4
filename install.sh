

# General setup instructions for install relevant software

# install ParallelFold to run MSA generation and AlphaFold2 prediction as separate SLURM jobs
# so they can be run on different compute queues
# See: https://github.com/Zuricho/ParallelFold
# Note it is recommended to follow their instructions precisely and create a different conda environment


# install python packages
conda install --yes --file requirements.txt

# install MPLearn for the embed_umap functionality
# see https://github.com/maomlab/MPLearn

# install R packages
Rscript -e "install.packages(c('tidyverse', 'seqinr', 'arrow'))"
