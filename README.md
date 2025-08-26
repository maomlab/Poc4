# Deep homology and design of proteasome chaperone proteins in Candida auris
Jackson R. Rapala<sup>1,2</sup>, Mohammad Siddiq<sup>3,4</sup>, Patricia J. Wittkopp<sup>3,4</sup>, Matthew J. O’Meara<sup>2,5*</sup>, Teresa R. O’Meara<sup>1*</sup>
<i><sup>1</sup>Department of Microbiology and Immunology, University of Michigan Medical School;
Ann Arbor MI 48019 USA.
<sup>2</sup>Gilbert S. Omenn Department of Computational Medicine and Bioinformatics, University of
Michigan Medical School; Ann Arbor MI 48019 USA.
<sup>3</sup>Department of Ecology and Evolutionary Biology, University of Michigan; Ann Arbor MI
48019 USA.
<sup>4</sup>Department of Molecular, Cellular, and Developmental Biology, University of Michigan;
Ann Arbor MI 48019 USA.
<sup>5</sup>Department of Medicinal Chemistry, University of Michigan, Ann Arbor, MI 48109, USA
<sup>*</sup>Corresponding author. Email: tromeara@umich.edu, maom@umich.edu</i>

## Abstract
A central tenet of biology is that protein structure mediates the
sequence-function relationship. Recently, there has been excitement
about the promise of advances in protein structure modeling to
generate hypotheses about sequence-structure-function relationships
based on successes with controlled benchmarks. Here, we leverage
structural similarity to identify rapidly evolving proteasome assembly
chaperones and characterize their function in the emerging fungal
pathogen Candida auris. Despite the large sequence divergence, we
demonstrate conservation of structure and function across hundreds of
millions of years of evolution, representing a case of rapid neutral
evolution. Using the functional constraints on structure from these
naturally evolved sequences, we prospectively designed de novo
chaperones and demonstrate that these artificial proteins can rescue
complex biological processes in the context of the whole cell.

## Analysis Scripts
This repository contains the analysis scripts for the study of *C. auris* Poc4 sequence, structure, function:

1. Gather sequence and structure homology
2. Redesign sequences and predict their structure
3. Analyze sequence and structure

To run these scripts, check out this repository

    git clone https://github.com/maomlab/Poc4
    cd Poc4
    # follow the instructions in install.sh to setup the environment

then run the scripts in the `scripts` directory one at a time. This will use data in the `data/` directory, and create files in the `intermediate_data/` and `product/`directories. 
    
The computational workflow builds on the following tools:
  * [Frame2seq](https://github.com/dakpinaroglu/Frame2seq)
  * [ESM](https://github.com/facebookresearch/esm)
  * [AlphaFold](https://alphafoldserver.com/)
  * [FoldSeek](https://search.foldseek.com/search)
  * [ParallelFold](https://github.com/Zuricho/ParallelFold)
  * [MPLearn](https://github.com/maomlab/MPLearn)
  * [BioPython](https://biopython.org/)
  * [Apache Arrow](https://arrow.apache.org/docs/python/index.html)
  * [SLURM](https://slurm.schedmd.com/documentation.html)


### License and Attribution
The computational workflow is released under the MIT License. If you find this workflow useful, please cite:

    Rapala JR, Siddiq M, Wittkopp PJ, O'Meara MJ, O'Meara Deep TR,
    Deep homology and design of proteasome chaperone proteins in Candida auris,
    bioRxiv (2025), DOI: 10.1101/2025.05.14.654010, PMCID: PMC12132334






