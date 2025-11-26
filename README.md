This repository contains Supplementary Dataset for the manuscript:

“The dermal skeleton of stem-actinopterygian Moythomasia durgaringa and its implications for the nature of the ancestral osteichthyan”
by Xianren Shan et al. (submitted to Journal of Morphology)

These materials include all datasets and scripts required to reproduce the ancestral state estimation (ASE) analyses performed in the study, using Structured Markov Models (SMMs).

Tree files (*.tree). provide alternative phylogenetic hypotheses used for structured Markov model (SMM) ancestral state reconstruction.
• tree.tree is based on Lu et al. (2016).
• tree2.tree and tree3.tree follow the topologies from King et al. (2017).
Character matrices (*.nex). include discrete histological characters coded from dermal skeletal elements of the cranial and trunk regions.
• cranial.nex and trunk.nex corresponds to tree.tree.
• cranial2.nex and trunk2.nex corresponds to tree2.tree and tree3.tree.
Tip-age files (*.csv). contain fossil age information formatted for tip-dating procedures.
•tip.age.csv corresponds to tree.tree.
•tip.age2.csv corresponds to tree2.tree and tree3.tree.
R script (ase.R). performs data import, model specification, and ancestral state estimation using the structured Markov model framework. It also includes commands for loading the phylogenies, processing the character matrices, and visualizing reconstructed states.

If you use these datasets, please cite the manuscript above once published.
