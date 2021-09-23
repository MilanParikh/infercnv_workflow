# infercnv_terra
Terra workflow for inferCNV with cell caching

## Automation Script Jupyter Notebook
This notebook/script helps you set up the input files and bash scripts to run the R version of inferCNV on Terra from a Python (ScanPy) environment

Instructions contained within the comments show where to configure paths, Google buckets, Terra information, your email, etc.

There are two options:
1. Per Sample (*Recommended*) - This will split up each of your samples into it's own parallel Terra run, which takes less time, requires less resources, and costs much less
2. All Samples - This runs your whole count matrix at one time, this can take days, requires very high resources for larger matrices, and can cost quite a bit.

## Docker
Dockerfile - contains the R and Python setup required for version 1.8.1 of inferCNV, with a modified inferCNV.R script that enables the up_to_step flag

(Already built here: https://hub.docker.com/repository/docker/mparikhbroad/infercnv_caching and referenced in the WDL)

## WDL
The workflow is set up with 6 roughly evenly timed stages using the "up_to_step" flag in order to enable cell caching. This allows for a long running inferCNV process to not lose too much progress if preempted.

(Already uploaded here: https://portal.firecloud.org/?return=terra#methods/mparikh/infercnv_caching/3)

At the end of each step the outputs are compressed into a .tar.gz archive and then the follow step will untar and resume
