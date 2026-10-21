# Project notes

## Project setup

### Version control and environment tracking

This project uses git for version control.
To make it reproducible on any machine, it uses a conda environment.
This is where all the dependencies are listed and can be installed from.


### Project blueprint

- immutable raw data
- conda environment
- git tracking
- numbered scripts
- functions folder
- sanity tests folder
- project notes file (this one)
- notebook for exploration and analysis

### General workflow

1. Create and understand a simple script for each step of the project, optionally report here explanations etc.
2. Run it
3. Explore output in the quarto notebook
4. Write down any observation and notes in the notebook
5. Report notes and results here
6. If a problem is found, investigate it further and solve it
7. Commit to git once done

### Sample origin

The 22 samples available in `data/raw/` are taken from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313694

More precisely they are the 22 samples that start with "Spatial transcriptomic profiling of meningioma sample MNG_Sample_0X"

These are the same samples used in https://www.nature.com/articles/s41588-025-02475-w


## Decisions

### Entropy

As of 2026-10-10
Calculate Shannon entropy on raw counts using plugin entropy estimator.


## Analysis notes

### 2026-10-08

Project re-created and set up
Next: start with the entropy function creation and run it on the samples

### 2026-10-16

Todo: run entropy function through a new Seurat script, produce and analyze the H&E plots, run correlation tests; also create a file listing the dataset characteristics (Visium FFPE - formalin fixed; there should probably be more in the SC_Cohort part) like WHO grade (most are grade 2) and methylation class. I could also make one spanorm test for a couple of samples to compare, and also one with the Chao-Shen implementation?

### 2026-10-17

- Refined the entropy function itself, it's mostly unchanged and still works with sparse matrices to be fast.
- Created the first script, `01_load_and_compute_entropy.R`, which creates the Seurat objects, uses the entropy function on the spots and stores the values in a new column.
- It also generates the SpatialFeaturePlots, saving them alongside the H&E images for better visual inspection.
- The script was tested and ran on all the samples, so we now have Seurat objects, plots and a resume csv file for each sample.
- Results folder organization changed, and now each script creates its own subfolder with its outputs for better clarity.
- Thought about the implementation of script 02, but I haven't yet come to a conclusion; the idea is to make it a diagnostic script that tells me some info about correlation and technical confounding.

TODO: finish script 02 and run it on all samples. Based on the results, decide whether SpaNorm or Chao-Shen is the next step (could also have two parallel branches, but I'm still unsure if I should use both tools together or pick only one).

### 2026-10-18

- Completed script 2, evaluated high rho correlation
- Created the new functions to implement Chao-Shen and CWJ
- Created script 3 that applies them to the previously saved Seurat objects
- ran script 3

TODO: analyze script 3 output and create script 4 for correlation analysis

### 2026-10-19

- created script 4 which works much like script 2, and outputs the correlations for both entropy estimators
- the results show that these new estimators show way lower nCounts correlations with respect to plugin entropy; it's still impossible to assess whether it lowers the technical confounding issue or overcorrects the biology itself
- script 4 also shows that the correlation decrease is not consitent across samples, and sometimes is very low while sometimes its so high that it becomes negative, llms suggest a small downsampling experiment on selected samples to investigate this but I should consult with my supervisor
- created script 5 to normalize data using SpaNorm, in a raw version using no chunking and no GPU acceleration because the version in my conda environment (1.4.0) doesn't support it
- created a new branch (gpu-entropy-stemness) to implement a version that uses chunking plus testing cuda with some llm help; this will not be in the final workflow but it's useful for me to run on my limited hardware and test if gpu acceleration gives meaningful speed improvements; I might add this as a bonus in the final repo

TODO: create script 06 for entropy calculation on the normalized data; create script 07 for a comparison across all methods; revise all script functioning and clean up the code

## Project components analysis

### R/shannon_entropy.R

This is the shannon entropy function that is applied to the Seurat object and adds the entropy values to it.
TODO: explain how it works line by line

### scripts/01_load_and_compute_entropy.R

### scripts/02_entropy_depth_correlation.R

### R/chao_shen_entropy.R

### R/chao_wang_jost_entropy.R

### scripts/03_compute_alternative_entropies.Run

### scripts/04_alternative_entropy_depth_correlation.R

### scripts/05_spanorm_normalization

#### 05_spanorm_chunked_cpu

#### 05_spanorm_chunked_gpu
