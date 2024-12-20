
# DeepG <img src="man/figures/logo_small.png"  align="left" vspace="-1800px"/>

**deepG: toolbox for deep neural networks optimized for genomic
datasets** <!---
% <p><img alt="DeepG logo" height="70px" src="man/figures/logo_small.png" align="left" hspace="-1000px" vspace="-180px"></p>
-->

The goal of the package is to speed up the development of
bioinformatical tools for sequence classification, homology detection
and other bioinformatical tasks. It is developed for biologists and
advanced AI researchers. DeepG is a collaborative effort from the
McHardy Lab at the *Helmholtz Centre for Infection Research*, the Chair of
Statistical Learning and Data Science at the *Ludwig Maximilian
University of Munich* and the Huttenhower lab at *Harvard T.H. Chan
School of Public Health*.

[![DOI](https://zenodo.org/badge/387820006.svg)](https://zenodo.org/badge/latestdoi/387820006)

## Overview

The package offers several functions to create, train and evaluate
neural networks as well as data processing.

- **Data processing**
  - Create data generator to handle large collections of files.
  - Different options to encode fasta/fastq file (one-hot encoding,
    coverage or quality score encoding).
  - Different options to handle ambiguous nucleotides.
- **Deep learning architectures**
  - Create network architectures with single function call.
  - Custom loss and metric functions available.
- **Model training**
  - Automatically create model/data pipeline.
- **Visualizing training progress**
  - Visualize training progress and metrics in tensorboard.  
- **Model evaluation**
  - Evaluate trained models.
- **Model interpretability**
  - Use Integrated Gradient to visualize relationship of model’s
    predictions with regard to its input.

## Installation

Install the tensorflow python package

``` r
install.packages("tensorflow")
tensorflow::install_tensorflow()
```

and afterwards install the latest version of deepG from github

``` r
devtools::install_github("GenomeNet/deepG")
```

## Usage

See the Package website at <https://deepg.de> for documentation and
example code.

<!-- ## Examples  -->

<!-- ## Datasets -->
<!-- The library comes with mutiple different datasets for testing: -->
<!-- - The set `data(parenthesis)` contains 100k characters of the parenthesis synthetic language generated from a very simple counting language with a parenthesis and letter alphabet Σ = {( ) 0 1 2 3 4 }. The language is constrained to match parentheses, and nesting is limited to at most 4 levels deep. Each opening parenthesis increases and each closing parenthesis decreases the nesting level, respectively. Numbers are generated randomly, but are constrained to indicate the nesting level at their position. -->
<!-- - The set `data(crispr_full)` containing all CRISPR loci found in NCBI representative genomes with neighbor nucleotides up and downstream. -->
<!-- - The set `data(crispr_sample)` containing a subset of `data(crispr_full)`. -->
<!-- - The set `data(ecoli)` contains the *E. coli* genome, see [the genome sequence of Escherichia coli K-12](https://science.sciencemag.org/content/277/5331/1453.long). -->
<!-- - The set `data(ecoli_small)` contains a subset of `data(ecoli)`. -->
<!---
## Installation and Usage
&#10;Please see our [Wiki](https://github.com/hiddengenome/deepG/wiki) for further installation instructions. It covers also usage instructions for multi-GPU machines.
&#10;- [Installation on desktop machine](https://github.com/hiddengenome/deepG/wiki/Installation-of-deepG-on-desktop)
- [Installation on GPU server](https://github.com/hiddengenome/deepG/wiki/Installation-of-deepG-on-GPU-server)
- [Installation AWS](https://github.com/hiddengenome/deepG/wiki/Installation-AWS)
- [GPU Usage](https://github.com/hiddengenome/deepG/wiki/manage-GPU-usage)
- [Tensorboard Integration](https://github.com/hiddengenome/deepG/wiki/Tensorboard-integration)
&#10;See the help files `?deepG` to get started and for questions use the [FAQ](https://github.com/hiddengenome/deepG/wiki/FAQ).
-->
