# deepG <img src="man/figures/logo.png" width="131px" height="140px" align="right" style="padding-left:10px;background-color:white;" />

[![Build Status](https://travis-ci.org/hiddengenome/deepG.svg?branch=master)](https://travis-ci.org/hiddengenome/deepG)
[![codecov](https://codecov.io/gh/hiddengenome/deepG/branch/master/graph/badge.svg)](https://codecov.io/gh/hiddengenome/deepG)

## Overview

deepG is a package for generating CNN/LSTM models from genomic text and provides scripts for various common tasks such as the extraction of cell response. It also comes with example datasets of genomic and human-readable languages for testing.

## Installation and Usage

Please see our [Wiki](https://github.com/hiddengenome/deepG/wiki) for further installation instructions. It covers also usage instructions for multi-GPU machines.

- [Installation on desktop machine](https://github.com/hiddengenome/deepG/wiki/Installation-of-deepG-on-desktop)
- [Installation on GPU server](https://github.com/hiddengenome/deepG/wiki/Installation-of-deepG-on-GPU-server)
- [Installation AWS](https://github.com/hiddengenome/deepG/wiki/Installation-AWS)
- [GPU Usage](https://github.com/hiddengenome/deepG/wiki/manage-GPU-usage)
- [Tensorboard Integration](https://github.com/hiddengenome/deepG/wiki/Tensorboard-integration)

See the help files `?deepG` to get started and for questions use the [FAQ](https://github.com/hiddengenome/deepG/wiki/FAQ).

## Datasets

The library comes with mutiple different datasets for testing:

- The set `data(parenthesis)` contains 100k characters of the parenthesis synthetic language generated from a very simple counting language with a parenthesis and letter alphabet Σ = {( ) 0 1 2 3 4 }. The language is constrained to match parentheses, and nesting is limited to at most 4 levels deep. Each opening parenthesis increases and each closing parenthesis decreases the nesting level, respectively. Numbers are generated randomly, but are constrained to indicate the nesting level at their position.
- The set `data(crispr_full)` containing all CRISPR loci found in NCBI representative genomes with neighbor nucleotides up and downstream.
- The set `data(crispr_sample)` containing a subset of `data(crispr_full)`.
- The set `data(ecoli)` contains the *E. coli* genome, see [the genome sequence of Escherichia coli K-12](https://science.sciencemag.org/content/277/5331/1453.long).
- The set `data(ecoli_small)` contains a subset of `data(ecoli)`.

## Example

### Preprocessing

```r
library(deepG)
data("ecoli") # loads the nucleotide sequence of E. coli
preprocessed <- preprocessSemiRedundant(substr(ecoli, 2, 5000), maxlen = 250) # prepares the batches (one-hot encoding)
```

### Training a language model on CPU

Will generate the binary file `example_full_model.hdf5`. For more options see the Wiki [Training of GenomeNet](https://github.com/hiddengenome/deepG/wiki/Howto-train-GenomeNet).

```r
model <- create_model_lstm_cnn(maxlen = 250, layer_lstm = c(25, 25), layer_dense = c(4))
trainNetwork(model = model, dataset = preprocessed, batch.size = 500, epochs = 5, run.name = "example", tensorboard.log = "log", path.val = "", output = list(none = FALSE, checkpoints =FALSE, tensorboard = FALSE, log = FALSE, serialize_model = FALSE, full_model = TRUE))
```

### Generation of the states
We can now use the trained model to generate neuron responses (states) for a subset of the E coli genome. This will generate a binary file named `states.h5`    

```r
writeStates(model.path = "example_full_model.hdf5", sequence = substr(ecoli, 2, 5000), batch.size = 256, layer.depth = 1, filename = "states", vocabulary = c("a","g","c","t"), step = 1, padding = TRUE)
```

## License and Copyright
Copyright 2019 Philipp Münch

## Supported by

<p float="left">
  <img src="man/figures/hzi.jpg" width="200" />
  <img src="man/figures/dfg.jpg" width="200" />
  <img src="man/figures/bmbf.jpeg" width="100" /> 
  <img src="man/figures/aws.png" width="100" /> 
</p>
