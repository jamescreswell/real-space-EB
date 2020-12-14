# Real Space EB

This repository contains Python functions for computing the E- and B-family decomposition of polarized microwave sky data, as described in [1804.10382](https://arxiv.org/abs/1804.10382) and other papers.

## Idea and motivation

Separation of polarization maps into E/B modes can be realized with a "scalar-type" decomposition, in which E and B are real-valued functions on the sphere, and with a "vector-type" decomposition, in which E and B are each maps of tuples (Q, U) on the sphere.
In principle, the "vector-type" decomposition can be computed in QU space direclty without harmonic transformations via convolution with the appropriate kernel in real space.

## Code

`eb_functions.py` contains helper functions for performing E/B decomposition and generating the real-space E and B kernels.
`eb_functions_examples.ipynb` is a Jupyter containing some examples.

## Prerequisites

Python 3 with `numpy`, `matplotlib`, and `healpy`. Some of the code is based on [`ssht`](https://astro-informatics.github.io/ssht/).
