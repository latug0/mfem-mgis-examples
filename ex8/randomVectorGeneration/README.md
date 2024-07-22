# Random Orthogonal Vector Pair Generator

This project provides a simple Python script to generate a specified number of pairs of orthogonal vectors. In mfem-mgis-examples, these vectors are used to describe the principal directions of crystals. The script is used to randomly generate the principal directions of crystals within a polycrystal.

## Prerequisites

To run this script, you need to have Python installed on your system. The script requires the `Numpy` Python package

## Usage

The script generates pairs of orthonormal vectors. You can specify the number of vectors to be generated and the algorithm to use.

### Arguments

- nbVec (int): The number of vectors to generate.
- algoKey (str): The algorithm to use for generating vectors. It can be either 'default' or 'angles'. 'default' uses a standard method for generating orthonormal vectors, while 'angles' generates orthonormal vectors using specific angles.

### Algorithms

#### Default

The “default” algorithm generates random vectors by generating vector coordinates according to a normal distribution.

#### Angles

The “angles” algorithm generates random vector directions using Euler angles obtained from a uniform distribution. This method is based on the approach described in Appendix A of the paper *Approche micromécanique du comportement du combustible dioxyde d'uranium* by Julian Soulacroix (2014).