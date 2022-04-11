# Neural-information-processing-and-computations-of-two-input-synapses
This repository contains simulation code for reproducing figures from the following manuscript which is currently under peer review:

Neural information processing and computations of two-input synapses (2022) by Soon Ho Kim, JunHyuk Woo, Kiri Choi, MooYoung Choi, and Kyungreem Han, _Neural Computation, Under Review_.

## Required Installations

The code in ```Fig2BandFigC``` and ```Fig3C``` requires installation of [NEURON](https://neuron.yale.edu/neuron/download) (Python v. 3.7.10 and NEURON v. 7.8.1 was used). Code in ```Fig4``` requires [Tensorflow](https://www.tensorflow.org/install) (v. 2.0.0).

## Figures 2B, 3B and C

Simulations are based on the pyramidal neuron model developed by Wilmes et al. (2018) and uses code from ModelDB Accession number [187603](https://senselab.med.yale.edu/ModelDB/ShowModel?model=187603#tabs-1).

To run simulations, first compile the ```.mod``` files:

```bash
cd Neural-information-processing-and-computations-of-two-input-synapses/Fig2BandFigC
nrnivmodl mod_files
```

To run the simulation for Figure 3, use command

```bash
python3 run_Fig3BC.py 1 0 0.1 0
```

where arguments 1, 0, 0.1, and 0 in the above example are n_1, n_2, α, and the random seed number, respectively.

## Figure 3C

Information for Gidon et al. code

## Figure 4

Information for ANN.
