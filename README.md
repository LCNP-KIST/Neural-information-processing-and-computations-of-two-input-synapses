# Neural-information-processing-and-computations-of-two-input-synapses
This repository contains simulation code for reproducing figures from:

Kim, S. H., Woo, J.H., Choi, K., Choi, M.Y., & Han, K. (2022).
Neural information processing and computations of two-input synapses. _Neural Computation_.

## Required Installations

The code in ```Fig2BandFigC``` and ```Fig3C``` requires installation of [NEURON](https://neuron.yale.edu/neuron/download) (Python v. 3.7.10 and NEURON v. 7.8.1 was used) and [NeuroTools](https://pypi.org/project/NeuroTools/) (v. 0.3.0). Code in ```Fig4``` requires [Tensorflow](https://www.tensorflow.org/install) (v. 2.0.0).

## Figures 2B, 3B and 3C

Simulations use the pyramidal neuron model developed by Wilmes et al. (Wilmes, K. A., Sprekeler, H., & Schreiber, S. (2016). Inhibition as a binary switch for excitatory plasticity in pyramidal neurons. PLoS computational biology, 12(3), e1004768.) and uses code from ModelDB Accession number [187603](https://senselab.med.yale.edu/ModelDB/ShowModel?model=187603#tabs-1). Figure 3 additionally uses code from ModelDB Accession number [185350](https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=185350#tabs-1) to model synapses (Jedlicka, P., Benuskova, L., Abraham, W.C. (2015). A Voltage-Based STDP Rule Combined with Fast BCM-Like Metaplasticity Accounts for LTP and Concurrent "Heterosynaptic" LTD in the Dentate Gyrus In Vivo. PLoS Comput Biol 11, e1004588.

To run simulations, first compile the ```.mod``` files:

```bash
cd Fig2BandFigC
nrnivmodl mod_files
```

To run simulations for Figure 2B, use command

```bash
python3 run_Fig2B.py 1 0
```

where arguments 1and 0 in the above example are _n_ and the random seed number, respectively.
To run the simulation for Figure 3, use command

```bash
python3 run_Fig3BC.py 1 0 0.1 0
```

where arguments 1, 0, 0.1, and 0 in the above example are _n__1, _n__2, _α_, and the random seed number, respectively.

## Figure 2C

Simulations of the pyramidal neuron model developed by Gidon et al. (Gidon, A., Zolnik, T. A., Fidzinski, P., Bolduan, F., Papoutsi, A., Poirazi, P., ... & Larkum, M. E. (2020). Dendritic action potentials and computation in hum) and uses code from ModelDB Accession number [254217](https://senselab.med.yale.edu/ModelDB/ShowModel?model=254217#tabs-1).

To run simulations, first compile the ```.mod``` files:

```bash
cd Fig2C
nrnivmodl _mod
```

```bash
python3 run_Fig2C.py 1 0
```

## Figure 4

This code trains autoencoders based on the architecture of Hinton & Salakhutdinov (Hinton, G. E., & Salakhutdinov, R. R. (2006). Reducing the dimensionality of data with neural networks. science, 313(5786), 504-507.) with a mixture of neurons with two types of activation functions described in the paper. The model is trained on the [Fashion MNIST dataset](https://github.com/zalandoresearch/fashion-mnist). To train and evaluate the model, run the Juptyer notebook ```Fig4/Fig4B.ipynb```.
