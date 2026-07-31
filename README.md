# Climbing fiber evoked EPSC input analysis

R workflow for identifying climbing fiber input amplitude groups from
paired-pulse stimulation recordings.


## Contents

- `R/`: analysis script.
- `data-example/`: example input data.
- `results-example/`: example outputs.


## Citation



## Biological context

During the first three postnatal weeks in mice, cerebellar Purkinje neurons are 
initially multi-innervated by several climbing fiber (CF) inputs, which are 
progressively eliminated until reaching mono-innervation (one CF connected per 
Purkinje cell) by ~P21. 

It is therefore essential to ensure an objective estimation of the number of CF 
connected per cell, as well as the relative stregth of these inputs.

In this script, the number of CF inputs is estimated according to the density 
distribution of the amplitudes detected using paired-pulse stimulation. The 
functional development of the inputs is estimated by calculating the disparity 
index and disparity ratio per cell, as previously described in the literature 
(K. Hashimoto and M. Kano, Neuron 2003).


