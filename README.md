# RecursivePCALibUnmixing
Performs library pruning based hyperspectral unmixing using recursive PCA
%Library Pruning by Recursive PCA

This method identifies which spectral library elements are image endmembers. It performs spectral library pruning.

In this work, 
Data: HSI data of size NoPix,NoBand
Lib: Spectral library to be pruned (of size NoLibEle,NoBand)
NoEm: Number of spectral library elements to be pruned

FIrst estimate the number of endmembers (NoEm) by methods such as- HFC-VD [1], GapVD [2] etc. After pruning the library, run SUnSAL algorithm to estiamte the abundance of the pruned library endmembers.

[1] Estimation of number of spectrally distinct signal sources in hyperspectral imagery
[2] Noise robust estimation of number of endmembers in a hyperspectral image by eigenvalue based gap index
