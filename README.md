# Introduction

CIMPL is the abbreviation for Common Insertion site Mapping PLatfrom. This package identifies common insertion sites (CISs) given insertional mutagenesis data collected across a cohort of tumors. The implementation is based on the Gaussian Kernel Convolution (GKC) framework developed by de Ridder et al. The method was developed for analysis of retroviral insertional mutagenesis data. However, we have extended the method to enable the analysis of transposon insertional mutagenesis as well. In cancer research, well known used retro-viruses are Mouse Mammary Tumor Virus (MMTV) and Murine Leukemia Virus (MuLV). Sleeping Beauty (SB) and piggyBac (PB) are examples of transposon-based systems.

# Getting started

To get started, clone this repository
```
git clone https://github.com/NKI-CCB/cimpl.git
```
and run the R package installer
```
R CMD INSTALL cimpl
```

The package contains a vignette and manual pages with detailed usage information.

# Reference

De Ridder et al. Detecting statistically significant common insertion sites in retroviral insertional mutagenesis screens. PLoS Comput Biol, 2:e166, 2006. [doi:10.1371/journal.pcbi.0020166](http://dx.doi.org/10.1371/journal.pcbi.0020166)
