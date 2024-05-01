---
editor_options: 
  markdown: 
    wrap: 72
---

**Description :**

This repository contains all code and documentation for CEPHALOPOD -
Comprehensive Ensemble Pipeline for Habitat modelling Across Large-scale
Ocean Plankton Observation Datasets.

It corresponds to the backbone tool for the Blue-Cloud2026 project -
Ecosystem Workbench (WB3). The latter aims to improve the availability,
quality and interoperability of large collections of plankton
observations based on traditional counts, quantitative imaging and
genomic methods available from the OBIS, GBIF, AtlantECO and MATOU data
repositories. It provides a sustainable workflow for production of
ecosystem-level EOVs following clear QA/QC steps and according to best
practices in habitat modelling.

Further development stages will integrate data from EMODnet and ELIXIR
infrastructures. Ultimately, it will serve to the development of a
generic modelling workflow in Blue-Cloud 2026 to generate high-quality
interpolated maps of the global distribution of these plankton entities.

**Installation :**

To run CEPHALOPOD locally, follow the steps below:

\- open an R session on your machine or computing cluster

\- open a terminal at the path where CEPHALOPOD will be installed

\- run the following command:
`git -clone https://github.com/alexschickele/CEPHALOPOD.git`

\- open the file at: `/code/00_config.R`

\- run the R library installation line by line and answer "yes" to any
interactive prompts

\- now it should be ready to use !

**General structure** :

`/master.R` contains an example of script to use CEPHALOPOD. This file
is meant to be modified or duplicated every time you use CEPHALOPOD.

`/code/.` contains the R scripts corresponding to each step of
CEPHALOPOD

`/function/.` contains custom R functions used throughout the code

`/data/.` corresponds to the location of the input data and model cache.
In its current version, some input data are still local to the ETH
Zurich cluster (see further information below)

`/output/.` corresponds to the location where CEPHALOPOD outputs are
automatically stored

**Further information :**

-Technical documentation is provided in the function descriptions. A
detailed documentation will be created in a notebook file.

\- In its current state, the proportion and abundance query are specific
to local data.

\- In its current state, the list of suggested predictors is also
specific to local data. However, one can link to another predictor
collection in the same format (WOA grid, 1x1 degree, depth resolved
optional; variable = LAYER; dimensions = lon, lat, time, depth)
