About Tox21 Data
================


## Data Source
Tox21 data can be downloaded here:  
https://tripod.nih.gov/tox21/assays/


## Data Processing
The SMILES in the downloaded tox21 data are problematic, so we need to fetch
the correct SMILES via PubChem. This also indicates that compounds without
PubChemID will not be considered.

`Tox21_compounds.ipynb` shows how to process the tox21 data and convert the
PubChemID into `std_inchi_key`


## Overlap with Janssen's data
We take mitochondrial toxicity as an example to show the overlap. Directly using
PubChemID to match the overlap may underestimate the number of overlapped
compunds because some compounds are actually the same if we standardise them.
