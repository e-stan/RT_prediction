# RT_prediction

This package allows for the predicted retention time to be computed from smiles structures and training data.

# Installation

1. Clone git repository: https://github.com/e-stan/RT_prediction.git

2. Install RT_prediction package

```
 > cd <path_to_repo>/RT_prediction/src
 > pip install .
```

# Example Usage

1. Compute molecular properties for compounds of interest:

  For an example of how to compute molecular properties see scripts/compute_molecular_properties.R
  This script predicts molecular properties the compounds in data/example_smiles.csv. This csv file has two columnns. A compound name (or ID) and the smiles string.
  
2. Train model using valid molecular properties 
  
  For an example, see scripts/example_usage_train_predict.ipynb. This notebook uses the included pHILIC training
  data to train the model.
  
3. Predict retention times.

  The example from (2) also uses the properties computed in (1) to predict the retention times and export to a csv file.
  



