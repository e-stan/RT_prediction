
from RT_prediction.RT_prediction import *
import matplotlib.pyplot as plt
import numpy as np


#define predictor object, n_jobs sets the number of parralel processes to use
RT_pred = RT_predictor(n_jobs=1)

#get training data from pHILIC. This is molecular descriptors and RTs
training_data = RT_pred.get_training_data_pHILIC_IROA()

to_predict = pd.read_csv("../data/example_smiles.csv",index_col=0)

to_predict = RT_pred.getSmiles(to_predict.index.values)

print(to_predict)

props = RT_pred.computeMolecularDescriptors(to_predict)

props.to_csv("tmp.csv")

print(props)