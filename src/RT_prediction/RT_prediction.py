import pandas as pd
import sklearn.ensemble
import sklearn.tree
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score
import os

class RT_predictor():
    def __init__(self,n_jobs=2):
        self.n_jobs=n_jobs
        pass

    def fit(self,X_train,y_train,numTrees,max_feats,max_depth):
        self.model = sklearn.ensemble.RandomForestRegressor(oob_score=False, n_jobs=self.n_jobs, n_estimators=numTrees, max_features=max_feats,
                                                    max_depth=max_depth)

        self.model.fit(X_train, y_train)

    def predict(self,X_predict):
        return self.model.predict(X_predict)

    def fit_predict(self,X_train,y_train,X_predict,numTree,max_feats,max_depth):
        self.fit(X_train,y_train,numTree,max_feats,max_depth)
        return self.predict(X_predict)

    def get_good_desc(self,X_train):
        cols = X_train.columns.values
        goodCols = []
        for c in cols:
            tmp = X_train.loc[:, c]
            try:
                good = not any(pd.isna(x) for x in tmp.values)
                if good:
                    sdev = np.std([float(x) for x in tmp.values])
                    if sdev > 0:
                        goodCols.append(c)
            except:
                print(c,": non-real valued descriptor")

        return goodCols

    def crossValPrediction(self,X,y,numTree,max_feats,max_depth):
        loo = LeaveOneOut()
        loo.get_n_splits(X)
        y_true = []
        y_preds = []
        for train_index, test_index in loo.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            y_preds.append(self.fit_predict(X_train,y_train,X_test,numTree,max_feats,max_depth))
            y_true.append(y_test)
        return r2_score(y_true, y_preds)

    def get_training_data_pHILIC_IROA(self):
        application_path = os.path.dirname(__file__)
        fn = os.path.join(application_path, "molecular_properties.csv")
        return pd.read_csv(fn)
