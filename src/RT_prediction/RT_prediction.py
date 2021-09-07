import pandas as pd
import sklearn.ensemble
import sklearn.tree
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score
import os
import requests
import uuid
import sys
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import norm

class RT_predictor():
    def __init__(self,n_jobs=2):
        self.n_jobs=n_jobs
        pass

    def fit(self,X_train,y_train,numTrees,max_feats,max_depth,weights = "none"):
        if type(weights) == type(""):
            weights = [1 for _ in X_train]
        self.model = sklearn.ensemble.RandomForestRegressor(oob_score=False, n_jobs=self.n_jobs, n_estimators=numTrees, max_features=max_feats,
                                                    max_depth=max_depth)

        self.model.fit(X_train, y_train,sample_weight=weights)

    def predict(self,X_predict):
        return self.model.predict(X_predict)

    def fitPdf(self,xs):
        mean,std = norm.fit(xs)

        #prev = 0
        #for t in ts:
        #    tmp = len([x for x in xs if x <= t])
        #    pdf.append(tmp-prev)
        #    prev = tmp
        return mean,std



    def predict_pdf(self,X):
        fullPred = self.model.predict(X)
        pdfs = []

        for x in range(len(X)):
            tmp = np.array([X[x]])
            preds = []
            for pred in self.model.estimators_:
                preds.append(pred.predict(tmp)[0])
            pdf = self.fitPdf(preds)
            pdfs.append(pdf)
        return pdfs

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

    def computeMolecularDescriptors(self,queries):
        id = str(uuid.uuid1())
        queries.to_csv(id)
        outfn = id + '_computed'
        junk = id + "_jnk"
        junk2 = id + "_jnk2"

        dir = os.path.dirname(__file__)

        command = "Rscript --vanilla " + os.path.join(dir,"compute_molecular_properties.R") + " " + id + " " + outfn + " > " + junk + "2>" + junk2

        os.system(command)

        props = pd.read_csv(outfn)

        if os.path.isfile(id): os.remove(id)
        if os.path.isfile(outfn): os.remove(outfn)
        if os.path.isfile(junk): os.remove(junk)
        if os.path.isfile(junk2): os.remove(junk2)

        return props

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

    def getSmiles(self,inchikeys):
        outdict = {}
        for ick in inchikeys:

            smiles = InChiKeyToSMILES(ick.replace(" (M+1)",""))
            if len(smiles) > 0:
                outdict[ick] = {"SMILES":smiles}
        return pd.DataFrame.from_dict(outdict,orient="index")


def InChiKeyToInChi(inchikey):
    host = "http://www.chemspider.com"
    getstring = "/InChI.asmx/InChIKeyToInChI?inchi_key="

    try:
        r = requests.get('{}{}{}'.format(host, getstring, inchikey))
        if r.ok:
            res = str(
                r.text.replace('<?xml version="1.0" encoding="utf-8"?>\r\n<string xmlns="http://www.chemspider.com/">',
                               '').replace('</string>', '').strip())
        else:
            res = ""
        return res
    except:
        return ""


def InChiToSMILES(inchi):
    host = "http://www.chemspider.com"
    getstring = "/InChI.asmx/InChIToSMILES?inchi="

    try:
        r = requests.get('{}{}{}'.format(host, getstring, inchi))
        if r.ok:
            res = str(
                r.text.replace('<?xml version="1.0" encoding="utf-8"?>\r\n<string xmlns="http://www.chemspider.com/">',
                               '').replace('</string>', '').strip())
        else:
            res = ""
        return res
    except:
        return ""

def InChiKeyToSMILES(inchikey):
    inchi = InChiKeyToInChi(inchikey)
    if len(inchi) > 0:
        return InChiToSMILES(inchi)
    else:
        return ""