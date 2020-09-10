from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV
from scipy.special import inv_boxcox
import pickle
import numpy as np


with open("random_split.pickle", "rb") as f:
    datas = pickle.load(f)

models = {}
models["SVR_POLY2"] = GridSearchCV(SVR(kernel='poly', degree=2),
                                   param_grid={"C": [0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],
                                               "gamma": [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000],
                                               "epsilon": [0.01, 0.1, 0.5, 1, 2, 4]},
                                   cv=5, n_jobs=20)
models["RF"] = GridSearchCV(RFR(n_jobs=-1),
                            param_grid={"n_estimators":[10, 100, 1000, 10000],
                                        "max_features": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]},
                            cv=5, n_jobs=20)
models["LASSO"] = LassoCV(max_iter=100000, cv=5, n_jobs=20)

results = {}

for key in models.keys():
    model = models[key]
    for data in datas:
        model.fit(data["train_X"], data["train_Y"])
        data["predict_Y"] = model.predict(data["test_X"])
        data["predict_ddG"] = inv_boxcox(abs(data["predict_Y"]), data["lambda"]) * np.sign(data["predict_Y"]) - data["shift"]
    results[key] = datas

with open("results.pickle", "wb") as f:
    pickle.dump(results, f)
with open("models.pickle", "wb") as f:
    pickle.dump(models, f)