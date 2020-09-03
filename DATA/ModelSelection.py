#!/lustre1/lhlai_pkuhpc/wangsw/software/anaconda3/envs/tensorflow/bin/python
import numpy as np
import pickle
from scipy.special import inv_boxcox
from sklearn.ensemble import RandomForestRegressor as RFR
from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LassoCV, LassoLarsCV

# Load Data
with open("combination_split.pickle", "rb") as f:
    datas = pickle.load(f)

models = {"RF": RFR(n_estimators=100)}
MADs = {}
results = {}
lambda_ = datas['cat_data']['lambda']
shift = datas['cat_data']['shift']

models["RF"] = GridSearchCV(RFR(n_jobs=-1),
                            param_grid={"n_estimators":[10, 100, 1000, 10000],
                                        "max_features": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]},
                            cv=5, n_jobs=20)
models["LASSOLARS"] = LassoLarsCV(max_iter=5000, cv=5, n_jobs=-1)
models["SVR_POLY2"] = GridSearchCV(SVR(kernel='poly', degree=2),
                                   param_grid={"C": [0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],
                                               "gamma": [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000],
                                               "epsilon": [0.01, 0.1, 0.5, 1, 2, 4]},
                                   cv=5, n_jobs=20)

for key in models.keys():
    model = models[key]
    model.fit(datas['cat_data']['train_X'], datas['cat_data']['train_Y'])

    train_predict_Y = model.predict(datas["cat_data"]['train_X'])
    train_predict_ddG = inv_boxcox(abs(train_predict_Y), lambda_) * np.sign(train_predict_Y) - shift
    train_MAD = abs(datas['cat_data']['train_ddG'] - train_predict_ddG).mean()

    cat_predict_Y = model.predict(datas['cat_data']['test_X'])
    cat_predict_ddG = inv_boxcox(abs(cat_predict_Y), lambda_) * np.sign(cat_predict_Y) - shift
    cat_MAD = abs(datas['cat_data']['test_ddG'] - cat_predict_ddG).mean()

    subs_predict_Y = model.predict(datas['subs_data']['test_X'])
    subs_predict_ddG = inv_boxcox(abs(subs_predict_Y), lambda_) * np.sign(subs_predict_Y) - shift
    subs_MAD = abs(datas['subs_data']['test_ddG'] - subs_predict_ddG).mean()

    catsubs_predict_Y = model.predict(datas['catsubs_data']['test_X'])
    catsubs_predict_ddG = inv_boxcox(abs(catsubs_predict_Y), lambda_) * np.sign(catsubs_predict_Y) - shift
    catsubs_MAD = abs(datas['catsubs_data']['test_ddG'] - catsubs_predict_ddG).mean()

    MADs[key] = {"train": train_MAD, "cat": cat_MAD, "subs": subs_MAD, "catsubs_MAD": catsubs_MAD}
    results[key] = {"train_predict_Y": train_predict_Y, "train_predict_ddG": train_predict_ddG,
                    "cat_predict_Y": cat_predict_Y, "cat_predict_ddG": cat_predict_ddG,
                    "subs_predict_Y": subs_predict_Y, "subs_predict_ddG": subs_predict_ddG,
                    "catsubs_predict_Y": catsubs_predict_Y, "catsubs_predict_ddG": catsubs_predict_ddG}

# Save Models
with open("models.pickle", "wb") as f:
    pickle.dump(models, f)
# Save MADs
with open("MADs.txt", "w+") as f:
    f.write(str(MADs))
# Save Results
with open("results.pickle", "wb") as f:
    pickle.dump(results, f)