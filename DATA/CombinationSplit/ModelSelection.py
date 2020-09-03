import numpy as np
import pandas as pd
import pickle
from scipy.special import inv_boxcox
from sklearn.ensemble import RandomForestRegressor as RFR
from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LassoCV, LassoLarsCV, RidgeCV, ElasticNetCV
from sklearn.kernel_ridge import KernelRidge

# Load Data
with open("combination_split.pickle", "rb") as f:
    datas = pickle.load(f)

results = {}

# boxcox-shift params
lambda_ = datas['cat_data']['lambda']
shift = datas['cat_data']['shift']

# models
models = {}
models["RF"] = GridSearchCV(RFR(n_jobs=-1),
                            param_grid={"n_estimators":[10, 100, 1000, 10000],
                                        "max_features": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]},
                            cv=5, n_jobs=20)
models["LASSO"] = LassoCV(max_iter=100000, cv=5, n_jobs=20)
models["RIDGE"] = RidgeCV(cv=5)
models["LASSOLARS"] = LassoLarsCV(max_iter=5000, cv=5, n_jobs=-1)
models["SVR_POLY2"] = GridSearchCV(SVR(kernel='poly', degree=2),
                                   param_grid={"C": [0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],
                                               "gamma": [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000],
                                               "epsilon": [0.01, 0.1, 0.5, 1, 2, 4]},
                                   cv=5, n_jobs=20)
models["SVR_RBF"] = GridSearchCV(SVR(kernel='rbf'),
                                 param_grid={"C": [0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],
                                             "gamma": [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000],
                                             "epsilon": [0.01, 0.1, 0.5, 1, 2, 4]},
                                 cv=5, n_jobs=20)
models["SVR_LINEAR"] = GridSearchCV(SVR(kernel='linear'),
                                    param_grid={"C": [0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],
                                                "epsilon": [0.01, 0.1, 0.5, 1, 2, 4]},
                                    cv=5, n_jobs=20)
models["EN"] = ElasticNetCV(max_iter=100000, cv=5, n_jobs=20)

models['KR_RBF'] = GridSearchCV(KernelRidge(kernel='rbf'),
                                param_grid={"alpha": [10, 1, 0.1, 1e-2, 1e-3],
                                            "gamma": [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000]},
                                cv=5, n_jobs=20)
models['KR_LINEAR'] = GridSearchCV(KernelRidge(kernel='linear'),
                                   param_grid={"alpha": [10, 1, 0.1, 1e-2, 1e-3]},
                                   cv=5, n_jobs=20)

for key in models.keys():
    model = models[key]
    model.fit(datas['cat_data']['train_X'], datas['cat_data']['train_Y'])

    train_predict_Y = model.predict(datas["cat_data"]['train_X'])
    train_predict_ddG = inv_boxcox(abs(train_predict_Y), lambda_) * np.sign(train_predict_Y) - shift
    train_observed_ddG = datas['cat_data']['train_ddG']
    train_combs = datas['cat_data']['train_combs']

    cat_predict_Y = model.predict(datas['cat_data']['test_X'])
    cat_predict_ddG = inv_boxcox(abs(cat_predict_Y), lambda_) * np.sign(cat_predict_Y) - shift
    cat_observed_ddG = datas['cat_data']['test_ddG']
    cat_test_combs = datas['cat_data']['test_combs']

    subs_predict_Y = model.predict(datas['subs_data']['test_X'])
    subs_predict_ddG = inv_boxcox(abs(subs_predict_Y), lambda_) * np.sign(subs_predict_Y) - shift
    subs_observed_ddG = datas['subs_data']['test_ddG']
    subs_test_combs = datas['subs_data']['test_combs']

    catsubs_predict_Y = model.predict(datas['catsubs_data']['test_X'])
    catsubs_predict_ddG = inv_boxcox(abs(catsubs_predict_Y), lambda_) * np.sign(catsubs_predict_Y) - shift
    catsubs_observed_ddG = datas['catsubs_data']['test_ddG']
    catsubs_test_combs = datas['catsubs_data']['test_combs']

    results[key] = {"train_observed_ddG": train_observed_ddG, "train_predict_ddG": train_predict_ddG, "train_combs": train_combs,
                    "cat_observed_ddG": cat_observed_ddG, "cat_predict_ddG": cat_predict_ddG, "cat_test_combs": cat_test_combs,
                    "subs_observed_ddG": subs_observed_ddG, "subs_predict_ddG": subs_predict_ddG, "subs_test_combs": subs_cat_combs,
                    "catsubs_observed_ddG": catsubs_observed_ddG, "catsubs_predict_ddG": catsubs_predict_ddG, "catsubs_test_combs": catsubs_test_combs}

# Save Models
with open("models.pickle", "wb") as f:
    pickle.dump(models, f)
    
# Save Results
with open("results.pickle", "wb") as f:
    pickle.dump(results, f)

# Save to excel
writer = pd.ExcelWriter("results.xlsx")
excel = {}
excel['train'] = pd.DataFrame(index=datas['cat_data']['train_combs'])
excel['test_cat'] = pd.DataFrame(index=datas['cat_data']['test_combs'])
excel['test_subs'] = pd.DataFrame(index=datas['subs_data']['test_combs'])
excel['test_catsubs'] = pd.DataFrame(index=datas['catsubs_data']['test_combs'])
for model in results:
    excel['train'][f"{model}_Observed"] = results[model]["train_observed_ddG"]
    excel['train'][f"{model}_Predict"] = results[model]["train_predict_ddG"]
    
    excel['test_cat'][f"{model}_Observed"] = results[model]["cat_observed_ddG"]
    excel['test_cat'][f"{model}_Predict"] = results[model]["cat_predict_ddG"]
    
    excel['test_subs'][f"{model}_Observed"] = results[model]["subs_observed_ddG"]
    excel['test_subs'][f"{model}_Predict"] = results[model]["subs_predict_ddG"]
    
    excel['test_catsubs'][f"{model}_Observed"] = results[model]["catsubs_observed_ddG"]
    excel['test_catsubs'][f"{model}_Predict"] = results[model]["catsubs_predict_ddG"]
for key in excel.keys():
    excel[key].to_excel(writer, sheet_name=key)
writer.save()