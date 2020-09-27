#!/Users/apple/opt/anaconda3/envs/deepchem/bin/python
"#!/lustre1/lhlai_pkuhpc/wyz/software/rdkit/bin/python"
from sklearn.linear_model import LassoCV, RidgeCV, ElasticNetCV
from sklearn.cross_decomposition import PLSRegression
from scipy.special import inv_boxcox
import pickle
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import cross_val_score
import time
import pandas as pd
import matplotlib.pyplot as plt
import sys


n = sys.argv[1]
start = time.time()
with open(f"{n}.pickle", "rb") as f:
    datas = pickle.load(f)

models = {}

models["LASSO"] = LassoCV(max_iter=10000, cv=5, n_jobs=-1)
models["RIDGE"] = RidgeCV(cv=5)
models["EN"] = ElasticNetCV(max_iter=10000, cv=5, n_jobs=-1)
models["PLS20"] = PLSRegression(n_components=20, scale=False)

results = {}

for key in models.keys():
    model = models[key]
    results[key] = {"metrics": {}, "data": {}}
    for i in range(len(datas)):
        data = datas[i]
        q2 = cross_val_score(model, data['train_X'], data['train_Y'], cv=5, scoring='r2').mean()
        model.fit(data["train_X"], data["train_Y"])
        predict_Y = model.predict(data["test_X"])
        predict_ddG = inv_boxcox(abs(predict_Y), data["lambda"]) * np.sign(predict_Y) - data["shift"]
        MAD = mean_absolute_error(data['test_ddG'], predict_ddG)
        r2 = r2_score(data['test_ddG'], predict_ddG)
        results[key]["metrics"][f"Trial {i}"] = {"q^2": q2, "MAD": MAD, "R^2": r2}
        results[key]["data"][f"Trial {i}"] = {"predict_ddG": predict_ddG.flatten(), "test_ddG": data['test_ddG']}

writer = pd.ExcelWriter(f"{n}_results.xlsx")
for key in results.keys():
    pd.DataFrame(results[key]["metrics"]).to_excel(writer, sheet_name=key)
writer.save()

with open(f"{n}_results.pickle", "wb") as f:
    pickle.dump(results, f)
# with open(f"{n}_models.pickle", "wb") as f:
#     pickle.dump(models, f)
    
# # Save and plot 
for key in results.keys():
    fig, axes = plt.subplots(5, 2, figsize=(13, 30))
    for i in range(10):
        ax = axes[i // 2][i % 2]
        data = results[key]["data"][f"Trial {i}"]
        R_2 = results[key]["metrics"][f"Trial {i}"]["R^2"]
        MAD = results[key]["metrics"][f"Trial {i}"]["MAD"]
        popt, pcov = curve_fit(lambda x, a, b: a + b * x, data['test_ddG'], data["predict_ddG"])
        regression = "y = {:-.4f}x{:+.4f}".format(popt[1], popt[0])
        ax.scatter(data["test_ddG"], data["predict_ddG"], marker="o", alpha=0.5)
        ax.plot(data["test_ddG"], popt[0] + popt[1] * data["test_ddG"], color="red")
        ax.set_title(f"{key} Trial {1+i}")
        ax.set_xlabel(r"Observed $\Delta \Delta G$")
        ax.set_ylabel(r"Predicted $\Delta \Delta G$")
        ax.text(np.percentile(data["test_ddG"], 90), np.percentile(data["predict_ddG"], 10),
                "${}$\n$MAD={:.4f}$\n$R^2$={:.4f}".format(regression, MAD, R_2))
    plt.savefig(f"{n}_{key}.png", dpi=300)

end = time.time()
print(f"Time: {end-start}s.")
