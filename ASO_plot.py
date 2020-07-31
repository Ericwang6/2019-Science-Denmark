# 画图
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import concurrent.futures


labels = pd.read_csv("./Catalysts/Catalyst.csv")["label"]
asos = pd.read_csv("./Catalysts/ASO_Catalysts.csv")
raw_asos = pd.read_csv("./BPAASO/combined_ASO_bpas.csv", header=None, index_col=0).stack()

for label in labels:
    print("{} Plotting...".format(label))
    aso = asos[label]
    raw_aso = raw_asos[label]
    fig = plt.figure(figsize=(8,5))
    plt.title(label)
    plt.xlim(-2000, 22000)
    plt.ylim(0, 1.05)
    plt.bar(np.arange(len(aso)), aso, color="C0")
    plt.bar(np.arange(len(raw_aso)), raw_aso, color="C1")
    plt.savefig("./ASO_graph/{}.png".format(label), dpi=300)
    plt.close(fig)
    print("{} Finished".format(label))