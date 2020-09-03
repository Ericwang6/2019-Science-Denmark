import numpy as np
import pandas as pd
from scipy.stats import boxcox
import random
from sklearn.decomposition import PCA
from sklearn.feature_selection import mutual_info_regression as MIR
from sklearn.feature_selection import SelectKBest, SelectPercentile

# Load reactions
reaction_data = pd.read_csv("Reaction_data.csv")
reactions = reaction_data["Reaction"]
reactions = reactions.apply(lambda x: x.replace("202_i", "202_vi"))    # correct a mistake in Table S1
ee = reaction_data["% ee"]
catalysts = list(set(reactions.apply(lambda x: x.split("_")[0] + "_" + x.split("_")[1])))
imines = [1, 2, 3, 4, 5]
thiols = ["A", "B", "C", "D", "E"]
products = [f"{i}_{t}" for i in imines for t in thiols]

    
# Split catalyst-imine-thiol combinations
def split_combinations(mode='random', save_to_file=""):
    global catalysts, imines, thiols, products, reaction_data, reactions
    if mode == "random":
        train_combinations = random.sample(list(reactions), 600)
        test_combinations = [c for c in reactions if c not in train_combinations]
    elif mode == "combination":
        train_catalysts = ['242_i', '246_vi', '76_vi', '72_i', '249_i', '205_vi',
                           '365_i', '157_i', '251_vi', '207_i', '253_i', '286_vi',
                           '202_vi', '276_i', '262_vi', '182_i', '181_i', '71_vi',
                           '99_vi', '328_vi', '99_i', '144_i', '230_i', '245_vi']
        test_catalysts = [c for c in catalysts if c not in train_catalysts]
        
        train_subs = [f"{i}_{t}" for i in imines[:-1] for t in thiols[:-1]]
        test_subs = [s for s in [f"{i}_{t}" for i in imines for t in thiols] if s not in train_subs]
        
        train_combinations = [f"{c}_{s}" for c in train_catalysts for s in train_subs]
        test_combinations = []
        # Catalysts Test
        test_combinations.append([f"{c}_{s}" for c in test_catalysts for s in train_subs])
        # Substrate Test
        test_combinations.append([f"{c}_{s}" for c in train_catalysts for s in test_subs])
        # Catalysts/Substrate Test
        test_combinations.append([f"{c}_{s}" for c in test_catalysts for s in test_subs])
        
    elif mode == "ee":
        train_combinations = reaction_data[reaction_data["% ee"] < 80]["Reaction"]
        test_combinations = [c for c in reactions if c not in train_combinations]
    else:
        raise ValueError("Please enter correct ways to split dataset: 'random', 'combination' or 'ee'")
    
    if save_to_file:
        with open(save_to_file, "w+") as f:
            f.write(str(train_combinations) + "\n" + str(test_combinations))
            
    return train_combinations, test_combinations


# Load my data
def load_my_data(save_to_file="", from_file=""):
    if not from_file:
        ASO_catalysts = pd.read_csv("../ASO/Catalysts/ASO_catalysts.csv")
        ASO_imines = pd.read_csv("../ASO/Imines/ASO_imines.csv")
        ASO_thiols = pd.read_csv("../ASO/Thiols/ASO_thiols.csv")
        ASO_products = pd.read_csv("../ASO/Products/ASO_products.csv")

        ESP_catalysts = pd.read_csv("../ESP/Catalysts/ESP_catalysts.csv", index_col=0)
        ESP_imines = pd.read_csv("../ESP/Imines/ESP_imines.csv", index_col=0)
        ESP_thiols = pd.read_csv("../ESP/Thiols/ESP_thiols.csv", index_col=0)

        ESP_thiols.columns = thiols
        ESP_imines.columns = imines
        
        NBO = pd.read_csv("../NBO/NBO.csv", index_col=0)
        
        # Data prepare
        data = {}
        for reaction in reactions:
            x = reaction.split("_")
            data[reaction] = pd.concat([ASO_catalysts[f"{x[0]}_{x[1]}"],
                                        ASO_imines[x[2]],
                                        ASO_thiols[x[3]],
                                        ASO_products[f"{x[2]}_{x[3]}"],
                                        ESP_catalysts[x[0]],
                                        ESP_imines[int(x[2])],
                                        ESP_thiols[x[3]],
                                        NBO[x[3]]])
        data = pd.DataFrame(data)

        # Remove var = 0
        data = data[data.apply(lambda x: x.var() != 0, axis=1)]

        # Scale to unit variance
        data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

        # Calculate ddG
        R = 1.98718
        T = 298.15
        er = (ee + 100) / (100 - ee)
        er.index = reactions
        ddG = R * T * np.log(er) / 1000

        # Box-cox transformation
        shift = abs(ddG.min()) + 1.0
        transformed_ddG, lambda_ = boxcox(ddG + shift)

        # Store ddG to data
        data.loc["ddG"] = ddG
        data.loc["transformed_ddG"] = transformed_ddG
        data.loc["lambda"] = lambda_
        data.loc["shift"] = shift
        
        # Save to .csv file
        if save_to_file:
            data.to_csv(save_to_file)
    # Load data from csv file
    else:
        data = pd.read_csv(from_file, index_col=0, dtype={"Unnamed: 0": np.str})
        
    return data


# Load Denmark data
def load_Denmark_data(save_to_file="", from_file=""):
    if not from_file:
        ASO_catalysts = pd.read_csv("../ASO/BPAASO/combined_ASO_bpas.csv", header=None, index_col=0).T
        ASO_imines = pd.read_csv("../ASO/BPAASO/combined_ASO_imines.csv", header=None, index_col=0).T
        ASO_products = pd.read_csv("../ASO/BPAASO/combined_ASO_products.csv", header=None, index_col=0).T

        ESP_catalysts = pd.read_csv("../ESP/Catalysts/ESP_catalysts.csv", index_col=0)
        ESP_imines = pd.read_csv("../ESP/Imines/ESP_imines.csv", index_col=0)
        ESP_thiols = pd.read_csv("../ESP/Thiols/ESP_thiols.csv", index_col=0)

        ESP_thiols.columns = thiols
        ESP_imines.columns = imines
        
        NBO = pd.read_csv("../NBO/NBO.csv", index_col=0)

        # Data prepare
        data = {}
        for reaction in reactions:
            x = reaction.split("_")
            data[reaction] = pd.concat([ASO_catalysts[f"{x[0]}_{x[1]}"],
                                        ASO_imines[int(x[2])],
                                        ASO_products[f"{x[2]}_{x[3]}"],
                                        ESP_catalysts[x[0]],
                                        ESP_imines[int(x[2])],
                                        ESP_thiols[x[3]],
                                        NBO[x[3]]])
        data = pd.DataFrame(data)
        data = data.dropna(axis=0, how='all')

        # Remove var = 0
        data = data[data.apply(lambda x: x.std() != 0, axis=1)]

        # Scale to unit variance
        data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

        # Calculate ddG
        R = 1.98718
        T = 298.15
        er = (ee + 100) / (100 - ee)
        er.index = reactions
        ddG = R * T * np.log(er) / 1000

        # Box-cox transformation
        shift = abs(ddG.min()) + 1.0
        transformed_ddG, lambda_ = boxcox(ddG + shift)

        # Store ddG to data
        data.loc["ddG"] = ddG
        data.loc["transformed_ddG"] = transformed_ddG
        data.loc["lambda"] = lambda_
        data.loc["shift"] = shift
    
        # Save to .csv file
        if save_to_file:
            data.to_csv(save_to_file)
    # Load from csv file
    else:
        data = pd.read_csv(from_file, index_col=0, dtype={"Unnamed: 0": np.str})
    
    return data


def load_data(data, train_combinations, test_combinations, feature_selection="", n_features=100):
    train = data[train_combinations].T
    test = data[test_combinations].T
    
    train_X = np.array(train.loc[:, :"ddG"])
    train_Y = np.array(train["transformed_ddG"])
    train_ddG = np.array(train["ddG"])

    test_X = np.array(test.loc[:, :"ddG"])
    test_Y = np.array(test["transformed_ddG"])
    test_ddG = np.array(test["ddG"])
    
    if feature_selection == "PCA":
        p = PCA(n_features)
        p.fit(train_X, train_Y)
        train_X = p.transform(train_X)
        test_X = p.transform(test_X)
    if feature_selection == "mutual_info_regression":
        if isinstance(n_features, float):
            selector = SelectPercentile(MIR, percentile=n_features*100)
        elif isinstance(n_features, int):
            selector = SelectKBest(MIR, k=n_features)
        else:
            raise ValueError("n_features must be integer or float")
        selector.fit(train_X, train_Y)
        train_X = selector.transform(train_X)
        test_X = selector.transform(test_X)
        
    
    return {"train_X": train_X, 
            "test_X": test_X,
            "train_Y": train_Y,
            "test_Y": test_Y,
            "train_ddG": train_ddG,
            "test_ddG": test_ddG,
            "lambda": train["lambda"][0],
            "shift": train["shift"][0],
            "train_combs": train_combinations,
            "test_combs": test_combinations}