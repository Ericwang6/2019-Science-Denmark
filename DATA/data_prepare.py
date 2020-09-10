#!/Users/apple/opt/anaconda3/envs/deepchem/bin/python
import numpy as np
import pandas as pd
from scipy.stats import boxcox
import random
from sklearn.decomposition import PCA
from sklearn.feature_selection import mutual_info_regression as MIR
from sklearn.feature_selection import SelectKBest, SelectPercentile, f_regression

# descriptors
catalyst_descriptors = {}
imine_descriptors = {}
thiol_descriptors = {}
product_descriptors = {}

# Load reactions
reaction_data = pd.read_csv("Reaction_data.csv")
reactions = reaction_data["Reaction"]
catalysts = list(set(reactions.apply(lambda x: x.split("_")[0] + "_" + x.split("_")[1])))
imines = ["1", "2", "3", "4", "5"]
thiols = ["A", "B", "C", "D", "E"]
products = [f"{i}_{t}" for i in imines for t in thiols]
ddG = pd.Series(list(reaction_data["ddG"]), index=list(reaction_data["Reaction"]))
    
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


# Load descriptor
def load_descriptor(descriptor, origin="Denmark"):
    global catalyst_descriptors, imine_descriptors, thiol_descriptors, product_descriptors
    # ESP descriptor
    if descriptor == "ESP_catalysts" or descriptor == "ESPMAX_catalysts" or descriptor == "ESPMIN_catalysts":
        raw_ESP_catalysts = pd.read_csv("../ESP/Catalysts/ESP_catalysts.csv", index_col=0)
        ESP_catalysts = pd.DataFrame()
        for col in raw_ESP_catalysts:
            ESP_catalysts[f"{col}_i"] = raw_ESP_catalysts[col]
            ESP_catalysts[f"{col}_vi"] = raw_ESP_catalysts[col]
        if descriptor == "ESPMAX_catalysts":
            catalyst_descriptors[descriptor] = ESP_catalysts.loc["ESPMAX"]
        elif descriptor == "ESPMIN_catalysts":
            catalyst_descriptors[descriptor] = ESP_catalysts.loc["ESPMIN"]
        elif descriptor == "ESP_catalysts":
            catalyst_descriptors[descriptor] = ESP_catalysts
            
    elif descriptor == "ESP_imines" or descriptor == "ESPMAX_imines" or descriptor == "ESPMIN_imines":
        ESP_imines = pd.read_csv("../ESP/Imines/ESP_imines.csv", index_col=0)
        ESP_imines.columns = [str(x) for x in ESP_imines.columns]
        if descriptor == "ESPMAX_imines":
            imine_descriptors[descriptor] = ESP_imines.loc["ESPMAX"]
        elif descriptor == "ESPMIN_imines":
            imine_descriptors[descriptor] = ESP_imines.loc["ESPMIN"]
        elif descriptor == "ESP_imines":
            imine_descriptors[descriptor] = ESP_imines
            
    elif descriptor == "ESP_thiols" or descriptor == "ESPMAX_thiols" or descriptor == "ESPMIN_thiols":
        ESP_thiols = pd.read_csv("../ESP/Thiols/ESP_thiols.csv", index_col=0)
        if descriptor == "ESPMAX_thiols":
            thiol_descriptors[descriptor] = ESP_thiols.loc["ESPMAX"]
        elif descriptor == "ESPMIN_thiols":
            thiol_descriptors[descriptor] = ESP_thiols.loc["ESPMIN"]
        elif descriptor == "ESP_thiols":
            thiol_descriptors[descriptor] = ESP_thiols
    
    elif descriptor == "ESP" or descriptor == "ESPMAX" or descriptor == "ESPMIN":
        for item in ["catalysts", "imines", "thiols"]:
            load_descriptor(f"{descriptor}_{item}")
    
    # ASO descriptor
    elif "ASO" in descriptor and origin == "Denmark":
        if descriptor == 'ASO_catalysts':
            ASO_catalysts = pd.read_csv("../ASO/BPAASO/combined_ASO_bpas.csv", header=None, index_col=0).T
            catalyst_descriptors[descriptor] = ASO_catalysts
        elif descriptor == "ASO_imines":
            ASO_imines = pd.read_csv("../ASO/BPAASO/combined_ASO_imines.csv", header=None, index_col=0).T
            ASO_imines.columns = [str(x) for x in ASO_imines.columns]
            imine_descriptors[descriptor] = ASO_imines
        elif descriptor == "ASO_products":
            ASO_products = pd.read_csv("../ASO/BPAASO/combined_ASO_products.csv", header=None, index_col=0).T
            product_descriptors[descriptor] = ASO_products
        elif descriptor == "ASO":
            for item in ["ASO_catalysts", "ASO_imines", "ASO_products"]:
                load_descriptor(item, origin="Denmark")
        else:
            raise ValueError(f"Descriptor {descriptor} not found.")
    
    elif "ASO" in descriptor and origin == "My":
        if descriptor == 'ASO_catalysts':
            ASO_catalysts = pd.read_csv("../ASO/Catalysts/ASO_catalysts.csv")
            catalyst_descriptors[descriptor] = ASO_catalysts
        elif descriptor == "ASO_imines":
            ASO_imines = pd.read_csv("../ASO/Imines/ASO_imines.csv")
            ASO_imines.columns = [str(x) for x in ASO_imines.columns]
            imine_descriptors[descriptor] = ASO_imines
        elif descriptor == "ASO_thiols":
            ASO_thiols = pd.read_csv("../ASO/Thiols/ASO_thiols.csv")
            thiol_descriptors[descriptor] = ASO_thiols
        elif descriptor == "ASO_products":
            ASO_products = pd.read_csv("../ASO/Products/ASO_products.csv")
            product_descriptors[descriptor] = ASO_products
        elif descriptor == "ASO":
            for item in ["ASO_catalysts", "ASO_imines", "ASO_thiols", "ASO_products"]:
                load_descriptor(item, origin="My")
        else:
            raise ValueError(f"Descriptor {descriptor} not found.")
    
    # NBO descriptor
    elif descriptor == "NBO":
        NBO = pd.read_csv("../NBO/NBO.csv", index_col=0)
        thiol_descriptors[descriptor] = NBO
    
    # SIF descriptor
    elif descriptor == 'SIF_catalysts':
        SIF_catalysts = pd.read_csv("../SIF/Catalysts/SIF_catalysts.csv")
        catalyst_descriptors[descriptor] = SIF_catalysts
    elif descriptor == "SIF_imines":
        SIF_imines = pd.read_csv("../SIF/Imines/SIF_imines.csv")
        SIF_imines.columns = [str(x) for x in SIF_imines.columns]
        imine_descriptors[descriptor] = SIF_imines
    elif descriptor == "SIF_thiols":
        SIF_thiols = pd.read_csv("../SIF/Thiols/SIF_thiols.csv")
        thiol_descriptors[descriptor] = SIF_thiols
    elif descriptor == "SIF_products":
        SIF_products = pd.read_csv("../SIF/Products/SIF_products.csv")
        product_descriptors[descriptor] = SIF_products
    elif descriptor == "SIF":
        for item in ["SIF_catalysts", "SIF_imines", "SIF_thiols", "SIF_products"]:
            load_descriptor(item)
            
    else:
        raise ValueError(f"Descriptor {descriptor} not found.")


# Load data according to given descriptors
def load_raw_data(descriptors=None, save_to_file=None, from_file=None, origin="Denmark"):
    global reactions, ee
    if not from_file:
        if descriptors:
            data = {}
            for des in descriptors:
                load_descriptor(des, origin=origin)
            for reaction in reactions:
                catalyst = reaction.split("_")[0] + "_" + reaction.split("_")[1]
                imine = reaction.split("_")[2]
                thiol = reaction.split("_")[3]
                product = f"{imine}_{thiol}"

                vector = []
                vector += [pd.Series(catalyst_descriptors[key][catalyst]) for key in catalyst_descriptors.keys()]
                vector += [pd.Series(imine_descriptors[key][imine]) for key in imine_descriptors.keys()]
                vector += [pd.Series(thiol_descriptors[key][thiol]) for key in thiol_descriptors.keys()]
                vector += [pd.Series(product_descriptors[key][product]) for key in product_descriptors.keys()]

                data[reaction] = pd.concat(vector)

            data = pd.DataFrame(data)
            data = data.dropna(axis=0, how='all')

            # Remove var = 0
            data = data[data.apply(lambda x: x.var() != 0, axis=1)]

            # Scale to unit variance
            data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

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
            
            return data
        else:
            raise ValueError("descriptors none")
    else:
        data = pd.read_csv(from_file, index_col=0, dtype={"Unnamed: 0": np.str})
        return data


# split data and conduct feature_selection
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
    elif feature_selection == "mutual_info_regression":
        if isinstance(n_features, float):
            selector = SelectPercentile(MIR, percentile=n_features*100)
        elif isinstance(n_features, int):
            selector = SelectKBest(MIR, k=n_features)
        else:
            raise ValueError("n_features must be integer or float")
        selector.fit(train_X, train_Y)
        train_X = selector.transform(train_X)
        test_X = selector.transform(test_X)
    elif feature_selection == "f_regression":
        if isinstance(n_features, float):
            selector = SelectPercentile(f_regression, percentile=n_features*100)
        elif isinstance(n_features, int):
            selector = SelectKBest(f_regression, k=n_features)
        else:
            raise ValueError("n_features must be integer or float")
        selector.fit(train_X, train_Y)
        train_X = selector.transform(train_X)
        test_X = selector.transform(test_X)
    elif feature_selction == "":
        pass
    else:
        raise ValueError(f"Invalid feature selection method: '{feature_selection}'.")
    
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