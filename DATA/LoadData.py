import pickle
from data_prepare import *

# Load data
my_data = load_my_data(save_to_file="DATA_my.csv")
Denmark_data = load_Denmark_data(save_to_file="DATA_Denmark.csv")

train_combs, test_combs = split_combinations(mode="combination")
test_cat, test_subs, test_catsubs = tuple(test_combs)

cat_data = load_data(Denmark_data, train_combs,
                     test_cat, feature_selection="mutual_info_regression", n_features=0.25)
subs_data = load_data(Denmark_data, train_combs,
                      test_subs, feature_selection="mutual_info_regression", n_features=0.25)
catsubs_data = load_data(Denmark_data, train_combs,
                         test_catsubs, feature_selection="mutual_info_regression", n_features=0.25)

datas = {}
datas["cat_data"] = cat_data
datas["subs_data"] = subs_data
datas["catsubs_data"] = catsubs_data

with open("./CombinationSplit/combination_split.pickle", "wb") as f:
    pickle.dump(datas, f)