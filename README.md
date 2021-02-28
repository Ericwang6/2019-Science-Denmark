# 2019-Science-Denmark

```ccheminfolib```是Denmark组开发的软件包

```ASO.py```是计算ASO的代码

43种BPA催化剂的Smiles存储在```Catalysts```文件夹的```Catalysts.csv```中

```Raw_ASO_data.zip```是Denmark组计算的所有化合物的ASO值


**存在的问题**：

+ 和Denmark组计算的ASO有差异


## 2020.08.10 更新

+ 将```ASO.py```和```ASO_plot.py```合并成```ASO.ipynb```文件，计算43种Catalysts，5种Imine，5种Thiol和25种Product的ASO，存储在```ASO_catalysts.csv, ASO_imines.csv, ASO_thiols.csv, ASO_products.csv```中

+ 增加了```ESP.ipynb```文件，计算43种Catalyst，5种Imine，5种Thiol的ESPMAX，ESPMIN值，存储在```ESP_catalysts.csv, ESP_imines.csv, ESP_thiols.csv```中


**说明**：

```.grid```文件中是每个Grid点的ESP值（单位为a.u.）

```.nw```是使用NWChem计算ESP的输入文件

```.xyz```文件是经过DFT优化后的分子坐标，```xxx-xxx.xyz```是优化过程中产生的历史文件


## 2020.08.19 更新

+ 修正了328_vi的Smiles错误

+ 增加了```Automated_SubsExtract.ipynb```文件，用于自动从Catalysts中提取取代基并与NMe4+基团bonding


**说明**：

使用```Chem.ReplaceCore```函数，将催化剂的骨架移除，得到有虚原子（dummy atom）标记取代位点的取代基集合

使用```Chem.GetMolFrags```函数，得到所有取代基


## 2020.08.30 更新

+ 计算了由不同起始构象计算出的ESPMAX的值，用以探究ESPMAX对构象的敏感程度，结果储存在```/ESP/ESP_with_conf```目录下


## 2020.09.03 更新

+ 增加了5种Thiol的NBO分析结果

+ 增加了```data_prepare.py```：定义了一系列用于处理数据的函数

    + ```split_combinations(mode='random', save_to_file="")```：可根据文献中提到的不同split方式对'Catalyst-Substrate'组合进行划分，并将划分结果储存在```save_to_file```中
    
    + ```load_my_data(save_to_file="", from_file="")```和```load_Denmark_data(save_to_file="", from_file="")```：加载数据集并去除方差等于0的维度、进行scale并进行boxcox变换，最后保存在```save_to_file```中。或从```from_file```中加载所需的数据集
    
    + ```load_data(data, train_combinations, test_combinations, feature_selection="", n_features=100)```：按照指定的方式对数据集进行feature selection并划分成training set和test set

+ 增加了```LoadData.py```：按照SI中**Fig S9**的方式对数据集进行split，划分好的数据集以二进制序列文件形式保存在```/CombinationSplit/combination_split.pickle```中

+ 增加了```ModelSelection.py```“：按照SI中**Fig S9**的方式对数据集进行split，使用**Random Forerst, Support Vector Machine(kernel=rbf, poly_2, linear), LassoCV, LassoLarsCV, RidgeCV, ElasticNet, KernelRidge(kernel=rbf, linear)**等模型进行数据处理，训练好的模型以二进制序列保存在```/DATA/CombinationSplit/models.pickle```中，预测值分别以二进制序列文件保存在```/DATA/CombinationSplit/results.pickle```和```/DATA/CombinationSplit/results.pickle```中。

+ 增加了```CaseStudy.ipynb```文件，对Fig 5B和Fig S10两个Case进行了复现：

    + Fig 5B：比较19种Test Catalysts催化的25个Imine-Thiol组合的平均ddG值
    
    + Fig S10：使用mutual_info_regression方法对数据进行降维，取25%的维度，然后比较不同模型在三种Test Set上的MAD和R^2

**问题**：

+ 两个Case Study表明RF的表现异常地好，而SVR的表现很差，与SI中报告的结果不符

## 2020.09.11 更新

+ 增加了```Reaction_data.xlsx```和```Denmark_data_original.xlsx```两份原始数据文件

+ 增加了```SIF```文件夹，其中是复现JACS中提到的**SIF**描述符（只考虑1个构象的ASO描述符）

+ 增加了```ddG_check.ipynb```，用于比较ddG的差异。结果：相差最多0.0015kcal/mol

+ 重写了```data_prepare.py```，在函数中增加了选择描述符的API，以便后续工作中筛选/比较不同描述符的结果

+ 增加了```LoadData.ipynb```，加载不同CaseStudy所需要的数据

+ 在```CaseStudy.ipynb```中增加了Fig 6和JACS-CaseStudy1的复现结果

+ 增加了```RandomSplit```文件夹，其中是根据```Denmark_data_original.xlsx```对数据进行split并测试的代码和结果（文件过大，进行了压缩）

### ！！勘误！！

1. SI中Fig S1将Imine-5的分子式画错了，取代基应是(2',4')二氯苯基

2. ```Reaction_data.xlsx```中所有的**202_i**都应该是**202_vi**

3. ```Denmark_data_original.xlsx```中原本应是**328_vi**的SMILES误写成了**328_i**的SMILES

**结果和问题**

1. RandomSplit的结果表明RF的表现确实过于优秀，但SVR的表现也符合预期，比CombinationSplit的结果好很多。究其原因，可能是SVR有较大的过拟合倾向，后续考虑进一步降低Input的维度。

2. RF在使用SIF+ESPMAX的结果也十分的突出，甚至与ASO+ESPMAX表现无异

**工作计划**

1. 继续进行JACS文章中CaseStudy的复现

2. 寻找新的反应，并在新的反应上测试各种描述符的表现

3. 学习量化软件的使用，复现Sigman的Nature文章

## 2020.11.01 更新
复现工作结束。将MIF-based descriptor计算代码整理成工具包，发布在https://github.com/Ericwang6/MIF-based-Descriptor上，上传了对Zahrt et.al, Science, 2019和Henle et.al, 2020, JACS两篇文章的复现结果```20201024.pptx```
