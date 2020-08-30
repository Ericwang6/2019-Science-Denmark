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