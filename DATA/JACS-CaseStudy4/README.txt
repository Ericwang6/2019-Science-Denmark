.png 文件是10次测试的可视化结果
.xlsx 是10次测试的r2, MAD, q2等结果
results.pickle等文件是预测的ddG和原始ddG结果的二进制序列文件，具体请看ModelSelection.py文件
因为输入数据量过大，故没有上传，通过LoadData.ipynb中的JacsCase4模块可以生成所需要的.pickle数据文件

ModelSelection.py是需要命令行执行，带一个参数n，表示要进行预测的数据文件
job.sh是自动化运行所有