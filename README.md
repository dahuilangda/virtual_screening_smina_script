# virtual_screening_smina_script

#### 安装说明
#### 1. 安装anaconda, 创建和激活虚拟环境
#### conda create -n env_vs python=3.7
#### conda activate env_vs

#### 2. 安装依赖 pandas, rdkit, openbabel, numpy, tqdm
#### conda install pandas numpy
#### conda install -c openbabel openbabel
#### conda install -c rdkit rdkit
#### pip install tqdm

#### 3. 安装mgltools，默认默认路径（mgltools需要python2.7）
#### 4. 安装smina、qvina（可选）、vina（可选），确保添加到全局环境

#### 5. 运行
####                                    受体蛋白文件   化合物库       化合物库中每个化合物的名称字段               参考化合物文件
#### python vs_smina_V2.0.py --receptor Rec.pdb --db drugbank.sdf --key_value drugbank_id --crystal_ligand crystal_lig.pdb

####                                    受体蛋白文件   化合物库       化合物库中每个化合物的名称字段  口袋坐标
#### python vs_smina_V2.0.py --receptor Rec.pdb --db drugbank.sdf --key_value drugbank_id --center_x -2.44 --center_y -11.425 --center_z 17.436

####                                    受体蛋白文件   化合物库                       化合物库中每个化合物的名称字段                参考化合物文件
#### python vs_smina_V2.0.py --receptor Rec.pdb --db drugbank.sdf -o results.xlsx --key_value drugbank_id --crystal_ligand crystal_lig.pdb --dock smina --parallel 4 --n_cpu 60 --toolkit rdk --size_x 20 --size_y 20 --size_z 20
