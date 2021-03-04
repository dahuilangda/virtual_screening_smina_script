from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import pandas as pd
import os
import multiprocessing as mp
from tqdm import tqdm
import numpy as np

from argparse import ArgumentParser

import warnings
warnings.filterwarnings('ignore')

def prepareReceptor(filename, output_path='./tmp'):
    if not os.path.exists(filename):
        raise Exception("受体文件不存在!")

    if filename.split('.')[-1] != 'pdb':
        raise Exception("受体请使用pdb格式!")

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    output_file_tmp = filename.split('/')[-1].split('.')[0] + '.pdbqt'
    rec = os.path.join(output_path, output_file_tmp)
    cmd = '/usr/bin/python2.7 /usr/lib/python2.7/dist-packages/AutoDockTools/Utilities24/prepare_receptor4.py -r %s -A hydrogens -o %s > /dev/null 2>&1' %(filename, rec)
    print('Step 1')
    print('###############################')
    print('正在准备受体文件!!!')
    os.system(cmd)
    
    print('受体文件准备完毕!!!')
    print('###############################')
    return rec

def getPocketsFromCrystalLigand(filename):
    # 计算共结晶配体的中心点
    if not os.path.exists(filename):
        raise Exception("配体文件不存在!")

    if filename.split('.')[-1] != 'pdb':
        raise Exception("共结晶配体请使用pdb格式!")

    with open(filename, 'r') as f:
        lines = f.readlines()
    l = []
    for line in lines:
        if ('HETATM' in line) or ('ATOM' in line):
            a = [float(x) for x in line.split()[5:8]]
            l.append(a)
    l = np.array(l)
    center_x = round(np.average(l[:,0]), 3)
    center_y = round(np.average(l[:,1]), 3)
    center_z = round(np.average(l[:,2]), 3)
    print('')
    print('Step 2')
    print('###############################')
    print('选择%s,%s,%s作为口袋中心点!!!'%(center_x, center_y, center_z))
    print('###############################')
    return center_x, center_y, center_z

def loadDrugDatabase(filename):
    '''
    输入文件支持sdf和csv格式
    如果是csv格式必须存存在smiles列,例如 smiles, id, name, ....
                                    CCC(CO)NC, 1, HD0001, ...
    '''
    if not os.path.exists(filename):
        raise Exception("数据库文件不存在!")

    if not ((filename.split('.')[-1] == 'csv') or (filename.split('.')[-1] == 'sdf')):
        raise Exception("数据库文件仅支持sdf和csv格式!")

    print('')
    print('Step 3')
    print('###############################')
    print('正在加在数据库!!!')

    if filename.split('.')[-1] == 'sdf':
        drugs = Chem.SDMolSupplier(filename)
        # drugs = [x for x in drugs]
        descr = []
        smiles = []
        for x in drugs:
            descr.append(x.GetPropsAsDict())
            smiles.append(Chem.MolToSmiles(x))
        col = [i for i in x.GetPropsAsDict().keys()]
        df = pd.DataFrame(columns=[col])
        df = df.from_dict(descr)

        df['smiles_tmp'] = smiles
        print('正在去除化合物中的盐离子!!!')
        df['smiles'] = [removeSalt(smi) for smi in df.smiles_tmp]
        print('正在产生结构文件!!!')
        PandasTools.AddMoleculeColumnToFrame(df, 'smiles', 'mol')

    elif filename.split('.')[-1] == 'csv':
        df = pd.read_csv(filename)
        df['smiles_tmp'] = df['smiles']

        print('正在去除化合物中的盐离子!!!')
        df['smiles'] = [removeSalt(smi) for smi in df.smiles_tmp]
        print('正在产生结构文件!!!')
        PandasTools.AddMoleculeColumnToFrame(df, 'smiles', 'mol')
        # drugs = [x for x in df.mol]

    else:
        raise Exception('化合物数据库必须是SDF或CSV格式！')
    
    print('共加载%s个化合物!!!' %df.shape[0])
    print('###############################')
    return df

def removeSalt(smi):
    # 去除salt
    smi_list = smi.split('.')
    max_list = []
    for step, l in enumerate(smi_list):
        max_list.append(len(l))
    ndx = max_list.index(max(max_list))
    return smi_list[ndx]

def smi2sdf(smi, output_name):
    try:
        smi = removeSalt(smi)
        mol = Chem.MolFromSmiles(smi)
        hmol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(hmol,AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(hmol,1000)
        writer = Chem.SDWriter(output_name)
        writer.write(hmol)
        writer.close()
        return 1
    except Exception as e:
        return 0

def spilt2SingleMolecule(df, input_path, output_path, toolkit='rdk', key='name'):
    if not os.path.exists(input_path):
        os.makedirs(input_path)

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if not key in df.columns:
        raise Exception('关键字输入错误！')

    filenames = []
    outputs = []
    print('')
    print('Step 4')
    print('###############################')
    print('正在将分割化合物到%s文件夹，请耐心等待!!!' %input_path)

    for step, smi in enumerate(tqdm(df.smiles)):
        try:
            filename = os.path.join(input_path, df[key].iloc[step]+'.sdf')
            output_name = os.path.join(output_path, df[key].iloc[step] + '_docked.sdf')
            filenames.append(filename)
            outputs.append(output_name)
            if not os.path.exists(filename):
                if len(smi) <= 100: # 少于100个原子
                    if toolkit == 'ob':
                        smi = removeSalt(smi)
                        cmd = 'echo \'%s\' > tmp.smi' %smi
                        os.system(cmd)
                        cmd = '''babel -ismi tmp.smi -o%s %s -p 7.4 --gen3D --title %s > /dev/null 2>&1''' %(filename.split('.')[-1], filename, mol.GetPropsAsDict()['drugbank_id'])
                        os.system(cmd)
                    elif toolkit == 'rdk':
                        state = smi2sdf(smi, filename)
                        if state == 0:
                            filenames.remove(filename)
                            outputs.remove(output_name)
                    else:
                        raise Exception('只能选择rdk或ob')
                else:
                    filenames.remove(filename)
                    outputs.remove(output_name)
        except Exception as e:
            pass
    print('共分离出%s/%s个化合物' %(len(filenames), df.shape[0]))
    print('###############################')
    return filenames, outputs

# 定义得到smina命令行方法
def get_vina_command(rec, filename, output, center_x, center_y, center_z, size_x=20, size_y=20, size_z=20, parallel=8, dock='smina'):
    score_file = output.split('.')[0]+'_score.txt'
    if dock == 'smina':
        cmd = 'smina.static --receptor %s --ligand %s --center_x %s --center_y %s --center_z %s --size_x %s --size_y %s --size_z %s --out %s --cpu %s > %s' %(rec, filename, center_x, center_y, center_z, size_x, size_y, size_z, output, parallel, score_file)
    elif dock == 'qvina':
        convert_filename = filename.split('.')[0]+'.pdbqt'
        if not os.path.exists(convert_filename):
            cmd = 'babel -isdf %s -opdbqt %s > /dev/null 2>&1' %(filename, convert_filename)
            os.system(cmd)
        cmd = 'qvina-w --receptor %s --ligand %s --center_x %s --center_y %s --center_z %s --size_x %s --size_y %s --size_z %s --out %s --cpu %s > %s' %(rec, convert_filename, center_x, center_y, center_z, size_x, size_y, size_z, output, parallel, score_file)
    elif dock == 'vina':
        convert_filename = filename.split('.')[0]+'.pdbqt'
        if not os.path.exists(convert_filename):
            cmd = 'babel -isdf %s -opdbqt %s > /dev/null 2>&1' %(filename, convert_filename)
            os.system(cmd)
        cmd = 'vina --receptor %s --ligand %s --center_x %s --center_y %s --center_z %s --size_x %s --size_y %s --size_z %s --out %s --cpu %s > %s' %(rec, convert_filename, center_x, center_y, center_z, size_x, size_y, size_z, output, parallel, score_file)
    else:
        raise Exception('只能选择smina, qvina, vina三种对接方法！')
    return cmd

def run_vina(cmd):
    if not os.path.exists(cmd.split()[-1]):
        os.system(cmd)
    return cmd

def docking(rec, filenames, outputs, center_x, center_y, center_z, size_x=20, size_y=20, size_z=20, parallel=8, n_cpu=60, dock='smina'):
    print('')
    print('Step 5')
    print('###############################')
    print('正在对接，请耐心等待!!!')
    command_list = [get_vina_command(rec, filename, output, center_x, center_y, center_z, size_x, size_y, size_z, parallel, dock) for filename, output in zip(filenames, outputs)]
    m = int(n_cpu//parallel)
    with mp.Pool(m) as p:
        results = list(tqdm(p.imap(run_vina, command_list),total=len(command_list)))
    print('共完成%s个化合物对接!!!' %len(results))
    print('###############################')

def getScoreFromSmina(filename):
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        line_split = []
        for step, line in enumerate(lines):
            if '-----+------------+----------+----------' in line:
                line_split = lines[step+1]

        return float(line_split.split()[1])
    except Exception as e:
        return 999

def sumScore(df, key, outputs, output_file_name):
    print('')
    print('Step 6')
    print('###############################')
    print('正在统计化合物打分!!!')
    df_tmp = pd.DataFrame(columns=[key, 'score'])
    scores = []
    drugbank_id_list = []
    for output in outputs:
        output_result_name = output.split('.')[0] + '_score.txt'
        score = getScoreFromSmina(output_result_name)
        drugbank_id = output.split('/')[-1].split('_')[0]
        scores.append(score)
        drugbank_id_list.append(drugbank_id)
        
    df_tmp[key] = drugbank_id_list
    df_tmp['score'] = scores
    df_final = pd.merge(df,df_tmp, how='left', left_on=key, right_on=key)
    df_final.score.fillna(999, inplace=True)
    df_final.sort_values('score', ascending=True, inplace=True)
    print('正在生成化合物结构!!!')
    PandasTools.SaveXlsxFromFrame(df_final, output_file_name, molCol='mol')
    print('共得到%s个化合物得分!!!' %(df_final.shape[0]))
    print('###############################')
    return df_final

def getParser():
    parser = ArgumentParser()
    group = parser.add_argument_group('Input File Options')
    group.add_argument('-r', '--receptor', dest='receptor', metavar='<PDB_FILE>', required=True,
                    help='''受体蛋白文件，支持PDB和PDBQT格式''')
    group.add_argument('--db', dest='db', metavar='FILE', required=True,
                    help='''化合物数据库，支持SDF和CSV格式''')
    group.add_argument('--key_value', dest='key_value', metavar='STRING', required=True,
                    help='''化合物库分割命名''', default='name')
    group.add_argument('--input_path', dest='input_path', metavar='STRING',
                    help='''化合物库分割命名''', default='CPD')
    group.add_argument('--toolkit', dest='toolkit', metavar='STRING',
                    help='''使用rdkit(rdk)或openbabel(ob)生成化合物3D结构''', default='rdk')  
    group.add_argument('-c', '--crystal_ligand', dest='crystal_ligand', metavar='<PDB_FILE>',
                    help='''参考配体文件，支持PDB格式''')
    group.add_argument('--center_x', dest='center_x', metavar='FLOAT',
                    help='''对接口袋的x坐标''')
    group.add_argument('--center_y', dest='center_y', metavar='FLOAT',
                    help='''对接口袋的y坐标''')
    group.add_argument('--center_z', dest='center_z', metavar='FLOAT',
                    help='''对接口袋的z坐标''')
    group.add_argument('--size_x', dest='size_x', metavar='INT',
                    help='''口袋的x轴大小''', default=20)
    group.add_argument('--size_y', dest='size_y', metavar='INT',
                    help='''口袋的y轴大小''', default=20)
    group.add_argument('--size_z', dest='size_z', metavar='INT',
                    help='''口袋的z轴大小''', default=20)
    group.add_argument('--parallel', dest='parallel', metavar='INT',
                    help='''每个对接占用的CPU数量''', default=8)
    group.add_argument('--n_cpu', dest='n_cpu', metavar='INT',
                    help='''一共使用的CPU数量''', default=60)
    group.add_argument('--dock', dest='dock', metavar='STRING',
                    help='''可选择smina, vina, qvina''', default='smina')

    group = parser.add_argument_group('Output File Options')
    group.add_argument('--output_path', dest='output_path', metavar='FOLDER',
                    help='''数据库分割的储存地址''', default='output')
    group.add_argument('-o', '--output_file_name', dest='output_file_name', metavar='FILE',
                    help='''结果文件''', default='Final_Results.xlsx')
    opt = parser.parse_args()
    return opt

if __name__ == '__main__':

    opt = getParser()
    receptor_file = opt.receptor
    crystal_ligand = opt.crystal_ligand
    drug_db = opt.db
    key = opt.key_value
    toolkit = opt.toolkit

    size_x = opt.size_x
    size_y = opt.size_y
    size_z = opt.size_z

    parallel = opt.parallel
    n_cpu = opt.n_cpu
    dock = opt.dock

    input_path = opt.input_path
    output_path = opt.output_path

    # 生成受体pdbqt文件
    rec = prepareReceptor(receptor_file)

    if opt.crystal_ligand:
        center_x, center_y, center_z = getPocketsFromCrystalLigand(crystal_ligand)
    elif opt.center_x and opt.center_y and opt.center_z:
        center_x = opt.center_x
        center_y = opt.center_y
        center_z = opt.center_z
        print('')
        print('Step 2')
        print('###############################')
        print('选择%s,%s,%s作为口袋中心点!!!'%(center_x, center_y, center_z))
        print('###############################')
    else:
        raise Exception('必须输入参考化合物或口袋中心坐标!')

    # 载入化合物库，返回DataFrame
    df = loadDrugDatabase(drug_db)
    # 将分割化合物库，去除盐离子，生成3D构象
    filenames, outputs = spilt2SingleMolecule(df, input_path, output_path, toolkit, key)
    # 对接
    docking(rec, filenames, outputs, center_x, center_y, center_z, size_x, size_y, size_z, parallel, n_cpu, dock)
    # 结果处理
    sumScore(df, key, outputs, output_file_name)


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
