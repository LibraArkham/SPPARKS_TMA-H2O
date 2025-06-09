# -*- coding: utf-8 -*-

from ase import Atoms
from ase.io import read, write
import numpy as np
import pandas as pd
from ase.visualize import view
from copy import deepcopy
import os
import argparse
import glob

def translation(row):
    """将分子移动到指定位置"""
    new_position = np.array(row['position'])
    atoms = deepcopy(row['ligand'][0])
    center_atom_index = row['ligand'][1]
    vector = new_position - atoms.positions[center_atom_index]
    atoms.positions = atoms.positions + vector
    return atoms

def run(file_path, dict_ligand, lattice_vector, output_dir='./xyz', time_range=None):
    """处理KMC dump文件并生成可视化文件"""
    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 解析时间范围参数
    start_time = 0
    end_time = float('inf')
    time_step = 1
    
    if time_range is not None:
        if len(time_range) >= 1:
            start_time = time_range[0]
        if len(time_range) >= 2:
            end_time = time_range[1]
        if len(time_range) >= 3:
            time_step = time_range[2]
    
    print(f"时间范围设置: 开始={start_time}, 结束={end_time}, 步长={time_step}")
    
    with open(file_path, 'r', encoding='utf-8') as file:
        flag = 0
        row = 0
        natoms = -1
        time = -1
        data = []
        saved_count = 0
        
        for i, line in enumerate(file):
            line = line.strip()
            
            if i == 3:
                natoms = int(line)
                continue
        
            if 'ITEM: ATOMS' in line:
                time += 1
                flag = 1
                row = 0
                columns = line.split()[2:]
                
                if time < start_time:
                    continue
                if time > end_time:
                    break
                if (time - start_time) % time_step != 0:
                    continue
                    
                print(f"处理时间步 {time}...")
                continue
            
            if time < start_time or time > end_time or (time - start_time) % time_step != 0:
                continue
            
            if flag == 1:
                str_list = line.split()
                values = [float(v) for v in str_list]
                data.append(values)
                row += 1
            
            if row == natoms:
                df = pd.DataFrame(data)
                df.columns = columns
                flag = 0
                row = 0
                
                atoms = get_atoms(dict_ligand, df, lattice_vector)
                output_file = os.path.join(output_dir, f'frame_{time:06d}.cif')
                write(output_file, atoms)
                saved_count += 1
                
                del data
                data = []
        
        print(f"总共保存了 {saved_count} 个时间步的文件")

def get_atoms(dict_ligand, df, lattice_vector):
    """根据KMC数据构建原子结构"""
    df['position'] = ((df['x'].values[:, np.newaxis] * lattice_vector[0]) + 
                      (df['y'].values[:, np.newaxis] * lattice_vector[1]) + 
                      (df['z'].values[:, np.newaxis] * lattice_vector[2])).tolist()
    
    df['ligand'] = df['i1'].map(dict_ligand)
    df = df.dropna()
    df['all_atoms'] = df.apply(translation, axis=1)
    
    atoms = df['all_atoms'].sum()  
    atoms.set_cell(lattice_vector)
    atoms.set_pbc(True)
    
    return atoms

def load_molecules_from_directory(mol_dir='./mol'):
    """从目录中自动加载所有分子文件"""
    molecules = {}
    
    if not os.path.exists(mol_dir):
        print(f"错误: 分子目录 {mol_dir} 不存在!")
        return molecules
    
    # 查找所有支持的分子文件格式
    supported_formats = ['*.xyz', '*.mol', '*.cif', '*.pdb', '*.json']
    mol_files = []
    
    for fmt in supported_formats:
        mol_files.extend(glob.glob(os.path.join(mol_dir, fmt)))
    
    if not mol_files:
        print(f"警告: 在 {mol_dir} 中未找到任何分子文件")
        return molecules
    
    for filepath in mol_files:
        filename = os.path.basename(filepath)
        name = os.path.splitext(filename)[0]  # 去掉文件扩展名
        
        try:
            mol = read(filepath)
            molecules[name] = mol
            print(f"加载分子: {name} <- {filename}")
        except Exception as e:
            print(f"警告: 无法加载 {filepath}: {e}")
    
    print(f"总共加载了 {len(molecules)} 个分子")
    return molecules

def create_species_mapping(molecules):
    """根据KMC代码中的enum定义创建物种映射"""
    # enum{VACANCY,O,OH,Ala,OHAlaX3,OAlaX2,OAlaX2H2O,OAlaXOH,OAlaX,OAlaOH,OAlaOH2,AlaOH,AlaOH2,Alb,OHAlbX3,OAlbX2,OAlbX2H2O,OAlbXOH,OAlbX,OAlbOH,OAlbOH2,AlbOH,AlbOH2,OAla,OAlb,H2O}
    
    # 物种名称映射表（将KMC中的名称映射到分子文件名）
    species_to_molecule = {
        1: 'O',           # O
        2: 'OH',          # OH
        3: 'Al',          # Ala -> Al (统一使用Al)
        4: 'OHAlX3',      # OHAlaX3
        5: 'OAlX2',       # OAlaX2
        6: 'OAlX2H2O',    # OAlaX2H2O
        7: 'OAlXOH',      # OAlaXOH
        8: 'OAlX',        # OAlaX
        9: 'OAlOH',       # OAlaOH
        10: 'OAlOH2',     # OAlaOH2
        11: 'AlOH',       # AlaOH -> AlOH
        12: 'AlOH2',      # AlaOH2 -> AlOH2
        13: 'Al',         # Alb -> Al (统一使用Al)
        14: 'OHAlX3',     # OHAlbX3 -> 使用相同的OHAlX3
        15: 'OAlX2',      # OAlbX2 -> 使用相同的OAlX2
        16: 'OAlX2H2O',   # OAlbX2H2O -> 使用相同的OAlX2H2O
        17: 'OAlXOH',     # OAlbXOH -> 使用相同的OAlXOH
        18: 'OAlX',       # OAlbX -> 使用相同的OAlX
        19: 'OAlOH',      # OAlbOH -> 使用相同的OAlOH
        20: 'OAlOH2',     # OAlbOH2 -> 使用相同的OAlOH2
        21: 'AlOH',       # AlbOH -> 使用相同的AlOH
        22: 'AlOH2',      # AlbOH2 -> 使用相同的AlOH2
        23: 'OAl',        # OAla
        24: 'OAl',        # OAlb -> 使用相同的OAl
        25: 'H2O',        # H2O
    }
    
    dict_ligand = {}
    missing_molecules = []
    
    for species_id, mol_name in species_to_molecule.items():
        if mol_name in molecules:
            # 假设每个分子的第一个原子是中心原子
            dict_ligand[species_id] = (molecules[mol_name], 0)
        else:
            missing_molecules.append(mol_name)
            # 使用单个原子作为备用
            dict_ligand[species_id] = (Atoms('X', positions=[(0, 0, 0)]), 0)
    
    if missing_molecules:
        print(f"警告: 以下分子文件未找到: {missing_molecules}")
        print("将使用占位符原子代替")
    
    return dict_ligand

def parse_time_range(time_str):
    """解析时间范围字符串"""
    if not time_str:
        return None
    
    parts = time_str.split(':')
    
    if len(parts) == 1:
        return (int(parts[0]), float('inf'), 1)
    
    start = 0 if parts[0] == '' else int(parts[0])
    end = float('inf') if parts[1] == '' else int(parts[1])
    step = 1 if len(parts) < 3 or parts[2] == '' else int(parts[2])
    
    return (start, end, step)

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='KMC模拟结果可视化工具')
    parser.add_argument('--dump', '-d', default='./dump.ald', 
                        help='dump文件路径 (默认: ./dump.ald)')
    parser.add_argument('--output', '-o', default='./xyz', 
                        help='输出目录 (默认: ./xyz)')
    parser.add_argument('--time', '-t', default=None, 
                        help='时间范围，格式: start:end:step (例如: 100:500:10)')
    parser.add_argument('--mol-dir', '-m', default='./mol',
                        help='分子文件目录 (默认: ./mol)')
    parser.add_argument('--lattice', '-l', nargs=9, type=float,
                        default=[24.0251361, 0, 0, -12.012568049999997, 20.806378191978606, 0, 0, 0, 65.08907438726766],
                        help='晶格向量 (9个数字，按行展开)')
    
    args = parser.parse_args()
    
    # 解析时间范围
    time_range = parse_time_range(args.time)
    
    # 重构晶格向量
    lattice_vector = np.array(args.lattice).reshape(3, 3)
    
    # 加载分子
    print(f"从目录加载分子: {args.mol_dir}")
    molecules = load_molecules_from_directory(args.mol_dir)
    
    if not molecules:
        print("错误: 未加载到任何分子，请检查分子文件目录")
        return
    
    # 创建物种映射
    dict_ligand = create_species_mapping(molecules)
    
    # 处理dump文件
    if os.path.exists(args.dump):
        print(f"开始处理dump文件: {args.dump}")
        if time_range:
            print(f"时间范围: {time_range[0]}:{time_range[1]}:{time_range[2]}")
        else:
            print("处理全部时间步")
        
        run(args.dump, dict_ligand, lattice_vector, args.output, time_range)
        print(f"可视化文件已保存到: {args.output}")
        print("使用VESTA, Ovito或其他软件打开CIF文件进行可视化")
    else:
        print(f"错误: 找不到dump文件 {args.dump}")
        print("请确保SPPARKS生成了dump文件")

if __name__ == "__main__":
    main()