# -*- coding: utf-8 -*-
"""
分子结构生成脚本
用于生成KMC模拟中使用的准确分子结构
"""

from ase import Atoms
from ase.io import write
from ase.build import molecule
import numpy as np
import os
from math import sin, cos, pi, radians

def create_trimethylaluminum():
    """
    创建三甲基铝 Al(CH3)3 分子
    使用准确的键长和键角
    """
    al_c_bond = 1.95  # Al-C键长 (Å)
    c_h_bond = 1.09   # C-H键长 (Å)
    
    positions = []
    symbols = []
    
    # 铝原子在中心
    positions.append([0, 0, 0])
    symbols.append('Al')
    
    # 三个甲基基团，平面三角形排列，120°间隔
    for i in range(3):
        angle = i * 2 * pi / 3  # 120度间隔
        
        # 碳原子位置
        c_x = al_c_bond * cos(angle)
        c_y = al_c_bond * sin(angle)
        c_z = 0
        positions.append([c_x, c_y, c_z])
        symbols.append('C')
        
        # 每个碳原子连接三个氢原子
        # 一个氢在z轴正方向，两个氢在下方形成四面体
        
        # H1: 沿着Al-C键方向稍微偏移
        h1_r = c_h_bond
        h1_angle = angle  # 与Al-C键同向
        h1_x = c_x + h1_r * cos(h1_angle) * 0.3
        h1_y = c_y + h1_r * sin(h1_angle) * 0.3
        h1_z = c_z + h1_r * 0.9
        positions.append([h1_x, h1_y, h1_z])
        symbols.append('H')
        
        # H2, H3: 在碳原子周围形成四面体结构
        for j in range(2):
            h_angle = angle + (j + 1) * 2 * pi / 3
            h_x = c_x + c_h_bond * cos(h_angle) * 0.6
            h_y = c_y + c_h_bond * sin(h_angle) * 0.6
            h_z = c_z - c_h_bond * 0.5
            positions.append([h_x, h_y, h_z])
            symbols.append('H')
    
    return Atoms(symbols=symbols, positions=positions)

def create_oh_al_methyl3():
    """
    创建 OH-Al(CH3)3 表面吸附复合物
    OH基团与三甲基铝结合
    """
    # 先创建TMA
    tma = create_trimethylaluminum()
    
    # 在铝原子上方添加OH基团
    oh_positions = [
        [0, 0, 1.8],     # O原子，Al-O键长约1.8Å
        [0, 0.96, 1.8]   # H原子，O-H键长约0.96Å
    ]
    oh_symbols = ['O', 'H']
    
    # 合并原子
    all_positions = list(tma.positions) + oh_positions
    all_symbols = list(tma.symbols) + oh_symbols
    
    return Atoms(symbols=all_symbols, positions=all_positions)

def create_o_al_methyl2():
    """
    创建 O-Al(CH3)2 (失去一个甲基的产物)
    """
    al_c_bond = 1.95
    al_o_bond = 1.8
    c_h_bond = 1.09
    
    positions = []
    symbols = []
    
    # 铝原子
    positions.append([0, 0, 0])
    symbols.append('Al')
    
    # 氧原子
    positions.append([0, 0, al_o_bond])
    symbols.append('O')
    
    # 两个甲基基团，夹角约120°
    for i in range(2):
        angle = i * 2 * pi / 3 + pi/6  # 起始角度偏移
        
        # 碳原子
        c_x = al_c_bond * cos(angle)
        c_y = al_c_bond * sin(angle)
        c_z = 0
        positions.append([c_x, c_y, c_z])
        symbols.append('C')
        
        # 氢原子
        for j in range(3):
            h_angle = angle + j * 2 * pi / 3
            h_x = c_x + c_h_bond * cos(h_angle) * 0.6
            h_y = c_y + c_h_bond * sin(h_angle) * 0.6
            h_z = c_z + (-1)**j * c_h_bond * 0.5
            positions.append([h_x, h_y, h_z])
            symbols.append('H')
    
    return Atoms(symbols=symbols, positions=positions)

def create_o_al_methyl1():
    """
    创建 O-Al(CH3) (失去两个甲基的产物)
    """
    al_c_bond = 1.95
    al_o_bond = 1.8
    c_h_bond = 1.09
    
    positions = []
    symbols = []
    
    # 铝原子
    positions.append([0, 0, 0])
    symbols.append('Al')
    
    # 氧原子
    positions.append([0, 0, al_o_bond])
    symbols.append('O')
    
    # 一个甲基基团
    c_x = al_c_bond
    c_y = 0
    c_z = 0
    positions.append([c_x, c_y, c_z])
    symbols.append('C')
    
    # 三个氢原子
    for i in range(3):
        angle = i * 2 * pi / 3
        h_x = c_x + c_h_bond * cos(angle) * 0.6
        h_y = c_y + c_h_bond * sin(angle) * 0.6
        h_z = c_z + (-1)**i * c_h_bond * 0.5
        positions.append([h_x, h_y, h_z])
        symbols.append('H')
    
    return Atoms(symbols=symbols, positions=positions)

def create_water_complexes():
    """
    创建含水的各种复合物
    """
    complexes = {}
    
    # OAlX2H2O: O-Al(CH3)2·H2O
    oalx2 = create_o_al_methyl2()
    h2o_pos = [[2.5, 0, 1.0], [3.2, 0.6, 1.0], [3.2, -0.6, 1.0]]  # H2O分子
    h2o_syms = ['O', 'H', 'H']
    
    all_pos = list(oalx2.positions) + h2o_pos
    all_syms = list(oalx2.symbols) + h2o_syms
    complexes['OAlX2H2O'] = Atoms(symbols=all_syms, positions=all_pos)
    
    # OAlXOH: O-Al(CH3)(OH)
    oalx = create_o_al_methyl1()
    oh_pos = [[0, 1.8, 0], [0, 2.76, 0]]  # OH基团
    oh_syms = ['O', 'H']
    
    all_pos = list(oalx.positions) + oh_pos
    all_syms = list(oalx.symbols) + oh_syms
    complexes['OAlXOH'] = Atoms(symbols=all_syms, positions=all_pos)
    
    return complexes

def create_hydroxyl_complexes():
    """
    创建各种羟基复合物
    """
    complexes = {}
    
    # OAlOH: O-Al-OH
    positions = [
        [0, 0, 0],      # Al
        [0, 0, 1.8],    # O (上方)
        [1.8, 0, 0],    # O (侧面)
        [2.76, 0, 0]    # H
    ]
    symbols = ['Al', 'O', 'O', 'H']
    complexes['OAlOH'] = Atoms(symbols=symbols, positions=positions)
    
    # OAlOH2: O-Al-(OH)2
    positions = [
        [0, 0, 0],       # Al
        [0, 0, 1.8],     # O (上方)
        [1.8, 0, 0],     # O (侧面1)
        [2.76, 0, 0],    # H
        [0, 1.8, 0],     # O (侧面2)
        [0, 2.76, 0]     # H
    ]
    symbols = ['Al', 'O', 'O', 'H', 'O', 'H']
    complexes['OAlOH2'] = Atoms(symbols=symbols, positions=positions)
    
    # AlOH: Al-OH
    positions = [
        [0, 0, 0],       # Al
        [1.8, 0, 0],     # O
        [2.76, 0, 0]     # H
    ]
    symbols = ['Al', 'O', 'H']
    complexes['AlOH'] = Atoms(symbols=symbols, positions=positions)
    
    # AlOH2: Al-(OH)2
    positions = [
        [0, 0, 0],       # Al
        [1.8, 0, 0],     # O1
        [2.76, 0, 0],    # H1
        [0, 1.8, 0],     # O2
        [0, 2.76, 0]     # H2
    ]
    symbols = ['Al', 'O', 'H', 'O', 'H']
    complexes['AlOH2'] = Atoms(symbols=symbols, positions=positions)
    
    return complexes

def generate_all_molecules():
    """
    生成所有需要的分子结构
    """
    molecules = {}
    
    print("生成基本分子...")
    # 基本分子
    molecules['O'] = Atoms('O', positions=[(0, 0, 0)])
    molecules['OH'] = Atoms('OH', positions=[(0, 0, 0), (0.96, 0, 0)])
    molecules['H2O'] = Atoms('H2O', positions=[
        (0, 0, 0), (0.76, 0.59, 0), (-0.76, 0.59, 0)
    ])
    molecules['Al'] = Atoms('Al', positions=[(0, 0, 0)])
    
    print("生成铝甲基化合物...")
    # 铝甲基化合物
    molecules['TMA'] = create_trimethylaluminum()  # Al(CH3)3
    molecules['OHAlX3'] = create_oh_al_methyl3()   # OH-Al(CH3)3
    molecules['OAlX2'] = create_o_al_methyl2()     # O-Al(CH3)2
    molecules['OAlX'] = create_o_al_methyl1()      # O-Al(CH3)
    molecules['OAl'] = Atoms('OAl', positions=[(0, 0, 0), (0, 0, 1.8)])
    
    print("生成水合物...")
    # 水合物
    water_complexes = create_water_complexes()
    molecules.update(water_complexes)
    
    print("生成羟基化合物...")
    # 羟基化合物
    hydroxyl_complexes = create_hydroxyl_complexes()
    molecules.update(hydroxyl_complexes)
    
    return molecules

def save_molecules(molecules, output_dir='./mol', format='xyz'):
    """
    保存分子到文件
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建目录: {output_dir}")
    
    saved_count = 0
    for name, mol in molecules.items():
        try:
            filename = f"{name}.{format}"
            filepath = os.path.join(output_dir, filename)
            write(filepath, mol, format=format)
            print(f"保存: {filename}")
            saved_count += 1
        except Exception as e:
            print(f"错误: 无法保存 {name}: {e}")
    
    print(f"总共保存了 {saved_count} 个分子文件")

def main():
    """
    主函数
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='生成KMC模拟用分子结构')
    parser.add_argument('--output', '-o', default='./mol', 
                        help='输出目录 (默认: ./mol)')
    parser.add_argument('--format', '-f', default='xyz', 
                        choices=['xyz', 'cif', 'pdb', 'json'],
                        help='输出格式 (默认: xyz)')
    
    args = parser.parse_args()
    
    print("开始生成分子结构...")
    molecules = generate_all_molecules()
    
    print(f"生成了 {len(molecules)} 个分子:")
    for name in sorted(molecules.keys()):
        print(f"  - {name}")
    
    print(f"\n保存分子到: {args.output}")
    save_molecules(molecules, args.output, args.format)
    
    print("\n完成!")

if __name__ == "__main__":
    main()