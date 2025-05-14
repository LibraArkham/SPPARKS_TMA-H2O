import graphviz as gv
import pandas as pd
import os
import datetime
import hashlib

def parse_reactions(file_path):
    reactions = []
    unique_reactions = {}  # 用于存储唯一反应
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('event'):
                parts = line.strip().split()
                if len(parts) >= 8:
                    event_type = parts[1]
                    
                    # 根据不同的事件类型正确解析反应物、产物和能量
                    if event_type == '1':
                        # Type 1: 单分子反应
                        reactant = parts[2]
                        product = parts[3]
                        energy = parts[6]  # 能量障碍
                        label = parts[-1] if len(parts) > 8 else ""
                        reaction_key = f"{event_type}_{reactant}_{product}_{energy}"
                        reaction_data = (event_type, reactant, product, energy, label)
                    elif event_type in ['2', '3', '4']:
                        # Type 2,3,4: 多分子反应
                        reactant = parts[2]
                        reactant2 = parts[4]
                        product = parts[3]
                        product2 = parts[5]
                        energy = parts[8]  # 修正: 使用正确的能量位置 (第9个参数)
                        label = parts[-1] if len(parts) > 10 else ""
                        reaction_key = f"{event_type}_{reactant}_{reactant2}_{product}_{product2}_{energy}"
                        reaction_data = (event_type, f"{reactant}+{reactant2}", f"{product}+{product2}", energy, label)
                    
                    # 如果这是一个新的唯一反应，添加到字典中
                    if reaction_key not in unique_reactions:
                        unique_reactions[reaction_key] = reaction_data
                        reactions.append(reaction_data)
    
    print(f"解析了 {len(reactions)} 个唯一反应（忽略坐标差异）")
    return reactions

def create_reaction_diagram(reactions):
    # 提高分辨率设置，横版布局
    g = gv.Digraph(format='png')
    g.attr(rankdir='LR', size='30,20', dpi='1200')  # 横版布局，更大尺寸和更高DPI
    g.attr('node', shape='box', style='filled', fillcolor='lightblue', fontname='Arial', fontsize='16')
    g.attr('edge', fontname='Arial', fontsize='14', penwidth='2.0')
    
    # 按反应类型分组
    type1 = [r for r in reactions if r[0] == '1']
    type2 = [r for r in reactions if r[0] == '2']
    type3 = [r for r in reactions if r[0] == '3']
    type4 = [r for r in reactions if r[0] == '4']
    
    # 添加子图和节点
    with g.subgraph(name='cluster_type1') as c:
        c.attr(label='Type 1', style='filled', fillcolor='lightyellow', fontsize='20')
        for r in type1:
            c.node(r[1])
            c.node(r[2])
            c.edge(r[1], r[2], label=f'E={r[3]} eV\n{r[4]}')
    
    with g.subgraph(name='cluster_type2') as c:
        c.attr(label='Type 2', style='filled', fillcolor='lightgreen', fontsize='20')
        for r in type2:
            c.node(r[1])
            c.node(r[2])
            c.edge(r[1], r[2], label=f'E={r[3]} eV\n{r[4]}')
    
    with g.subgraph(name='cluster_type3') as c:
        c.attr(label='Type 3', style='filled', fillcolor='lightpink', fontsize='20')
        for r in type3:
            c.node(r[1])
            c.node(r[2])
            c.edge(r[1], r[2], label=f'E={r[3]} eV\n{r[4]}')
    
    with g.subgraph(name='cluster_type4') as c:
        c.attr(label='Type 4', style='filled', fillcolor='lightcyan', fontsize='20')
        for r in type4:
            c.node(r[1])
            c.node(r[2])
            c.edge(r[1], r[2], label=f'E={r[3]} eV\n{r[4]}')
    
    return g

def create_oh_cycle_diagram(reactions):
    """创建从OH出发的物种循环图"""
    g = gv.Digraph(format='png')
    g.attr(rankdir='LR', size='40,20', dpi='1200')  # 横版布局，更大尺寸和更高DPI
    g.attr('node', shape='box', style='filled', fontname='Arial', fontsize='18')
    g.attr('edge', fontname='Arial', fontsize='16', penwidth='2.0')
    
    # 创建一个有向图来表示反应网络
    reaction_network = {}
    
    # 添加所有反应到网络
    for r in reactions:
        event_type = r[0]
        
        if event_type == '1':
            # 单分子反应
            reactant = r[1]
            product = r[2]
            energy = r[3]
            label = r[4]
            
            if reactant not in reaction_network:
                reaction_network[reactant] = []
            reaction_network[reactant].append((product, energy, label))
            
        elif event_type in ['2', '3', '4']:
            # 多分子反应，需要解析复合物
            reactants = r[1].split('+')
            products = r[2].split('+')
            energy = r[3]
            label = r[4]
            
            for reactant in reactants:
                if reactant not in reaction_network:
                    reaction_network[reactant] = []
                
                for product in products:
                    reaction_network[reactant].append((product, energy, label))
    
    # 找出所有包含OH的反应
    oh_reactions = []
    o_reactions = []
    
    # 添加OH和O节点，设置特殊颜色
    g.node('OH', fillcolor='gold', style='filled')
    g.node('O', fillcolor='gold', style='filled')
    
    # 从OH开始的路径
    visited = set()
    
    def dfs_path(current, path=None, depth=0):
        if path is None:
            path = []
        
        if depth > 10:  # 限制深度以避免无限循环
            return
        
        if current in visited:
            return
        
        visited.add(current)
        new_path = path + [current]
        
        # 如果回到了O或OH，并且路径长度大于1，则添加这条路径
        if (current == 'O' or current == 'OH') and len(new_path) > 2:
            cycle_paths.append(new_path)
            return
        
        # 继续搜索
        if current in reaction_network:
            for next_species, energy, label in reaction_network[current]:
                if next_species not in new_path or next_species in ['O', 'OH']:
                    dfs_path(next_species, new_path, depth + 1)
        
        visited.remove(current)
    
    # 存储找到的循环路径
    cycle_paths = []
    
    # 从OH开始搜索
    dfs_path('OH')
    
    # 为每个找到的循环添加节点和边
    added_edges = set()  # 避免重复添加边
    
    for path in cycle_paths:
        for i in range(len(path) - 1):
            current = path[i]
            next_species = path[i + 1]
            
            # 为节点设置颜色
            if current not in ['O', 'OH']:
                g.node(current, fillcolor='lightblue', style='filled')
            
            if next_species not in ['O', 'OH']:
                g.node(next_species, fillcolor='lightblue', style='filled')
            
            # 查找这条边的能量和标签
            energy = "?"
            label = ""
            
            if current in reaction_network:
                for target, e, l in reaction_network[current]:
                    if target == next_species:
                        energy = e
                        label = l
                        break
            
            # 添加边（如果尚未添加）
            edge_key = f"{current}->{next_species}"
            if edge_key not in added_edges:
                g.edge(current, next_species, label=f'E={energy} eV\n{label}')
                added_edges.add(edge_key)
    
    # 如果没有找到循环路径，添加一个提示节点
    if not cycle_paths:
        g.node('No_Cycles_Found', label='未找到从OH到O/OH的循环路径', shape='plaintext', fontsize='24')
    
    return g

def create_output_directory():
    """创建带有时间戳的输出目录"""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"reaction_diagrams_{timestamp}"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建输出目录: {output_dir}")
    
    return output_dir, timestamp

def main():
    file_path = 'in.ald'
    
    # 创建输出目录
    output_dir, timestamp = create_output_directory()
    
    # 解析反应
    reactions = parse_reactions(file_path)
    
    # 创建反应图
    diagram = create_reaction_diagram(reactions)
    
    # 设置输出文件名和选项
    output_file = os.path.join(output_dir, f'ald_reaction_diagram_{timestamp}')
    diagram.render(output_file, view=True, cleanup=True)
    
    # 创建OH循环图
    oh_cycle_diagram = create_oh_cycle_diagram(reactions)
    oh_cycle_file = os.path.join(output_dir, f'oh_cycle_diagram_{timestamp}')
    oh_cycle_diagram.render(oh_cycle_file, view=True, cleanup=True)
    
    # 保存反应数据到CSV文件
    reaction_df = pd.DataFrame(reactions, columns=['Type', 'Reactants', 'Products', 'Energy', 'Label'])
    csv_file = os.path.join(output_dir, f'reactions_data_{timestamp}.csv')
    reaction_df.to_csv(csv_file, index=False)
    
    print(f"反应路径图已生成: {output_file}.png")
    print(f"OH循环图已生成: {oh_cycle_file}.png")
    print(f"反应数据已保存: {csv_file}")
    print(f"图像分辨率: 1200 DPI")
    print(f"所有文件已保存到目录: {output_dir}")

if __name__ == '__main__':
    main()