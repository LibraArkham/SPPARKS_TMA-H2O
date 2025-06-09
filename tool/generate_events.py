import pandas as pd

def excel_to_events(input_excel, output_file, sheet_name='Sheet1'):
    try:
        # 读取Excel文件并处理空值
        df = pd.read_excel(input_excel, sheet_name=sheet_name, dtype=str)
        df = df.fillna('')
        if 'expression' not in df.columns:
            df['expression'] = ''
        
        # 为Type I事件添加默认的coord2列（如果不存在）
        if 'coord2' not in df.columns:
            df['coord2'] = ''

        # 生成基础命令列表
        commands = []
        for _, row in df.iterrows():
            event_type = row['type']
            old1 = row['old1']
            new1 = row['new1']
            old2 = row.get('old2', 'VAC')
            new2 = row.get('new2', 'VAC')
            A = row['A']
            n = row['n']
            E = row['E']
            coord_str = str(row['coord']).strip()
            coord2_str = str(row.get('coord2', '')).strip()
            pressureOn = row['pressureOn']
            expression = row['expression']

            # 处理coord字段：支持单个值、逗号分隔的多个值、或"all"
            coord_values = parse_coord_field(coord_str)
            
            # 处理coord2字段（仅对Type II, III, IV有效）
            if event_type in ['2', '3', '4'] and coord2_str:
                coord2_values = parse_coord_field(coord2_str)
            else:
                coord2_values = ['']  # Type I事件或未指定coord2时使用空值

            # 生成所有coord和coord2的组合
            for coord in coord_values:
                for coord2 in coord2_values:
                    # 构建基础命令
                    if event_type == '1':
                        # Type I事件：只使用coord，忽略coord2
                        cmd = f"event {event_type} {old1} {new1} {A} {n} {E} {coord} {pressureOn} {expression}"
                    else:
                        # Type II, III, IV事件：使用coord和coord2
                        cmd = f"event {event_type} {old1} {new1} {old2} {new2} {A} {n} {E} {coord} {coord2} {pressureOn} {expression}"
                    
                    cmd = cmd.strip()  # 去除末尾空格
                    commands.append((int(event_type), cmd))

        # 按事件类型排序（1 > 2 > 3 > 4）
        sorted_commands = sorted(commands, key=lambda x: x[0])

        # 生成带编号的命令（编号在末尾）
        counters = {'1': 1, '2': 1, '3': 1, '4': 1}
        prefix_map = {'1': 's', '2': 'd', '3': 'v', '4': 'f'}
        final_commands = []
        
        for event_type, cmd in sorted_commands:
            type_str = str(event_type)
            prefix = prefix_map.get(type_str, 'x')
            count = counters[type_str]
            # 在命令末尾添加编号（自动处理末尾空格）
            numbered_cmd = f"{cmd} {prefix}{count}"
            final_commands.append(numbered_cmd)
            counters[type_str] += 1

        # 写入文件
        with open(output_file, 'w') as f:
            f.write('\n'.join(final_commands))
        
        # 打印统计信息
        total_events = len(final_commands)
        type_counts = {f"类型{k}": v-1 for k, v in counters.items()}
        print(f"成功生成命令文件：{output_file}")
        print(f"总事件数：{total_events}")
        print(f"各类型事件数量：{type_counts}")

    except FileNotFoundError:
        print(f"错误：Excel文件 {input_excel} 不存在！")
    except Exception as e:
        print(f"发生错误：{str(e)}")

def parse_coord_field(coord_str):
    """
    解析coord字段，支持：
    1. 单个数值：'1' -> ['1']
    2. 多个数值：'1,2,3' -> ['1', '2', '3']
    3. 'all'参数：'all' -> ['all']
    4. 混合格式：'1,all,3' -> ['1', 'all', '3']
    """
    if not coord_str or coord_str.lower() == 'nan':
        return ['']
    
    if ',' in coord_str:
        # 多个coord值，分割并去除空格
        coord_values = []
        for val in coord_str.split(','):
            val = val.strip()
            if val:
                # 检查是否为'all'参数（不区分大小写）
                if val.lower() == 'all':
                    coord_values.append('all')
                else:
                    coord_values.append(val)
        return coord_values if coord_values else ['']
    else:
        # 单个coord值
        if coord_str.lower() == 'all':
            return ['all']
        else:
            return [coord_str] if coord_str else ['']

def test_excel_to_events():
    """
    测试函数 - 创建示例Excel文件并测试新功能
    """
    import pandas as pd
    
    # 创建测试数据，包含各种coord和coord2组合
    test_data = {
        'type': ['1', '1', '2', '2', '3', '3', '4', '4'],
        'old1': ['OH', 'Ala', 'OH', 'OAlaX', 'OH', 'VAC', 'OAlaX2', 'OAlbX'],
        'new1': ['H2O', 'AlaOH', 'H2O', 'OAlaXOH', 'OHAlaX3', 'Ala', 'OAlaX', 'OAlbOH'],
        'old2': ['VAC', 'VAC', 'VAC', 'VAC', 'VAC', 'TMA', 'OAlbX2', 'OAlaX2'],
        'new2': ['VAC', 'VAC', 'TMA', 'H2O', 'VAC', 'VAC', 'OAlbX', 'OAlaX'],
        'A': ['1.0e10', '2.0e10', '3.0e10', '4.0e10', '5.0e10', '6.0e10', '7.0e10', '8.0e10'],
        'n': ['0', '1', '0', '1', '0', '1', '0', '1'],
        'E': ['0.5', '1.0', '0.3', '0.8', '0.4', '1.2', '0.6', '0.9'],
        'coord': ['1', 'all', '1,2', 'all', '1', '2,3', 'all', '1,all,3'],     # 测试各种coord格式
        'coord2': ['', '', 'all', '2', '1,2', 'all', '1', '2,all'],            # 测试coord2
        'pressureOn': ['1', '1', '1', '2', '1', '1', '2', '2'],
        'expression': ['', '', 'test1', '', 'test2', '', 'test3', 'test4']
    }
    
    # 保存为Excel文件
    df = pd.DataFrame(test_data)
    df.to_excel('test_events_new.xlsx', index=False)
    print("已创建测试Excel文件：test_events_new.xlsx")
    print("\n测试数据说明：")
    print("- Type 1事件：测试单coord和'all'参数")
    print("- Type 2事件：测试coord和coord2的各种组合")
    print("- Type 3事件：测试多coord值组合")
    print("- Type 4事件：测试混合格式（数字+all）")
    
    # 测试转换
    excel_to_events('test_events_new.xlsx', 'test_output_new.txt')

def validate_excel_format(input_excel, sheet_name='Sheet1'):
    """
    验证Excel文件格式是否正确
    """
    try:
        df = pd.read_excel(input_excel, sheet_name=sheet_name)
        
        required_columns = ['type', 'old1', 'new1', 'A', 'n', 'E', 'coord', 'pressureOn']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            print(f"错误：缺少必需的列：{missing_columns}")
            return False
        
        # 检查Type 2,3,4事件是否有old2和new2
        for _, row in df.iterrows():
            if row['type'] in ['2', '3', '4']:
                if 'old2' not in df.columns or 'new2' not in df.columns:
                    print(f"错误：Type {row['type']}事件需要old2和new2列")
                    return False
        
        print("Excel文件格式验证通过！")
        return True
        
    except Exception as e:
        print(f"验证过程中发生错误：{str(e)}")
        return False

# 主程序
if __name__ == "__main__":
    # 验证并转换事件文件
    input_file = 'events.xlsx'
    output_file = 'sorted_events.txt'
    
    print("开始验证Excel文件格式...")
    if validate_excel_format(input_file):
        print(f"开始转换事件文件：{input_file} -> {output_file}")
        excel_to_events(input_file, output_file)
    else:
        print("Excel文件格式不正确，请检查后重试。")
    
    # 如果要测试功能，可以取消下面的注释
    # print("\n" + "="*50)
    # print("运行测试功能...")
    # test_excel_to_events()