import pandas as pd

def excel_to_events(input_excel, output_file, sheet_name='Sheet1'):
    try:
        # 读取Excel文件并处理空值
        df = pd.read_excel(input_excel, sheet_name=sheet_name, dtype=str)
        df = df.fillna('')
        if 'expression' not in df.columns:
            df['expression'] = ''

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
            pressureOn = row['pressureOn']
            expression = row['expression']

            # 处理coord字段：支持单个值或逗号分隔的多个值
            if ',' in coord_str:
                # 多个coord值，分割并去除空格
                coord_values = [val.strip() for val in coord_str.split(',') if val.strip()]
            else:
                # 单个coord值
                coord_values = [coord_str] if coord_str else ['']

            # 为每个coord值生成一个事件
            for coord in coord_values:
                # 构建基础命令
                if event_type == '1':
                    cmd = f"event {event_type} {old1} {new1} {A} {n} {E} {coord} {pressureOn} {expression}"
                else:
                    cmd = f"event {event_type} {old1} {new1} {old2} {new2} {A} {n} {E} {coord} {pressureOn} {expression}"
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

def test_excel_to_events():
    """
    测试函数 - 创建示例Excel文件并测试多coord功能
    """
    import pandas as pd
    
    # 创建测试数据
    test_data = {
        'type': ['1', '2', '3', '4'],
        'old1': ['A', 'B', 'C', 'D'],
        'new1': ['A1', 'B1', 'C1', 'D1'],
        'old2': ['VAC', 'E', 'F', 'G'],
        'new2': ['VAC', 'E1', 'F1', 'G1'],
        'A': ['1.0e10', '2.0e10', '3.0e10', '4.0e10'],
        'n': ['0', '1', '2', '3'],
        'E': ['0.5', '1.0', '1.5', '2.0'],
        'coord': ['1,2,3', '4', '5,6', '7,8,9,10'],  # 测试多coord值
        'pressureOn': ['1', '2', '3', '4'],
        'expression': ['', 'test1', '', 'test2']
    }
    
    # 保存为Excel文件
    df = pd.DataFrame(test_data)
    df.to_excel('test_events.xlsx', index=False)
    print("已创建测试Excel文件：test_events.xlsx")
    
    # 测试转换
    excel_to_events('test_events.xlsx', 'test_output.txt')

# 主程序
if __name__ == "__main__":
    # 调用示例
    excel_to_events('events.xlsx', 'sorted_events.txt')
    
    # 如果要测试功能，可以取消下面的注释
    # test_excel_to_events()