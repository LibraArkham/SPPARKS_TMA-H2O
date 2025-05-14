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
            coord = row['coord']
            pressureOn = row['pressureOn']
            expression = row['expression']

            # 构建基础命令
            if event_type == '1':
                cmd = f"event {event_type} {old1} {new1} {A} {n} {E} {coord} {pressureOn} {expression}"
            else:
                cmd = f"event {event_type} {old1} {new1} {old2} {new2} {A} {n} {E} {coord} {pressureOn} {expression}"
            cmd = cmd.strip()  # 去除末尾空格
            commands.append( (int(event_type), cmd) )

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
        print(f"成功生成命令文件：{output_file}")

    except FileNotFoundError:
        print(f"错误：Excel文件 {input_excel} 不存在！")
    except Exception as e:
        print(f"发生错误：{str(e)}")

# 调用示例
excel_to_events('events.xlsx', 'sorted_events.txt')