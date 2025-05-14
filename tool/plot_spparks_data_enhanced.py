import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import re
import os
import json
import argparse
import time
from datetime import datetime
from tqdm import tqdm  # Import tqdm for progress bars
from multiprocessing import Pool, cpu_count
import functools
import signal
import sys
import atexit

# 导入Dash相关库
import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc

def parse_model_dimensions(model_file_path):
    """Parse model dimensions from the ALD data file"""
    print("Parsing model dimensions...")
    
    # Default values in case parsing fails
    x_dim = 1.0
    y_dim = 1.0
    
    try:
        with open(model_file_path, 'r') as f:
            lines = f.readlines()
            
            for i, line in enumerate(lines):
                if "xlo xhi" in line:
                    parts = line.strip().split()
                    x_dim = float(parts[1]) - float(parts[0])
                    # Extract vector information if available
                    if "#" in line and "[" in line:
                        vector_info = line.split("#")[1].strip()
                        vector_values = re.findall(r"[-+]?\d*\.\d+|\d+", vector_info)
                        if len(vector_values) >= 1:
                            x_dim = float(vector_values[0])
                
                if "ylo yhi" in line:
                    parts = line.strip().split()
                    y_dim = float(parts[1]) - float(parts[0])
                    # Extract vector information if available
                    if "#" in line and "[" in line:
                        vector_info = line.split("#")[1].strip()
                        vector_values = re.findall(r"[-+]?\d*\.\d+|\d+", vector_info)
                        if len(vector_values) >= 2:
                            y_dim = float(vector_values[1])
        
        # Convert from Angstrom to cm
        area_cm2 = (x_dim * y_dim) / (10**16)  # 1 Å² = 10^-16 cm²
        
        print(f"Model dimensions: X = {x_dim} Å, Y = {y_dim} Å")
        print(f"Surface area: {area_cm2} cm²")
        
        return x_dim, y_dim, area_cm2
    
    except Exception as e:
        print(f"Error parsing model dimensions: {str(e)}")
        # Return default values
        return x_dim, y_dim, (x_dim * y_dim) / (10**16)

def parse_log_file(file_path, chunk_size=10000):
    """Parse SPPARKS log file and extract data - Optimized version"""
    print("Parsing log file...")
    
    # Find the line with headers and the start line for data
    header_line = None
    start_line = 0
    
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if "Time" in line and "Naccept" in line:
                header_line = line
                start_line = i + 1
                break
    
    if header_line is None:
        raise ValueError("Could not find data header in log file")
    
    # Parse column headers
    headers = header_line.strip().split()
    
    # Pre-allocate memory, collect all row data
    all_rows = []
    
    # Create progress bar
    print("Processing data rows...")
    
    # Calculate total lines for progress bar
    with open(file_path, 'r') as f:
        total_lines = sum(1 for _ in f) - start_line
    
    # Read file in chunks
    with open(file_path, 'r') as f:
        # Skip already read lines
        for _ in range(start_line):
            next(f)
        
        # Process in chunks
        chunk = []
        for line in tqdm(f, total=total_lines, desc="Parsing data", unit="rows"):
            # Skip non-data lines
            if not re.match(r'\s*[\d\.]+\s+\d+', line):
                continue
            
            # Split data
            values = line.strip().split()
            if len(values) < len(headers):
                continue
            
            # Add data to row list
            row_data = {}
            for i, header in enumerate(headers):
                if i < len(values):
                    try:
                        # Try to convert to numeric type
                        if '.' in values[i]:
                            row_data[header] = float(values[i])
                        else:
                            row_data[header] = int(values[i])
                    except ValueError:
                        row_data[header] = values[i]
            
            chunk.append(row_data)
            
            # When chunk reaches specified size, add to all_rows and clear chunk
            if len(chunk) >= chunk_size:
                all_rows.extend(chunk)
                chunk = []
        
        # Process the last chunk
        if chunk:
            all_rows.extend(chunk)
    
    # Create DataFrame at once
    print("Creating DataFrame...")
    data = pd.DataFrame(all_rows)
    
    # Ensure numeric columns have correct types
    numeric_cols = ['Time', 'Naccept', 'Nreject', 'Nsweeps', 'CPU', 'events', 'QCM']
    for col in numeric_cols:
        if col in data.columns:
            # Try to convert column to numeric type, set errors to NaN
            data[col] = pd.to_numeric(data[col], errors='coerce')
            # Replace NaN values with 0
            data[col] = data[col].fillna(0)
    
    return data

def convert_qcm_to_ng_cm2(data, area_cm2):
    """Convert QCM from atomic mass units to ng/cm²"""
    if 'QCM' in data.columns:
        # Atomic mass unit (amu) to grams: 1 amu = 1.66053886 × 10^-24 grams
        amu_to_g = 1.66053886e-24
        
        # Convert to ng/cm²: (QCM * amu_to_g * 10^9) / area_cm2
        # 10^9 converts g to ng
        data['QCM_ng_cm2'] = (data['QCM'] * amu_to_g * 1e9) / area_cm2
        
        print(f"Converted QCM to ng/cm²")
        print(f"QCM range: {data['QCM_ng_cm2'].min()} to {data['QCM_ng_cm2'].max()} ng/cm²")
    
    return data

def create_output_directory(base_dir):
    """Create time-categorized output directory"""
    # Create timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create new output directory
    output_dir = os.path.join(base_dir, f"run_{timestamp}")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    # Create subdirectories
    plots_dir = os.path.join(output_dir, "plots")
    exports_dir = os.path.join(output_dir, "exports")
    
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
    if not os.path.exists(exports_dir):
        os.makedirs(exports_dir)
    
    return output_dir, plots_dir, exports_dir, timestamp

def write_html_with_config(fig, file_path):
    """Write figure to HTML with enhanced configuration options"""
    config = {
        'editable': True,  # 允许编辑图表
        'toImageButtonOptions': {
            'format': 'png',  # 导出图像格式
            'filename': 'custom_image',
            'height': 800,
            'width': 1200,
            'scale': 2  # 高分辨率
        },
        'modeBarButtonsToAdd': [
            'drawline',
            'drawopenpath',
            'drawclosedpath',
            'drawcircle',
            'drawrect',
            'eraseshape'
        ],
        'displaylogo': False,  # 隐藏Plotly logo
        'responsive': True,    # 响应式布局
    }
    
    fig.write_html(file_path, config=config)
    return file_path

def add_style_controls(fig):
    """添加交互式样式控制到图表"""
    # 添加按钮和滑块等控件来控制图表样式
    # 这里提供一个简单的实现，您可以根据需要扩展
    
    # 添加范围选择器
    fig.update_layout(
        xaxis=dict(
            rangeslider=dict(visible=True),
            type="linear"
        )
    )
    
    # 添加注释按钮
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                active=0,
                x=0.1,
                y=1.15,
                buttons=[
                    dict(
                        label="线性比例",
                        method="relayout",
                        args=[{"yaxis.type": "linear"}]
                    ),
                    dict(
                        label="对数比例",
                        method="relayout",
                        args=[{"yaxis.type": "log"}]
                    )
                ]
            )
        ]
    )
    
    return fig

# 修改图表保存函数
def create_events_plot(data, plots_dir, timestamp, style_config):
    """Create interactive chart of event count over time"""
    # Use complete data
    plot_data = data
    
    fig_events = go.Figure()
    fig_events.add_trace(go.Scatter(
        x=plot_data['Time'], 
        y=plot_data['Naccept'], 
        mode='lines+markers',
        name='Accepted Events',
        line=dict(color='blue', width=2),
        marker=dict(size=5)
    ))
    
    fig_events.update_layout(
        title='Accepted Events vs Time',
        xaxis_title='Time',
        yaxis_title='Event Count',
        template='plotly_white',
        hovermode='x unified'
    )
    
    # 添加交互式样式控件
    fig_events = add_style_controls(fig_events)
    
    # Save as HTML file
    events_path = os.path.join(plots_dir, f'events_vs_time_{timestamp}.html')
    write_html_with_config(fig_events, events_path)
    return events_path

def create_qcm_plot(data, plots_dir, timestamp, style_config=None):
    """Create interactive chart of QCM mass over time"""
    if 'QCM_ng_cm2' not in data.columns:
        return None
    
    # Use complete data
    plot_data = data
    
    fig_qcm = go.Figure()
    fig_qcm.add_trace(go.Scatter(
        x=plot_data['Time'], 
        y=plot_data['QCM_ng_cm2'], 
        mode='lines+markers',
        name='QCM Mass',
        line=dict(color='purple', width=2),
        marker=dict(size=5)
    ))
    
    fig_qcm.update_layout(
        title='QCM Mass vs Time',
        xaxis_title='Time',
        yaxis_title='Mass (ng/cm²)',
        template='plotly_white',
        hovermode='x unified'
    )
    
    # 添加交互式样式控件
    fig_qcm = add_style_controls(fig_qcm)
    
    # Save as HTML file
    qcm_path = os.path.join(plots_dir, f'qcm_vs_time_{timestamp}.html')
    write_html_with_config(fig_qcm, qcm_path)
    return qcm_path

def create_cpu_plot(data, plots_dir, timestamp, style_config=None):
    """Create interactive chart of CPU time vs simulation time"""
    # Use complete data
    plot_data = data
    
    fig_cpu = go.Figure()
    fig_cpu.add_trace(go.Scatter(
        x=plot_data['Time'], 
        y=plot_data['CPU'], 
        mode='lines+markers',
        name='CPU Time',
        line=dict(color='red', width=2),
        marker=dict(size=5)
    ))
    
    fig_cpu.update_layout(
        title='CPU Time vs Simulation Time',
        xaxis_title='Simulation Time',
        yaxis_title='CPU Time (seconds)',
        template='plotly_white',
        hovermode='x unified'
    )
    
    # 添加交互式样式控件
    fig_cpu = add_style_controls(fig_cpu)
    
    # Save as HTML file
    cpu_path = os.path.join(plots_dir, f'cpu_vs_time_{timestamp}.html')
    write_html_with_config(fig_cpu, cpu_path)
    return cpu_path

def create_main_elements_plot(data, plots_dir, timestamp, style_config=None):
    """Create interactive chart of main element counts over time"""
    # Get all element columns (exclude non-element data columns)
    exclude_cols = ['Time', 'Naccept', 'Nreject', 'Nsweeps', 'CPU', 'events', 'QCM', 'QCM_ng_cm2']
    element_cols = [col for col in data.columns if col not in exclude_cols]
    
    # Main elements
    main_elements = ['VAC', 'O', 'OH', 'Ala', 'Alb', 'OHAlaX3', 'OHAlbX3']
    main_elements = [elem for elem in main_elements if elem in data.columns]
    
    # Use complete data
    plot_data = data
    
    fig_main = go.Figure()
    for element in main_elements:
        fig_main.add_trace(go.Scatter(
            x=plot_data['Time'], 
            y=plot_data[element], 
            mode='lines',
            name=element
        ))
    
    fig_main.update_layout(
        title='Main Element Counts vs Time',
        xaxis_title='Time',
        yaxis_title='Element Count',
        template='plotly_white',
        hovermode='x unified',
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1)
    )
    
    # 添加交互式样式控件
    fig_main = add_style_controls(fig_main)
    
    # Save as HTML file
    main_path = os.path.join(plots_dir, f'main_elements_vs_time_{timestamp}.html')
    write_html_with_config(fig_main, main_path)
    return main_path

def create_dashboard(data, plots_dir, timestamp, style_config=None):
    """Create comprehensive dashboard"""
    # Get all element columns (exclude non-element data columns)
    exclude_cols = ['Time', 'Naccept', 'Nreject', 'Nsweeps', 'CPU', 'events', 'QCM', 'QCM_ng_cm2']
    element_cols = [col for col in data.columns if col not in exclude_cols]
    
    # Main elements
    main_elements = ['VAC', 'O', 'OH', 'Ala', 'Alb', 'OHAlaX3', 'OHAlbX3']
    main_elements = [elem for elem in main_elements if elem in data.columns]
    
    # Ala related elements
    ala_elements = [col for col in element_cols if 'Ala' in col]
    
    # Alb related elements
    alb_elements = [col for col in element_cols if 'Alb' in col]
    
    # Use complete data
    plot_data = data
    
    fig_dashboard = make_subplots(
        rows=4, cols=2,
        subplot_titles=(
            'Accepted Events vs Time', 'CPU Time vs Simulation Time',
            'QCM Mass vs Time', 'Main Element Counts vs Time',
            'O and OH Element Counts vs Time', 'Ala Related Elements',
            'Alb Related Elements', ''
        ),
        vertical_spacing=0.1,
        horizontal_spacing=0.05
    )
    
    # Add event count plot
    fig_dashboard.add_trace(
        go.Scatter(x=plot_data['Time'], y=plot_data['Naccept'], mode='lines', name='Accepted Events'),
        row=1, col=1
    )
    
    # Add CPU time plot
    fig_dashboard.add_trace(
        go.Scatter(x=plot_data['Time'], y=plot_data['CPU'], mode='lines', name='CPU Time', line=dict(color='red')),
        row=1, col=2
    )
    
    # Add QCM mass plot
    if 'QCM_ng_cm2' in plot_data.columns:
        fig_dashboard.add_trace(
            go.Scatter(x=plot_data['Time'], y=plot_data['QCM_ng_cm2'], mode='lines', name='QCM Mass (ng/cm²)', line=dict(color='purple')),
            row=2, col=1
        )
    
    # Add main elements plot
    for element in main_elements:
        fig_dashboard.add_trace(
            go.Scatter(x=plot_data['Time'], y=plot_data[element], mode='lines', name=element),
            row=2, col=2
        )
    
    # Add O and OH elements plot
    if 'O' in plot_data.columns and 'OH' in plot_data.columns:
        fig_dashboard.add_trace(
            go.Scatter(x=plot_data['Time'], y=plot_data['O'], mode='lines', name='O', line=dict(color='green')),
            row=3, col=1
        )
        fig_dashboard.add_trace(
            go.Scatter(x=plot_data['Time'], y=plot_data['OH'], mode='lines', name='OH', line=dict(color='purple')),
            row=3, col=1
        )
    
    # Add Ala related elements plot
    for i, element in enumerate(ala_elements[:5]):  # Limit to first 5 elements to avoid overcrowding
        if plot_data[element].max() > 0:
            fig_dashboard.add_trace(
                go.Scatter(x=plot_data['Time'], y=plot_data[element], mode='lines', name=element),
                row=3, col=2
            )
    
    # Add Alb related elements plot
    for i, element in enumerate(alb_elements[:5]):  # Limit to first 5 elements to avoid overcrowding
        if plot_data[element].max() > 0:
            fig_dashboard.add_trace(
                go.Scatter(x=plot_data['Time'], y=plot_data[element], mode='lines', name=element),
                row=4, col=1
            )
    
    # Update layout
    fig_dashboard.update_layout(
        title_text='SPPARKS KMC Simulation Data Analysis Dashboard',
        height=900,
        width=1200,
        template='plotly_white',
        showlegend=True,
        legend=dict(orientation='h', yanchor='bottom', y=-0.2, xanchor='center', x=0.5)
    )
    
    # 添加交互式样式控件
    fig_dashboard = add_style_controls(fig_dashboard)
    
    # Save as HTML file
    dashboard_path = os.path.join(plots_dir, f'dashboard_{timestamp}.html')
    write_html_with_config(fig_dashboard, dashboard_path)
    return dashboard_path

def create_plot(plot_type, data, plots_dir, timestamp, style_config=None):
    """Create a single plot function, for parallel processing"""
    if plot_type == 'events':
        return ('events', create_events_plot(data, plots_dir, timestamp, style_config))
    elif plot_type == 'qcm' and 'QCM_ng_cm2' in data.columns:
        return ('qcm', create_qcm_plot(data, plots_dir, timestamp, style_config))
    elif plot_type == 'cpu':
        return ('cpu', create_cpu_plot(data, plots_dir, timestamp, style_config))
    elif plot_type == 'main':
        return ('main', create_main_elements_plot(data, plots_dir, timestamp, style_config))
    else:
        return (plot_type, None)

def create_interactive_plots_parallel(data, output_dir, plot_types=None, custom_elements=None, style_config=None):
    """Create interactive plots in parallel"""
    plots_dir = os.path.join(output_dir, "plots")
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
    
    # Create timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # If plot types not specified, default to all plots
    if plot_types is None:
        plot_types = ['events', 'cpu', 'qcm', 'main', 'dashboard']
    
    # Exclude dashboard, handle separately
    parallel_plot_types = [pt for pt in plot_types if pt != 'dashboard']
    
    output_files = {}
    
    if parallel_plot_types:
        print("Creating plots in parallel...")
        # Prepare parameters
        args = [(plot_type, data, plots_dir, timestamp, style_config) for plot_type in parallel_plot_types]
        
        # Determine number of processes, not exceeding CPU cores and task count
        num_processes = min(cpu_count(), len(args))
        
        # Parallel processing
        with Pool(processes=num_processes) as pool:
            results = list(tqdm(
                pool.starmap(create_plot, args),
                total=len(args),
                desc="Generating plots",
                unit="plots"
            ))
        
        # Process results
        for plot_type, path in results:
            if path:
                output_files[plot_type] = path
    
    # Handle dashboard separately (as it depends on other plots)
    if 'dashboard' in plot_types:
        print("Creating dashboard...")
        dashboard_path = create_dashboard(data, plots_dir, timestamp, style_config)
        output_files['dashboard'] = dashboard_path
    
    print(f"Interactive plots saved to directory: {plots_dir}")
    return output_files

def create_events_figure(data):
    """Create interactive chart of event count over time"""
    # Use complete data
    plot_data = data
    
    fig_events = go.Figure()
    fig_events.add_trace(go.Scatter(
        x=plot_data['Time'], 
        y=plot_data['Naccept'], 
        mode='lines+markers',
        name='Accepted Events',
        line=dict(color='blue', width=2),
        marker=dict(size=5)
    ))
    
    fig_events.update_layout(
        title='Accepted Events vs Time',
        xaxis_title='Time',
        yaxis_title='Event Count',
        template='plotly_white',
        hovermode='x unified',
        margin=dict(t=50, b=50, l=50, r=25),
        height=400
    )
    
    return fig_events

def create_qcm_figure(data):
    """Create interactive chart of QCM mass over time"""
    if 'QCM_ng_cm2' not in data.columns:
        return go.Figure()
    
    # Use complete data
    plot_data = data
    
    fig_qcm = go.Figure()
    fig_qcm.add_trace(go.Scatter(
        x=plot_data['Time'], 
        y=plot_data['QCM_ng_cm2'], 
        mode='lines+markers',
        name='QCM Mass',
        line=dict(color='purple', width=2),
        marker=dict(size=5)
    ))
    
    fig_qcm.update_layout(
        title='QCM Mass vs Time',
        xaxis_title='Time',
        yaxis_title='Mass (ng/cm²)',
        template='plotly_white',
        hovermode='x unified',
        margin=dict(t=50, b=50, l=50, r=25),
        height=400
    )
    
    return fig_qcm

def create_cpu_figure(data):
    """Create interactive chart of CPU time vs simulation time"""
    # Use complete data
    plot_data = data
    
    fig_cpu = go.Figure()
    fig_cpu.add_trace(go.Scatter(
        x=plot_data['Time'], 
        y=plot_data['CPU'], 
        mode='lines+markers',
        name='CPU Time',
        line=dict(color='red', width=2),
        marker=dict(size=5)
    ))
    
    fig_cpu.update_layout(
        title='CPU Time vs Simulation Time',
        xaxis_title='Simulation Time',
        yaxis_title='CPU Time (seconds)',
        template='plotly_white',
        hovermode='x unified',
        margin=dict(t=50, b=50, l=50, r=25),
        height=400
    )
    
    return fig_cpu

def create_main_elements_figure(data):
    """Create interactive chart of main element counts over time"""
    # Get all element columns (exclude non-element data columns)
    exclude_cols = ['Time', 'Naccept', 'Nreject', 'Nsweeps', 'CPU', 'events', 'QCM', 'QCM_ng_cm2']
    element_cols = [col for col in data.columns if col not in exclude_cols]
    
    # Main elements
    main_elements = ['VAC', 'O', 'OH', 'Ala', 'Alb', 'OHAlaX3', 'OHAlbX3']
    main_elements = [elem for elem in main_elements if elem in data.columns]
    
    # Use complete data
    plot_data = data
    
    fig_main = go.Figure()
    for element in main_elements:
        fig_main.add_trace(go.Scatter(
            x=plot_data['Time'], 
            y=plot_data[element], 
            mode='lines',
            name=element
        ))
    
    fig_main.update_layout(
        title='Main Element Counts vs Time',
        xaxis_title='Time',
        yaxis_title='Element Count',
        template='plotly_white',
        hovermode='x unified',
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1),
        margin=dict(t=50, b=50, l=50, r=25),
        height=400
    )
    
    return fig_main

def create_o_oh_figure(data):
    """Create interactive chart of O and OH element counts over time"""
    # Use complete data
    plot_data = data
    
    fig = go.Figure()
    
    if 'O' in plot_data.columns:
        fig.add_trace(go.Scatter(
            x=plot_data['Time'], 
            y=plot_data['O'], 
            mode='lines',
            name='O',
            line=dict(color='green')
        ))
    
    if 'OH' in plot_data.columns:
        fig.add_trace(go.Scatter(
            x=plot_data['Time'], 
            y=plot_data['OH'], 
            mode='lines',
            name='OH',
            line=dict(color='purple')
        ))
    
    fig.update_layout(
        title='O and OH Element Counts vs Time',
        xaxis_title='Time',
        yaxis_title='Element Count',
        template='plotly_white',
        hovermode='x unified',
        margin=dict(t=50, b=50, l=50, r=25),
        height=400
    )
    
    return fig

def create_ala_elements_figure(data):
    """Create interactive chart of Ala related elements over time"""
    # Get all element columns (exclude non-element data columns)
    exclude_cols = ['Time', 'Naccept', 'Nreject', 'Nsweeps', 'CPU', 'events', 'QCM', 'QCM_ng_cm2']
    element_cols = [col for col in data.columns if col not in exclude_cols]
    
    # Ala related elements
    ala_elements = [col for col in element_cols if 'Ala' in col]
    
    # Use complete data
    plot_data = data
    
    fig = go.Figure()
    
    for element in ala_elements[:10]:  # Limit to first 5 elements to avoid overcrowding
        if plot_data[element].max() > 0:
            fig.add_trace(go.Scatter(
                x=plot_data['Time'], 
                y=plot_data[element], 
                mode='lines',
                name=element
            ))
    
    fig.update_layout(
        title='Ala Related Elements vs Time',
        xaxis_title='Time',
        yaxis_title='Element Count',
        template='plotly_white',
        hovermode='x unified',
        margin=dict(t=50, b=50, l=50, r=25),
        height=400
    )
    
    return fig

def create_alb_elements_figure(data):
    """Create interactive chart of Alb related elements over time"""
    # Get all element columns (exclude non-element data columns)
    exclude_cols = ['Time', 'Naccept', 'Nreject', 'Nsweeps', 'CPU', 'events', 'QCM', 'QCM_ng_cm2']
    element_cols = [col for col in data.columns if col not in exclude_cols]
    
    # Alb related elements
    alb_elements = [col for col in element_cols if 'Alb' in col]
    
    # Use complete data
    plot_data = data
    
    fig = go.Figure()
    
    for element in alb_elements[:10]:  # Limit to first 5 elements to avoid overcrowding
        if plot_data[element].max() > 0:
            fig.add_trace(go.Scatter(
                x=plot_data['Time'], 
                y=plot_data[element], 
                mode='lines',
                name=element
            ))
    
    fig.update_layout(
        title='Alb Related Elements vs Time',
        xaxis_title='Time',
        yaxis_title='Element Count',
        template='plotly_white',
        hovermode='x unified',
        margin=dict(t=50, b=50, l=50, r=25),
        height=400
    )
    
    return fig

def export_data(data, output_dir, timestamp):
    """Export data as CSV and Excel formats"""
    # Create export directory
    export_dir = os.path.join(output_dir, 'exports')
    if not os.path.exists(export_dir):
        os.makedirs(export_dir)
    
    print("Exporting data...")
    
    # Export as CSV
    csv_path = os.path.join(export_dir, f'spparks_data_{timestamp}.csv')
    data.to_csv(csv_path, index=False)
    
    # Export as Excel
    excel_path = os.path.join(export_dir, f'spparks_data_{timestamp}.xlsx')
    data.to_excel(excel_path, index=False)
    
    # Export element statistics
    exclude_cols = ['Time', 'Naccept', 'Nreject', 'Nsweeps', 'CPU', 'events', 'QCM', 'QCM_ng_cm2']
    element_cols = [col for col in data.columns if col not in exclude_cols]
    
    # Calculate statistics for each element
    stats = {}
    
    # Use NumPy arrays for calculation, faster
    for col in tqdm(element_cols, desc="Calculating statistics", unit="elements"):
        # Convert to NumPy array and to tuple for caching
        data_array = data[col].to_numpy()
        data_tuple = tuple(data_array)
        
        # Calculate statistics
        element_stats = calculate_element_stats(col, data_tuple)
        if element_stats:
            stats[col] = element_stats
    
    # Export statistics as JSON
    stats_path = os.path.join(export_dir, f'element_stats_{timestamp}.json')
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=4)
    
    print(f"Data exported to directory: {export_dir}")
    print(f"CSV file: {csv_path}")
    print(f"Excel file: {excel_path}")
    print(f"Statistics: {stats_path}")
    
    return csv_path, excel_path, stats_path

def calculate_element_stats(element, data_tuple):
    """Calculate and cache element statistics"""
    # Convert tuple back to NumPy array
    data_array = np.array(data_tuple)
    
    if np.max(data_array) > 0:
        return {
            'min': float(np.min(data_array)),
            'max': float(np.max(data_array)),
            'mean': float(np.mean(data_array)),
            'std': float(np.std(data_array)),
            'change': float(np.max(data_array) - np.min(data_array))
        }
    return None

# 创建Dash应用函数
def create_dash_app(data, output_dir):
    """Create and run a Dash application for interactive data visualization"""
    # 创建Dash应用
    app = dash.Dash(__name__, 
                   external_stylesheets=[dbc.themes.BOOTSTRAP],
                   meta_tags=[{'name': 'viewport', 
                              'content': 'width=device-width, initial-scale=1.0'}])
    
    # 设置应用标题
    app.title = "SPPARKS KMC Simulation Data Analysis"
    
    # 获取元素列
    exclude_cols = ['Time', 'Naccept', 'Nreject', 'Nsweeps', 'CPU', 'events', 'QCM', 'QCM_ng_cm2']
    element_cols = [col for col in data.columns if col not in exclude_cols]
    
    # 创建应用布局
    app.layout = dbc.Container([
        dbc.Row([
            dbc.Col([
                html.H1("SPPARKS KMC Simulation Data Analysis", 
                        className="text-center my-4"),
                html.Hr()
            ], width=12)
        ]),
        
        # 控制面板
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("控制面板"),
                    dbc.CardBody([
                        html.H5("时间范围选择"),
                        dcc.RangeSlider(
                            id='time-slider',
                            min=data['Time'].min(),
                            max=data['Time'].max(),
                            value=[data['Time'].min(), data['Time'].max()],
                            marks={
                                str(round(data['Time'].min(), 2)): str(round(data['Time'].min(), 2)),
                                str(round(data['Time'].max(), 2)): str(round(data['Time'].max(), 2))
                            },
                            tooltip={"placement": "bottom", "always_visible": True}
                        ),
                        
                        html.Div(className="my-3"),
                        
                        html.H5("Y轴比例"),
                        dbc.ButtonGroup([
                            dbc.Button("线性", id="linear-scale", color="primary", outline=True, active=True, className="me-1"),
                            dbc.Button("对数", id="log-scale", color="primary", outline=True, className="me-1")
                        ], className="mb-3"),
                        
                        html.Div(className="my-3"),
                        
                        html.H5("元素选择"),
                        dcc.Dropdown(
                            id='element-dropdown',
                            options=[{'label': elem, 'value': elem} for elem in element_cols],
                            value=element_cols[:5] if len(element_cols) > 5 else element_cols,
                            multi=True
                        ),
                        
                        html.Div(className="my-3"),
                        
                        dbc.Button("导出数据", id="export-button", color="success", className="mt-2")
                    ])
                ], className="mb-4")
            ], width=12)
        ]),
        
        # 图表行 - 第一行
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("接受事件 vs 时间"),
                    dbc.CardBody([
                        dcc.Graph(id='events-graph', figure=create_events_figure(data))
                    ])
                ])
            ], width=6),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("CPU时间 vs 模拟时间"),
                    dbc.CardBody([
                        dcc.Graph(id='cpu-graph', figure=create_cpu_figure(data))
                    ])
                ])
            ], width=6)
        ], className="mb-4"),
        
        # 图表行 - 第二行
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("QCM质量 vs 时间"),
                    dbc.CardBody([
                        dcc.Graph(id='qcm-graph', figure=create_qcm_figure(data))
                    ])
                ])
            ], width=6),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("主要元素计数 vs 时间"),
                    dbc.CardBody([
                        dcc.Graph(id='main-elements-graph', figure=create_main_elements_figure(data))
                    ])
                ])
            ], width=6)
        ], className="mb-4"),
        
        # 图表行 - 第三行
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("O和OH元素计数 vs 时间"),
                    dbc.CardBody([
                        dcc.Graph(id='o-oh-graph', figure=create_o_oh_figure(data))
                    ])
                ])
            ], width=6),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Ala相关元素"),
                    dbc.CardBody([
                        dcc.Graph(id='ala-elements-graph', figure=create_ala_elements_figure(data))
                    ])
                ])
            ], width=6)
        ], className="mb-4"),
        
        # 图表行 - 第四行
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Alb相关元素"),
                    dbc.CardBody([
                        dcc.Graph(id='alb-elements-graph', figure=create_alb_elements_figure(data))
                    ])
                ])
            ], width=6),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("自定义元素图表"),
                    dbc.CardBody([
                        dcc.Graph(id='custom-elements-graph')
                    ])
                ])
            ], width=6)
        ], className="mb-4"),
        
        # 底部信息
        dbc.Row([
            dbc.Col([
                html.Hr(),
                html.P("SPPARKS KMC Simulation Data Analysis Dashboard", className="text-center"),
                html.P(f"Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", className="text-center")
            ], width=12)
        ])
    ], fluid=True)
    
    # 回调函数 - 时间范围选择器
    @app.callback(
        [Output('events-graph', 'figure'),
         Output('cpu-graph', 'figure'),
         Output('qcm-graph', 'figure'),
         Output('main-elements-graph', 'figure'),
         Output('o-oh-graph', 'figure'),
         Output('ala-elements-graph', 'figure'),
         Output('alb-elements-graph', 'figure')],
        [Input('time-slider', 'value'),
         Input('linear-scale', 'n_clicks'),
         Input('log-scale', 'n_clicks')]
    )
    def update_time_range(time_range, linear_clicks, log_clicks):
        # 确定当前上下文
        ctx = dash.callback_context
        button_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None
        
        # 确定Y轴类型
        y_type = 'linear'
        if button_id == 'log-scale':
            y_type = 'log'
        
        # 过滤数据
        filtered_data = data[(data['Time'] >= time_range[0]) & (data['Time'] <= time_range[1])]
        
        # 创建图表
        events_fig = create_events_figure(filtered_data)
        cpu_fig = create_cpu_figure(filtered_data)
        qcm_fig = create_qcm_figure(filtered_data)
        main_fig = create_main_elements_figure(filtered_data)
        o_oh_fig = create_o_oh_figure(filtered_data)
        ala_fig = create_ala_elements_figure(filtered_data)
        alb_fig = create_alb_elements_figure(filtered_data)
        
        # 更新Y轴类型
        events_fig.update_layout(yaxis_type=y_type)
        cpu_fig.update_layout(yaxis_type=y_type)
        qcm_fig.update_layout(yaxis_type=y_type)
        main_fig.update_layout(yaxis_type=y_type)
        o_oh_fig.update_layout(yaxis_type=y_type)
        ala_fig.update_layout(yaxis_type=y_type)
        alb_fig.update_layout(yaxis_type=y_type)
        
        return events_fig, cpu_fig, qcm_fig, main_fig, o_oh_fig, ala_fig, alb_fig
    
    # 回调函数 - 自定义元素图表
    @app.callback(
        Output('custom-elements-graph', 'figure'),
        [Input('element-dropdown', 'value'),
         Input('time-slider', 'value'),
         Input('linear-scale', 'n_clicks'),
         Input('log-scale', 'n_clicks')]
    )
    def update_custom_elements(selected_elements, time_range, linear_clicks, log_clicks):
        # 确定当前上下文
        ctx = dash.callback_context
        button_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None
        
        # 确定Y轴类型
        y_type = 'linear'
        if button_id == 'log-scale':
            y_type = 'log'
        
        # 过滤数据
        filtered_data = data[(data['Time'] >= time_range[0]) & (data['Time'] <= time_range[1])]
        
        # 创建图表
        fig = go.Figure()
        
        if selected_elements:
            for element in selected_elements:
                if element in filtered_data.columns:
                    fig.add_trace(go.Scatter(
                        x=filtered_data['Time'], 
                        y=filtered_data[element], 
                        mode='lines',
                        name=element
                    ))
        
        fig.update_layout(
            title='自定义元素计数 vs 时间',
            xaxis_title='时间',
            yaxis_title='元素计数',
            yaxis_type=y_type,
            template='plotly_white',
            hovermode='x unified',
            margin=dict(t=50, b=50, l=50, r=25),
            height=400
        )
        
        return fig
    
    # 回调函数 - 导出数据
    @app.callback(
        Output('export-button', 'children'),
        [Input('export-button', 'n_clicks')],
        [State('time-slider', 'value')]
    )
    def export_filtered_data(n_clicks, time_range):
        if n_clicks is None:
            return "导出数据"
        
        # 过滤数据
        filtered_data = data[(data['Time'] >= time_range[0]) & (data['Time'] <= time_range[1])]
        
        # 创建时间戳
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 导出数据
        export_data(filtered_data, output_dir, timestamp)
        
        return "数据已导出!"
    
    return app

# 添加信号处理和清理函数
def cleanup_server(server=None):
    """清理函数，确保服务器正常关闭"""
    print("正在关闭服务器...")
    if server:
        try:
            server.stop()
            print("服务器已成功关闭")
        except:
            print("服务器关闭过程中出现错误")

def signal_handler(sig, frame, server=None):
    """处理信号的函数"""
    print(f"接收到信号 {sig}，正在关闭...")
    cleanup_server(server)
    sys.exit(0)

def main():
    """Main function"""
    start_time = time.time()
    
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='SPPARKS KMC Simulation Data Visualization Tool')
    
    # 基本参数
    parser.add_argument('--log-file', type=str, 
                        default='/Users/kazami/spparks/spparks-6Sep23/SPPARKS-for-Al2O3_TMA_H2O/AB_Surface/120/log.spparks',
                        help='SPPARKS log file path')
    parser.add_argument('--model-file', type=str, 
                        default='/Users/kazami/spparks/spparks-6Sep23/SPPARKS-for-Al2O3_TMA_H2O/AB_Surface/120/data_10X10.ald',
                        help='SPPARKS model file path')
    parser.add_argument('--output-dir', type=str, 
                        default='/Users/kazami/spparks/spparks-6Sep23/SPPARKS-for-Al2O3_TMA_H2O/AB_Surface/120/plots',
                        help='Output directory for plots and data')
    parser.add_argument('--no-dash', action='store_true',
                        help='只生成图表文件，不启动Dash服务器')
    
    # 解析参数
    args = parser.parse_args()
    
    # 解析模型尺寸
    x_dim, y_dim, area_cm2 = parse_model_dimensions(args.model_file)
    
    # 解析日志文件
    data = parse_log_file(args.log_file)
    
    # 转换QCM单位
    data = convert_qcm_to_ng_cm2(data, area_cm2)
    
    # 创建输出目录
    output_dir, plots_dir, exports_dir, timestamp = create_output_directory(args.output_dir)
    
    # 导出数据
    export_data(data, output_dir, timestamp)
    
    # 创建交互式图表
    print("创建交互式图表...")
    plot_files = create_interactive_plots_parallel(data, output_dir)
    
    print(f"总处理时间: {time.time() - start_time:.2f} 秒")
    print(f"图表已保存到目录: {plots_dir}")
    print(f"仪表板文件: {plot_files.get('dashboard', '未创建')}")
    
    # 如果指定了--no-dash参数，则不启动Dash服务器
    if args.no_dash:
        print("已指定--no-dash参数，不启动Dash服务器")
        return
    
    # 创建并运行Dash应用
    app = create_dash_app(data, output_dir)
    
    # 尝试不同的端口，直到找到一个可用的
    port = 8050
    max_port = 8100  # 最大尝试端口号
    server = None
    
    while port <= max_port:
        try:
            print(f"尝试在端口 {port} 上启动Dash服务器...")
            print(f"启动后请在浏览器中打开 http://127.0.0.1:{port}/ 查看仪表板")
            print("按 Ctrl+C 停止服务器")
            
            # 注册信号处理器
            server = app.server
            signal.signal(signal.SIGINT, lambda sig, frame: signal_handler(sig, frame, server))
            signal.signal(signal.SIGTERM, lambda sig, frame: signal_handler(sig, frame, server))
            
            # 注册退出处理器
            atexit.register(lambda: cleanup_server(server))
            
            # 启动服务器
            app.run(debug=False, port=port)
            break  # 如果成功启动，跳出循环
        except OSError as e:
            if "Address already in use" in str(e):
                print(f"端口 {port} 已被占用，尝试下一个端口...")
                port += 1
            else:
                # 其他类型的OSError，直接抛出
                raise e
    
    if port > max_port:
        print("无法找到可用端口。请手动终止正在运行的Dash应用程序，然后重试。")
        print("您可以使用以下命令查找并终止占用端口的进程：")
        print("  macOS/Linux: lsof -i :8050 | grep LISTEN")
        print("  然后使用: kill <PID> 终止进程")

if __name__ == "__main__":
    main()