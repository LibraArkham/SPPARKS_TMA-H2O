import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd 
import argparse
import os

#read file
def read_lines(file_path, str):
  str_line = 0
  with open (file_path, 'r') as f:
       for line in f:
           str_line += 1
           if line.startswith(str):
               break
  return str_line

#save data as a table
def save_table(file_path, start_line, end_line,filename):
    with open (file_path, 'r') as f:
        lines = f.readlines()     
    with open (filename, 'w') as d:
        d.writelines(lines[start_line:end_line])

#draw and save a figure  
def draw_fig(file_path,T,target):
    data = pd.read_csv(file_path, delimiter='\s+')
    y1 = data[target].to_numpy()
    x1 = data['Time'].to_numpy()
    plt.plot(x1,y1*0.0332*0.75,linewidth=3.0)
    plt.title(str(T)+'K',fontweight='bold',size = 22)
    plt.xlabel("Time(s)",fontweight='bold',size = 20)
    plt.ylabel(target+"(ng/cm^2)",fontweight='bold',size = 20)
    plt.xticks(fontweight='bold',fontsize = 15)
    plt.yticks(fontweight='bold',fontsize = 15)
    plt.savefig(target+'.png')
    plt.show()
    
def draw_log(log,temperature,target):
    start_line = read_lines(log, "Memory")
    end_line = read_lines(log, "Loop")
    save_table(log, start_line, end_line-1,"log.data")
    draw_fig('log.data',temperature,target)

#analyse log.spparks
def anl_log(log):  
    pre_line = read_lines(log, '#precursor:')
    time_line = read_lines(log, 'pulse_time')
    temperature_line = read_lines(log, 'temperature')
    with open (log,'r') as l:
        lines = l.readlines()
    return lines[pre_line,time_line,(time_line+1),temperature_line]

#analyse data.ald
def anl_data(data):
    start_line = read_lines(data,'simulation box')
    with open (data,'r') as d:
        lines = d.readlines()
    return lines[start_line:(start_line+6)]

#def anl_dump(dump):
    
def anl(log,data,dump,ana):
    temperature = float(input("input temperature(K):"))
    target = str(input('QCM/events/species/reactions:'))
    #log_lines = anl_log(log)
    #data_lines = anl_data(data)
    #with open ('analyse.txt','w') as ana:
    #    ana.writelines(log_lines+'\n'+data_lines)
    draw_log(log,temperature,target)

        

anl("log.spparks","data.ald","dump.ald","analyse.txt")