import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from datetime import datetime

def read_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    intervals = []
    polynomials = []
    
    for i in range(0, len(lines), 2):
        interval = list(map(float, lines[i].strip().split(',')))
        polynomial = list(map(float, lines[i+1].strip().split()))
        intervals.append(interval)
        polynomials.append(polynomial)
    
    return intervals, polynomials

def original_function(x):
    return 1 / (1 + 25 * x**2)

def plot_polynomials(intervals, polynomials, output_path):
    plt.figure()
    
    # 绘制拟合多项式
    for interval, coeffs in zip(intervals, polynomials):
        x = np.linspace(interval[0], interval[1], 400)
        y = np.polyval(coeffs[::-1], x)
        plt.plot(x, y, label=f'Interval {interval}')
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
 
    plt.savefig(output_path)
    plt.close()

def log_message(message):
    log_file = 'logs/log.txt'
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f'{timestamp} : {message}\n')

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 src/Check/plot.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    if not os.path.exists(input_file):
        print(f"File {input_file} does not exist.")
        sys.exit(1)
    
    output_dir = './figure/check'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    file_name = os.path.basename(input_file)
    output_file = os.path.join(output_dir, file_name.replace('.txt', '.png'))
    
    intervals, polynomials = read_data(input_file)
    plot_polynomials(intervals, polynomials, output_file)
    log_message(f'Successfully saved {output_file}')

if __name__ == '__main__':
    main()