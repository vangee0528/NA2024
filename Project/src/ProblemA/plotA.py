import numpy as np
import matplotlib.pyplot as plt
import os
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

def plot_polynomials(intervals, polynomials, output_path, N):
    plt.figure()
    
    # 绘制原函数
    x = np.linspace(-1, 1, 1000)
    y = original_function(x)
    plt.plot(x, y, 'k--', label='Original function')
    
    # 绘制拟合多项式
    for interval, coeffs in zip(intervals, polynomials):
        x = np.linspace(interval[0], interval[1], 400)
        y = np.polyval(coeffs[::-1], x)
        plt.plot(x, y, label=f'Interval {interval}')
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Piecewise Polynomial Plot with N = {}'.format(N))
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

def log_message(message):
    log_file = 'logs/log.txt'
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f'{timestamp} : {message}\n')

def main():
    input_dir = './output/problemA'
    output_dir = './figure/problemA'
    N_values = [6, 11, 21, 41, 81]
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for N in N_values:
        input_file = os.path.join(input_dir, f'N_{N}.txt')
        output_file = os.path.join(output_dir, f'N_{N}.png')
        
        intervals, polynomials = read_data(input_file)
        plot_polynomials(intervals, polynomials, output_file, N)
        log_message(f'Successfully saved {output_file}')

if __name__ == '__main__':
    main()