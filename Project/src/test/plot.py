import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime

def read_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    intervals = []
    polynomials = []
    current_intervals = []
    current_polynomials = []
    
    for line in lines:
        if line.strip() == '====':
            if current_intervals and current_polynomials:
                intervals.append(current_intervals)
                polynomials.append(current_polynomials)
                current_intervals = []
                current_polynomials = []
        else:
            if ',' in line:
                interval = list(map(float, line.strip().split(',')))
                current_intervals.append(interval)
            else:
                polynomial = list(map(float, line.strip().split()))
                current_polynomials.append(polynomial)
    
    if current_intervals and current_polynomials:
        intervals.append(current_intervals)
        polynomials.append(current_polynomials)
    
    return intervals, polynomials

def plot_polynomials(intervals, polynomials, output_path):
    plt.figure()

    # 绘制拟合多项式
    for idx, (interval_set, polynomial_set) in enumerate(zip(intervals, polynomials)):
        for interval, coeffs in zip(interval_set, polynomial_set):
            x = np.linspace(interval[0], interval[1], 400)
            y = np.polyval(coeffs[::-1], x)
            plt.plot(x, y, label=f'Function {idx+1}, Interval {interval}')
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Piecewise Polynomial Plot')
    plt.grid(True)
    # plt.legend()
    plt.savefig(output_path)
    plt.close()

def log_message(message):
    log_file = 'logs/log.txt'
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f'{timestamp} : {message}\n')

def main():
    input_file = './output/test/func.txt'
    output_dir = './figure/test'
    output_file = os.path.join(output_dir, 'func_plot.png')
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    intervals, polynomials = read_data(input_file)
    plot_polynomials(intervals, polynomials, output_file)
    log_message(f'Successfully saved {output_file}')

if __name__ == '__main__':
    main()