import numpy as np
import matplotlib.pyplot as plt
import os

# 曲线r(t) = (1/(1+25t^2), sin(t))  
def r1(t):
    x = 1.0 / (1 + 25 * t * t)
    y = np.sin(t)
    return x, y



def read_polynomial(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    
    polynomials = {'x': [], 'y': []}
    current = None
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == 'x(t) =':
            current = 'x'
            i += 1
        elif line == 'y(t) =':
            current = 'y'
            i += 1
        elif line == 'z(t) =':
            current = 'z'
            i += 1
            polynomials[current] = []
        elif current:
            interval = tuple(map(float, line.split(',')))
            i += 1
            coefficients = list(map(float, lines[i].strip().split()))
            polynomials[current].append((interval, coefficients))
            i += 1
        else:
            i += 1
    
    return polynomials

def evaluate_polynomial(coefficients, t):
    return sum(c * t**i for i, c in enumerate(coefficients))

def plot_r(k):
    if k == 1:  
        t = np.linspace(-np.pi, np.pi, 1000)
        x, y = r1(t)
        plt.plot(x, y, linestyle='--', label='r1(t)')

def plot_polynomial(k, polynomials, label, output_file):
    t_values = []
    x_values = []
    y_values = []     
    
    for interval, coefficients in polynomials['x']:
        t = np.linspace(interval[0], interval[1], 100)
        x = [evaluate_polynomial(coefficients, ti) for ti in t]
        t_values.extend(t)
        x_values.extend(x)
    
    for interval, coefficients in polynomials['y']:
        t = np.linspace(interval[0], interval[1], 100)
        y = [evaluate_polynomial(coefficients, ti) for ti in t]
        y_values.extend(y)
           
    if k == 1 or k == 2:
        plt.figure()
        plt.plot(x_values, y_values, label=label)
        plot_r(1)
        plt.xlabel('x(t)')
        plt.ylabel('y(t)')
        plt.title(f'Parametric Curve {label}')
        plt.legend()
        plt.savefig(output_file)
        log_message("Successfully saved " + output_file)
        plt.close()

from datetime import datetime
def log_message(message):
    log_file = 'logs/log.txt'
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f'{timestamp} : {message}\n')

def main():
    # 读取并绘制之前的曲线
    files = [
        ('output/jsontest/output_bspline.txt', 1),
        ('output/jsontest/output_ppspline.txt', 2)
    ]
    output_dir = 'figure/jsontest'

    os.makedirs(output_dir, exist_ok=True)
    
    for file, curve_num in files:
        if os.path.exists(file):
            polynomials = read_polynomial(file)
            label = os.path.basename(file).replace('.txt', '')
            output_file = os.path.join(output_dir, f'{label}.png')
            plot_polynomial(curve_num, polynomials, label, output_file)
        else:
            log_message(f"File {file} does not exist.")

if __name__ == "__main__":
    main()