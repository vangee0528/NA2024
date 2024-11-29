import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def read_polynomial(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    
    polynomials = {'x': [], 'y': []}
    current = None
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == 'x(t) =' or line == 'y(t) =' or line == 'z(t) =':
            current = line.split('(')[0]
            i += 1
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

def plot_polynomial(file, polynomials, color):
    x_values = []
    y_values = []
    for interval, coefficients in polynomials['x']:
        t = np.linspace(interval[0], interval[1], 100)
        x = [evaluate_polynomial(coefficients, ti) for ti in t]
        x_values.extend(x)
    for interval, coefficients in polynomials['y']:
        t = np.linspace(interval[0], interval[1], 100)
        y = [evaluate_polynomial(coefficients, ti) for ti in t]
        y_values.extend(y)
    plt.plot(x_values, y_values, color=color, label=file)

from datetime import datetime
def log_message(message):
    log_file = 'logs/log.txt'
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f'{timestamp} : {message}\n')

def main():
    if len(sys.argv) != 3:
        log_message("Usage: python3 compare.py <file1> <file2>")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]

    if not os.path.exists(file1):
        log_message(f"File {file1} does not exist.")
        sys.exit(1)
    
    if not os.path.exists(file2):
        log_message(f"File {file2} does not exist.")
        sys.exit(1)
    
    polynomials1 = read_polynomial(file1)
    polynomials2 = read_polynomial(file2)

    plt.figure()
    plot_polynomial(file1, polynomials1, color='blue')
    plot_polynomial(file2, polynomials2, color='red')

    plt.xlabel('x(t)')
    plt.ylabel('y(t)')
    plt.title('Compared Parametric Curves')
    plt.legend()
    output_dir = 'figure/check'
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'P4_compare.png')
    plt.savefig(output_file)
    log_message(f'Successfully saved {output_file}')
    plt.close()

if __name__ == "__main__":
    main()