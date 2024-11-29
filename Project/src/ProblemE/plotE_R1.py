import numpy as np
import matplotlib.pyplot as plt
import os

def r1(t):
    x = np.sqrt(3) * np.cos(t)
    y = 2.0 / 3.0 * (np.sqrt(np.sqrt(3) * np.abs(np.cos(t))) + np.sqrt(3) * np.sin(t))
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

def plot_polynomial(polynomials, label, output_file):
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
    
    plt.figure()
    plt.plot(x_values, y_values, label=label)
    
    # 绘制目标函数 r1(t)
    t = np.linspace(0, 2 * np.pi, 1000)
    x_r1, y_r1 = r1(t)
    plt.plot(x_r1, y_r1, label='r1(t)', linestyle='--')
    
    plt.xlabel('x(t)')
    plt.ylabel('y(t)')
    plt.title(f'Parametric Curve {label}')
    plt.legend()
    plt.savefig(output_file)
    print("Successfully saved", output_file)
    plt.close()

def main():
    files = ['output/problemE/r1_unit_N10.txt', 'output/problemE/r1_unit_N40.txt', 'output/problemE/r1_unit_N160.txt', 'output/problemE/r1_chord_N10.txt', 'output/problemE/r1_chord_N40.txt', 'output/problemE/r1_chord_N160.txt']
    output_dir = 'figure/problemE'
    os.makedirs(output_dir, exist_ok=True)
    
    for file in files:
        if os.path.exists(file):
            polynomials = read_polynomial(file)
            label = os.path.basename(file).replace('.txt', '')
            output_file = os.path.join(output_dir, f'{label}.png')
            plot_polynomial(polynomials, label, output_file)
        else:
            print(f"File {file} does not exist.")

if __name__ == "__main__":
    main()