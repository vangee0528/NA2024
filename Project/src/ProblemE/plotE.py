import numpy as np
import matplotlib.pyplot as plt
import os

# 曲线r1(t) = (sqrt(3)cos(t), 2/3(sqrt(sqrt(3)|cos(t)|) + sqrt(3)sin(t)))
def r1(t):
    x = np.sqrt(3) * np.cos(t)
    y = 2.0 / 3.0 * (np.sqrt(np.sqrt(3) * np.abs(np.cos(t))) + np.sqrt(3) * np.sin(t))
    return x, y

# 曲线r2(t) = (sin(t) + tcos(t), cos(t) - tsin(t))
def r2(t):
    x = np.sin(t) + t * np.cos(t)
    y = np.cos(t) - t * np.sin(t)
    return x, y

# 曲线r3(t) = [sin(cos(t))*cos(sin(t)), sin(cos(t))*sin(sin(t)), cos(cos(t))]
def r3(t):
    x = np.sin(np.cos(t)) * np.cos(np.sin(t))
    y = np.sin(np.cos(t)) * np.sin(np.sin(t))
    z = np.cos(np.cos(t))
    return x, y, z

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
    elif k == 2:
        t = np.linspace(0, 6 * np.pi, 1000)
        x, y = r2(t)
        plt.plot(x, y, linestyle='--', label='r2(t)')

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
           
    if 'z' in polynomials:
        z_values = []
        for interval, coefficients in polynomials['z']:
            t = np.linspace(interval[0], interval[1], 100)
            z = [evaluate_polynomial(coefficients, ti) for ti in t]
            z_values.extend(z)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x_values, y_values, z_values, label=label)
        t = np.linspace(0, 2 * np.pi, 1000)
        x, y, z = r3(t)
        ax.plot(x, y, z, linestyle='--', label='r3(t)')
        ax.set_xlabel('x(t)')
        ax.set_ylabel('y(t)')
        ax.set_zlabel('z(t)')
        plt.title(f'Parametric Curve {label}')
        plt.legend()
        plt.savefig(output_file)
        print("Successfully saved", output_file)
        plt.close()

    else:
        plt.figure()
        plt.plot(x_values, y_values, label=label)
        plot_r(k)
        plt.xlabel('x(t)')
        plt.ylabel('y(t)')
        plt.title(f'Parametric Curve {label}')
        plt.legend()
        plt.savefig(output_file)
        print("Successfully saved", output_file)
        plt.close()

def main():
    curve_name = [1, 2, 3]
    N_name = ['N10', 'N40', 'N160']
    type_name = ['unit', 'chord']
    prefix = 'output/problemE/r'

    files = [[f'{prefix}{curve}_{type}_{N}.txt', curve] for curve in curve_name for N in N_name for type in type_name]
    output_dir = 'figure/problemE'
    
    os.makedirs(output_dir, exist_ok=True)
    
    for file, curve_num in files:
        if os.path.exists(file):
            polynomials = read_polynomial(file)
            label = os.path.basename(file).replace('.txt', '')
            output_file = os.path.join(output_dir, f'{label}.png')
            plot_polynomial(curve_num, polynomials, label, output_file)
        else:
            print(f"File {file} does not exist.")

if __name__ == "__main__":
    main()