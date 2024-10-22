import numpy as np
import matplotlib.pyplot as plt

# 目标函数 f(x) = 1 / (1 + x^2)
def f(x):
    return 1 / (1 + x**2)

# 读取插值数据
def read_interpolation_data(filename):
    data = {}
    with open(filename, 'r', encoding='utf-8') as file:  
        lines = file.readlines()
        
        n = None
        current_data = []
        
        for line in lines:
            if line.startswith('n ='):
                if n is not None:
                    data[n] = np.array(current_data)
                n = int(line.split('=')[1].strip())
                current_data = []
            elif line.strip() == "====":
                continue
            else:
                x, y = map(float, line.split())
                current_data.append([x, y])
        
        if n is not None:
            data[n] = np.array(current_data)
    
    return data



x_values = np.linspace(-5, 5, 400)
f_values = f(x_values)

interpolation_data = read_interpolation_data("data/ProblemB.txt")

plt.figure(figsize=(10, 6))

plt.plot(x_values, f_values, label='Original function $f(x) = \\frac{1}{1 + x^2}$', color='black', linewidth=2)

colors = ['red', 'blue', 'green', 'orange']
for i, n in enumerate(interpolation_data):
    interpolated_vals = interpolation_data[n]
    plt.plot(interpolated_vals[:, 0], interpolated_vals[:, 1], label=f'Newton interpolation (n={n})', color=colors[i])

plt.title('ProblemB : Runge Phenomenon with Newton Interpolation')
plt.xlabel('$x$')
plt.ylabel('$f(x)$ and Interpolated Values')
plt.legend()
plt.grid(True)

plt.savefig('figure/ProblemB_figure.png')
plt.show()
