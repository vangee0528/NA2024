import numpy as np
import matplotlib.pyplot as plt

def draw_plot(lines):
    # 解析方程
    equations = []
    for i in range(0, len(lines), 3):
        if i + 3 < len(lines):  # 确保有四行
            # 第一行为 x 取值，转化为浮点数
            x = float(lines[i].strip().split('=')[1].strip())
            # 第二行为 x(t) 的表达式
            x_eq = lines[i + 1].strip().split('=')[1].strip()
            # 第三行为 y(t) 的表达式

            y_eq = lines[i + 2].strip().split('=')[1].strip()
            equations.append((x, x_eq, y_eq))
        else:
            print(f"Warning: Missing y equation for x equation at line {i + 1}")

    # 绘制每组方程
    t = np.linspace(0, 1, 100)
    plt.figure(figsize=(10, 8))

    for x_allow,x_eq, y_eq in equations:
        # 使用 eval 来计算 x(t) 和 y(t)
        x = eval(x_eq)
        y = eval(y_eq)
        
        # 筛选出 x(t) 在允许范围内 [x_allow-0.1, x_allow+0.1] 的索引
        indices = [i for i in range(len(x)) if x_allow - 0.01 <= x[i] <= x_allow + 0.01]
        
        # 使用筛选出的索引来筛选 x 和 y
        x = [x[i] for i in indices]
        y = [y[i] for i in indices]
        

        plt.plot(x, y)

    #绘制心形曲线 x^(2)+(((3)/(2)) y-sqrt(abs(x)))^(2)=3
    for i in np.linspace(-1.732, 1.732, 1000):
        x = i
        y1 = 2.0 / 3.0 * (np.sqrt(3 - x * x) + np.sqrt(abs(x)))
        y2 = 2.0 / 3.0 * (-np.sqrt(3 - x * x) + np.sqrt(abs(x)))
        plt.plot(x, y1, 'ro')
        plt.plot(x, y2, 'ro')


    plt.title('Parametric Curves from Problem F')
    plt.xlabel('x(t)')
    plt.ylabel('y(t)')
    plt.grid()
    plt.axis('equal')
    plt.show()


with open('data/ProblemF_1.txt', 'r') as file:
    lines = file.readlines()
draw_plot(lines)

with open('data/ProblemF_2.txt', 'r') as file:
    lines = file.readlines()   
draw_plot(lines)

with open('data/ProblemF_3.txt', 'r') as file:
    lines = file.readlines()
draw_plot(lines)