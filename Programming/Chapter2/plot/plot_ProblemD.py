import matplotlib.pyplot as plt
import numpy as np

def generate_expression(polynomial_expression):
    # 将表达式中的 e 替换成 *10**
    e_index = polynomial_expression.find('e')
    if e_index != -1:
        i = e_index + 1
        while i < len(polynomial_expression) and polynomial_expression[i] != 'x':
            i += 1
        e_number = polynomial_expression[e_index + 1:i]
        polynomial_expression = polynomial_expression.replace(f'e{e_number}', f'*10**({int(e_number)})')

    # 将表达式中的 x 替换成 *x
    polynomial_expression = polynomial_expression.replace('x', '*x')

    def f(x):
        return eval(polynomial_expression)
    return f

def draw_plot(polynomial_expression,points,x_text,y_text,title,y=-1):

    f = generate_expression(polynomial_expression)

    x_values = np.linspace(0, 13, 400)
    f_values = f(x_values)

    plt.figure(figsize=(10, 6))
    plt.plot(x_values, f_values, label='Polynomial expression: ' + polynomial_expression, color='black', linewidth=2)
    plt.title(f'Problem D : {title}')
    plt.xlabel(x_text)
    plt.ylabel(y_text)

    for i in range(len(points)):
        plt.scatter(points[i][0], points[i][1], color='red', s=50)

    if y!=-1:
        plt.plot(x_values, y*np.ones_like(x_values), label='y='+str(y), color='red', linewidth=2)

    plt.savefig(f'figure/ProblemD_figure_{title}.png')

    plt.show()
    
with open('data/ProblemD.txt', 'r') as file:
    lines = file.readlines()

points1 = []
points2 = []
for i in range(5, 10):
    x, y1, y2 = map(float, lines[i].split())
    points1.append([x, y1])
    points2.append([x, y2])


polynomial_expression = lines[1].strip()
draw_plot(polynomial_expression,points1,'Time/s','Distance/feet','Distance-Time Curve')

polynomial_expression = lines[3].strip()
draw_plot(polynomial_expression,points2,'Time/s','Velocity/feet/s','Velocity-Time Curve',81)







