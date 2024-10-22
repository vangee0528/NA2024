import matplotlib.pyplot as plt
import numpy as np

def generate_expression(polynomial_expression):

    e_index = polynomial_expression.find('e')
    if e_index != -1:
        i = e_index + 1
        while i < len(polynomial_expression) and polynomial_expression[i] != 'x':
            i += 1
        e_number = polynomial_expression[e_index + 1:i]
        polynomial_expression = polynomial_expression.replace(f'e{e_number}', f'*10**({int(e_number)})')

    polynomial_expression = polynomial_expression.replace('x', '*x')

    def f(x):
        return eval(polynomial_expression)
    return f

def draw_plot(polynomial_expression_sp1,polynomial_expression_sp2,points1,points2,x_text,y_text,title,end):

    f = generate_expression(polynomial_expression_sp1)
    g = generate_expression(polynomial_expression_sp2)
    
    # 绘制多项式曲线
    x_values = np.linspace(0, end, 400)
    f1_values = f(x_values)
    f2_values = g(x_values)

    plt.figure(figsize=(10, 6))
    plt.plot(x_values, f1_values, label='Polynomial expression: ' + polynomial_expression_sp1, color='black', linewidth=2)
    plt.plot(x_values, f2_values, label='Polynomial expression: ' + polynomial_expression_sp2, color='blue', linewidth=2)
    plt.title(f'Problem E : {title}')
    plt.xlabel(x_text)
    plt.ylabel(y_text)
    plt.legend(['sp1','sp2'])


    for i in range(len(points1)):
        plt.scatter(points1[i][0], points1[i][1], color='red', s=50)
    for i in range(len(points2)):
        plt.scatter(points2[i][0], points2[i][1], color='blue', s=50)


    plt.savefig(f'figure/ProblemE_figure_{title}.png')

    plt.show()
    
with open('data/ProblemE.txt', 'r') as file:
    lines = file.readlines()


points1 = []
points2 = []
for i in range(5, 12):
    x, y1, y2 = map(float, lines[i].split())
    points1.append([x, y1])
    points2.append([x, y2])

polynomial_expression_sp1 = lines[1].strip()
polynomial_expression_sp2 = lines[3].strip()
draw_plot(polynomial_expression_sp1,polynomial_expression_sp2,points1,points2,'Day','Gram','Known',28)
draw_plot(polynomial_expression_sp1,polynomial_expression_sp2,points1,points2,'Day','Gram','Predicted',45)










