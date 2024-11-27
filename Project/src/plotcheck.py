import sys
import matplotlib.pyplot as plt
import numpy as np

def plot_checkpoint_1(data):
    # 解析数据
    intervals = []
    coefficients = []
    parsing_intervals = True

    for line in data:
        line = line.strip()
        if line == "===":
            parsing_intervals = False
            continue

        if line:
            parts = line.split()
            if len(parts) == 2:
                interval = parts[0]
                coef = parts[1]
                intervals.append(tuple(map(float, interval.split(','))))
                coefficients.append(list(map(float, coef.split())))

    # 绘制原函数
    x = np.linspace(-1, 1, 400)
    y = 1 / (1 + 25 * x * x)
    plt.plot(x, y, label='Original function f(x)')

    # 绘制样条
    for (a, b), coef in zip(intervals, coefficients):
        x = np.linspace(a, b, 100)
        y = np.polyval(coef[::-1], x)
        plt.plot(x, y, label=f'Spline [{a}, {b}]')

    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Plot for CheckPoint 1')
    plt.legend()
    plt.savefig('plot_checkpoint_1.png')
    plt.show()

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 plotcheck.py <data_file>")
        sys.exit(1)

    data_file = sys.argv[1]

    with open(data_file, 'r') as file:
        lines = file.readlines()

    checkpoint = lines[0].strip()
    data = lines[1:]

    if checkpoint == "CheckPoint 1":
        plot_checkpoint_1(data)
    else:
        print(f"Unknown checkpoint: {checkpoint}")

if __name__ == "__main__":
    main()