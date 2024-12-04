import numpy as np
import matplotlib.pyplot as plt
import os

def read_points(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    
    original_points = []
    plane_points = []
    plane_spline_points = []
    spherical_points = []
    reading_original = False
    reading_plane = False
    reading_plane_spline = False
    reading_spherical = False
    
    for line in lines:
        line = line.strip()
        if line == "original_points:":
            reading_original = True
            reading_plane = False
            reading_plane_spline = False
            reading_spherical = False
        elif line == "plane_points:":
            reading_original = False
            reading_plane = True
            reading_plane_spline = False
            reading_spherical = False
        elif line == "plane_spline:":
            reading_original = False
            reading_plane = False
            reading_plane_spline = True
            reading_spherical = False
        elif line == "spherical_points:":
            reading_original = False
            reading_plane = False
            reading_plane_spline = False
            reading_spherical = True
        elif reading_original:
            original_points.append(list(map(float, line.split(','))))
        elif reading_plane:
            plane_points.append(list(map(float, line.split(','))))
        elif reading_plane_spline:
            plane_spline_points.append(list(map(float, line.split(','))))
        elif reading_spherical:
            spherical_points.append(list(map(float, line.split(','))))
    
    return original_points, plane_points, plane_spline_points, spherical_points

def plot_points(original_points, plane_points, plane_spline_points, spherical_points, output_file):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # 绘制原始点
    original_points = np.array(original_points)
    ax.scatter(original_points[:, 0], original_points[:, 1], original_points[:, 2], c='r', label='Original Points')
    
    # 绘制平面点
    plane_points = np.array(plane_points)
    ax.scatter(plane_points[:, 0], plane_points[:, 1], plane_points[:, 2], c='g', label='Plane Points')
    
    # 绘制平面样条点
    plane_spline_points = np.array(plane_spline_points)
    ax.plot(plane_spline_points[:, 0], plane_spline_points[:, 1], plane_spline_points[:, 2], c='b', label='Plane Spline Points')
    
    # 绘制球面样条点
    spherical_points = np.array(spherical_points)
    ax.plot(spherical_points[:, 0], spherical_points[:, 1], spherical_points[:, 2], c='m', label='Spherical Spline Points')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('3D Spline Fitting')
    plt.legend()
    plt.savefig(output_file)
    plt.close()

def plot_sphere(ax, radius=1, center=(0, 0, 0)):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
    y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]
    ax.plot_surface(x, y, z, color='c', alpha=0.3)

def plot_spherical_points(original_points, spherical_points, output_file):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # 绘制球面
    plot_sphere(ax, radius=1, center=(0, 0, 1))
    
    # 绘制原始点
    original_points = np.array(original_points)
    ax.scatter(original_points[:, 0], original_points[:, 1], original_points[:, 2], c='r', label='Original Points')
    
    # 绘制球面样条点
    spherical_points = np.array(spherical_points)
    ax.plot(spherical_points[:, 0], spherical_points[:, 1], spherical_points[:, 2], c='m', label='Spherical Spline Points')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('Spherical Spline Fitting')
    plt.legend()
    plt.savefig(output_file)
    plt.close()

def main():
    input_file = 'output/Sphere/output.txt'
    output_dir = 'figure/Sphere'
    os.makedirs(output_dir, exist_ok=True)
    
    original_points, plane_points, plane_spline_points, spherical_points = read_points(input_file)
    
    # 绘制包含所有点的图像
    output_file_all = os.path.join(output_dir, 'spline_fitting_all.png')
    plot_points(original_points, plane_points, plane_spline_points, spherical_points, output_file_all)
    
    # 绘制只包含原始点和球面样条点的图像
    output_file_spherical = os.path.join(output_dir, 'spline_fitting_spherical.png')
    plot_spherical_points(original_points, spherical_points, output_file_spherical)

if __name__ == "__main__":
    main()