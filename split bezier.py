import numpy as np
import matplotlib.pyplot as plt

def bezier_curve(points, t):
    if len(points) == 1:
        return points[0]
    else:
        new_points = [((1 - t) * points[i] + t * points[i+1]) for i in range(len(points)-1)]
        return bezier_curve(new_points, t)

def continuous_spline(points, num_components):
    t_values = np.linspace(0, 1, num_components)
    spline_points = [bezier_curve(points, t) for t in t_values]
    return spline_points

# Example usage
control_points = [np.array([0, 0]), np.array([1, 3]), np.array([2,1]), np.array([3,3])]  # Replace with your own control points
num_components = 10  # Specify the number of components for the spline

spline_points = continuous_spline(control_points, num_components)
#print(spline_points)

# Plotting the results
x, y = zip(*control_points)
print(x)
plt.plot(x, y, 'ro-', label='Control Points')
x_spline, y_spline = zip(*spline_points)
print(x_spline)
plt.plot(x_spline, y_spline, 'b-', label=f'Continuous Spline ({num_components} components)')
plt.legend()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Bezier Curve and Continuous Spline')
plt.show()
