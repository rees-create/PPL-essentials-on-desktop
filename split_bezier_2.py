import numpy as np
import matplotlib.pyplot as plt

def bezier_curve(points, t):
    #De Casteljau's algorithm
    if len(points) == 1:
        return points[0]
    else:
        new_points = [((1 - t) * points[i] + t * points[i+1]) for i in range(len(points)-1)]
        return bezier_curve(new_points, t)

def chop_bezier_curve(control_points, num_segments, segment_resolution):
    t_values = np.linspace(0, 1, num_segments * segment_resolution)
    total_tvalues = len(t_values)
    curve_segments = []

    for i in range(total_tvalues):
        n = len(control_points[int(np.floor(i/segment_resolution))]) - 1
        t0, t1 = t_values[int(np.floor(i/segment_resolution))], t_values[int(np.floor((i/segment_resolution)+1))] if i < num_segments - 1 else 1.0
        segment_points = [bezier_curve(control_points[int(np.floor(i/segment_resolution))], t) for t in np.linspace(t0, t1, segment_resolution) ]
        curve_segments.append(segment_points)

    return curve_segments

# Example usage


#control_points = [np.array([0.0, 0.0]), np.array([1.0, 3.0]), np.array([2.0, 0.0]), np.array([3.0,3.0])]
control_points = [
    [np.array([0.0, 0.0]), np.array([1.0, 3.0]), np.array([1.5, 1.0]), np.array([2,0])],
    [np.array([2.0, 0.0]), np.array([3.0, 2.0]), np.array([4.0, 3.0])],
    [np.array([4.0, 3.0]), np.array([5.0, 1.0]), np.array([6.0, 2.0])],
    [np.array([6.0, 2.0]), np.array([7.0, 3.0]), np.array([8.0, 0.0])],
    [np.array([8.0, 0.0]), np.array([9.0, 3.0]), np.array([10.0, 0.0])]
                  ]  # Replace with your own control points

num_segments = 5  # Specify the number of segments for the spline
segment_resolution = 15

curve_segments = chop_bezier_curve(control_points, num_segments, segment_resolution)
#print(curve_segments)
# Plotting the results
for segment in curve_segments:
    x, y = zip(*segment)
    #print(x)
    plt.plot(x, y, 'b-', marker='o')
for control_rig in control_points:
    x, y = zip(*control_rig)
    #print(x)
    plt.plot(x, y, 'ro-', label='Control Points')
    
plt.legend()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Bezier Curve Split into Segments')
plt.show()
