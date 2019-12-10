import math as m
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as plt3d
import sympy.combinatorics.permutations as sp


def points_circle(r, n, sigma=0):
    """Return a list of n equidistant points on a circle of radius r in cartesian coordinates.
    Also add Gaussian noise with variance sigma.
    """
    points = []
    for i in range(n):
        points.append((r * m.sin(2 * m.pi * i / n) + np.random.normal(0, sigma),
                       r * m.cos(2 * m.pi * i / n) + np.random.normal(0, sigma)))
    return points


def points_circle_polar(r, n, sigma=0):
    """Return a list of n equidistant points on a circle of radius r in polar coordinates.
    Also add Gaussian noise with variance sigma.
    """
    points = []
    for i in range(n):
        points.append((max(r + np.random.normal(0, sigma), 0),
                       2 * m.pi * i / n + np.random.normal(0, sigma)))
    return points


def power_polar(points, eksp):
    """Map the function z -> z^eksp to a list of points in polar coordinates."""
    return [(r**eksp, eksp * fi) for (r, fi) in points]


def distance_polar(a, b):
    """Return the squared euclidean distance between two points in polar coordinates."""
    x = (a[0] * m.sin(a[1]), a[0] * m.cos(a[1]))
    y = (b[0] * m.sin(b[1]), b[0] * m.cos(b[1]))
    return distance(x, y)


def distance(a, b):
    """Return the squared euclidean distance between two points in cartesian coordinates."""
    dist = 0
    for i in range(len(a)):
        dist += (a[i] - b[i]) ** 2
    return dist


def points_torus_polar(R, r, N, n, sigma=0):
    """Return a list of N x n points on a torus with big circle of radius r and small circle of radius r
    in torus polar coordinates. Also add Gaussian noise with variance sigma.
    """
    points = []
    for i in range(N):
        for j in range(n):
            big = max(R + np.random.normal(0, sigma), 0)
            small = min(max(r + np.random.normal(0, sigma), 0), big)
            points.append((big, small, 2 * m.pi * i / N + np.random.normal(0, sigma),
                           2 * m.pi * j / n + np.random.normal(0, sigma)))

    return points


def rotate_torus_polar(points, eksp_big, eksp_small):
    """Map the function, which rotates the big and the small circle of torus eksp_big and eksp_small times,
    to a list of points in torus polar coordinates.
    """
    return [(R, r, fi * eksp_big, theta * eksp_small) for (R, r, fi, theta) in points]


def distance_torus_polar(a, b):
    """Return the squared euclidean distance between two points in torus polar coordinates."""
    x = ((a[0] + a[1] * m.cos(a[3])) * m.cos(a[2]), (a[0] + a[1] * m.cos(a[3])) * m.sin(a[2]), a[1] * m.sin(a[3]))
    y = ((b[0] + b[1] * m.cos(b[3])) * m.cos(b[2]), (b[0] + b[1] * m.cos(b[3])) * m.sin(b[2]), b[1] * m.sin(b[3]))
    return distance(x, y)


def points_circles_polar(radii, centres, numbers, sigmas=None):
    """Return a list of points in offset polar coordinates (centre, radius, angle). Also add Gaussian noise.

    Parameters:
    radii -- a list of radii of circles
    centres -- a list of centres of circles
    numbers -- a list of numbers of points to sample from each circle
    sigmas -- None (no noise) or a list of variances of Gaussian noise for each circle
    All lists have to be of equal length.
    """
    points = []
    for i in range(len(radii)):
        R = radii[i]
        c = centres[i]
        n = numbers[i]
        sig = 0 if sigmas is None else sigmas[i]
        for j in range(n):
            points.append((c, max(R + np.random.normal(0, sig), 0), 2 * m.pi * j / n + np.random.normal(0, sig)))

    return points


def rotate_circles_polar(points, factors, numbers):
    """Map a function which rotates each circle to a list of points in offset polar coordinates.

    Parameters:
    factors -- a list of factors to rotate for
    numbers -- a list of numbers of points which should be rotated with each factor
    Factors and numbers are lists of equal length. The sum of all numbers has to equal the length of points.
    """
    rotated_points = [None for _ in range(len(points))]
    index_correction = 0
    for i in range(len(numbers)):
        n = numbers[i]
        factor = factors[i]
        for j in range(n):
            c, R, fi = points[index_correction + j]
            rotated_points[index_correction + j] = (c, R, factor * fi)

        index_correction += n

    return rotated_points


def distance_circles_polar(a, b):
    """Return the squared euclidean distance between two points in offset polar coordinates."""
    (x1, y1), r1, fi1 = a
    (x2, y2), r2, fi2 = b
    first = (x1 + r1 * m.cos(fi1), y1 + r1 * m.sin(fi1))
    second = (x2 + r2 * m.cos(fi2), y2 + r2 * m.sin(fi2))
    return distance(first, second)


def rips_complex(points, distance, dim3=False):
    """Return the complete filtration of Vietoris-Rips complexes on a given list of points.

    Parameters:
    points -- a list of points (in any coordinates)
    distance -- a distance function between two points (takes as arguments points as given in the list points)
    dim2 -- True if we want to compute also 3D-simplices, False if we only want 0D, 1D and 2D-simplices

    Return:
    rips -- a list of pairs (ordered simplex, minimal filtration index) sorted by increasing minimal filtration index,
            followed by increasing dimension of simplices
    simp_dict -- a dictionary: key = ordered simplex, value = minimal filtration index
    filtration_index -- maximal filtration index
    epsilon -- a list of values of distances for each filtration index
    """
    n = len(points)
    simplices = []
    
    for i in range(n):
        simplices.append((0, 1, (i,)))
    for i in range(n):
        for j in range(i + 1, n):
            d = distance(points[i], points[j])
            simplices.append((d, 2, (i, j)))
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                d = max(distance(points[i], points[j]), distance(points[i], points[k]), distance(points[j], points[k]))
                simplices.append((d, 3, (i, j, k)))

    if dim3:
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    for l in range(k + 1, n):
                        d = max(distance(points[i], points[j]), distance(points[i], points[k]),
                                distance(points[i], points[l]), distance(points[j], points[k]),
                                distance(points[j], points[l]), distance(points[k], points[l]))
                        simplices.append((d, 4, (i, j, k, l)))
                
    simplices.sort()
    rips = []
    epsilon = [0]
    current_dist = 0
    filtration_index = 0
    simp_dict = dict()
    
    for (d, _, simp) in simplices:
        if d > current_dist:
            current_dist = d
            filtration_index += 1
            epsilon.append(m.sqrt(current_dist))
        rips.append((filtration_index, simp))
        simp_dict[simp] = filtration_index
        
    return rips, simp_dict, filtration_index, epsilon


def mapped_points(points, images, distance):
    """Return a list of indices of closest original point to for each point in images.

    Parameters:
    points -- a list of points (in any coordinates)
    images -- a list of points (in the same coordinates as in the list points)
    distance -- a distance function between two points (takes as arguments points as given in the list points)

    Return:
    mapped_pts -- mapped_pts[i] = j iff points[j] is the closest point in points to images[i]
    """
    mapped_pts = []
    for i in images:
        closest_dist = None
        closest_point = None
        
        for (index, p) in enumerate(points):
            d = distance(p, i)
            if closest_dist is None or d < closest_dist:
                closest_dist = d
                closest_point = index
                
        mapped_pts.append(closest_point)
        
    return mapped_pts


def domain(rips, simp_dict, mapped_pts):
    """Return the complete filtration of domains from the given filtration of complexes and list of mapped vertices.

    Parameters:
    rips -- as in the function rips_complex
    simp_dict -- as in the function rips_complex
    mapped_pts -- as in the function mapped_points

    Return:
    domain_filtration -- a list of pairs (ordered simplex, minimal domain filtration index) sorted by increasing minimal
                         domain filtration index, followed by increasing dimension of simplices
    mapped_simp -- a dictionary: key = ordered simplex, value = (mapped simplex, signature : 1 or -1, collapse : Bool)
                   signature tells us if we preserve or reverse orientation of the simplex
                   collapse tells us if the mapped simplex collapses in dimension
    domain_dict -- a dictionary, key = ordered simplex, value = minimal domain filtration index
    """
    mapped_simp = dict()
    domain_dict = dict()
    domain_filtration = []
    
    for (index, simp) in rips:
        image = []
        collapse = False

        for i in range(len(simp)):
            if mapped_pts[simp[i]] not in image:
                image.append(mapped_pts[simp[i]])
            else:
                collapse = True

        image = sorted(list(enumerate(image)), key=lambda x: (x[1], x[0]))
        image_simp = tuple([s for (_, s) in image])
        permutation = tuple([i for (i, _) in image])

        mapped_simp[simp] = (image_simp, sp.Permutation(permutation).signature(), collapse)
        domain_index = max(simp_dict[simp], simp_dict[image_simp])
        domain_dict[simp] = domain_index
        domain_filtration.append((domain_index, simp))
        
    return sorted(domain_filtration, key=lambda x: (x[0], len(x[1]), x[1])), mapped_simp, domain_dict


def visualize_circle_polar(points, images, cycles=False):
    """Plot a list of points and images in polar coordinates. Cycles (if not False) is a list of 1-cycles in
    homology, given as dictionary with keys = pairs of indices (index of start point, index of end point) and
    values = coefficient of simplex. We also plot the cycles.
    """
    points_x, points_y = [], []
    images_x, images_y = [], []

    for r, fi in points:
        points_x.append(r * m.cos(fi))
        points_y.append(r * m.sin(fi))

    for r, fi in images:
        images_x.append(r * m.cos(fi))
        images_y.append(r * m.sin(fi))

    # x_min, y_min = min(points_x + images_x), min(points_y + images_y)
    # x_max, y_max = max(points_x + images_x), max(points_y + images_y)
    #
    # factor = 1.2
    # lower = factor * min(x_min, y_min)
    # upper = factor * max(x_max, y_max)
    #
    # plt.xlim(lower, upper)
    # plt.ylim(lower, upper)
    plt.axes().set_aspect('equal')
    # plt.axis('equal')

    if cycles is not False:
        for cycle in cycles:
            for start, end in cycle:
                plt.plot([points_x[start], points_x[end]], [points_y[start], points_y[end]], c='g')

    plt.scatter(points_x, points_y, c='b')
    plt.scatter(images_x, images_y, c='r')


def visualize_circles_polar(points, images, cycles=False):
    """Plot a list of points and images in offset polar coordinates. Cycles (if not False) is a list of 1-cycles in
    homology, given as dictionary with keys = pairs of indices (index of start point, index of end point) and
    values = coefficient of simplex. We also plot the cycles.
    """
    points_x, points_y = [], []
    images_x, images_y = [], []

    for (a, b), r, fi in points:
        points_x.append(a + r * m.cos(fi))
        points_y.append(b + r * m.sin(fi))

    for (a, b), r, fi in images:
        images_x.append(a + r * m.cos(fi))
        images_y.append(b + r * m.sin(fi))

    plt.axes().set_aspect('equal')

    if cycles is not False:
        for cycle in cycles:
            for start, end in cycle:
                plt.plot([points_x[start], points_x[end]], [points_y[start], points_y[end]], c='g')

    plt.scatter(points_x, points_y, c='b')
    plt.scatter(images_x, images_y, c='r')


def visualize_torus_polar(points, images, cycles=False):
    """Plot a list of points and images in torus polar coordinates. Cycles (if not False) is a list of 1-cycles in
    homology, given as dictionary with keys = pairs of indices (index of start point, index of end point) and
    values = coefficient of simplex. We also plot the cycles.
    """
    points_x, points_y, points_z = [], [], []
    images_x, images_y, images_z = [], [], []

    for R, r, fi, theta in points:
        points_x.append((R + r * m.cos(theta)) * m.cos(fi))
        points_y.append((R + r * m.cos(theta)) * m.sin(fi))
        points_z.append(r * m.sin(theta))

    for R, r, fi, theta in images:
        images_x.append((R + r * m.cos(theta)) * m.cos(fi))
        images_y.append((R + r * m.cos(theta)) * m.sin(fi))
        images_z.append(r * m.sin(theta))

    x_min, y_min, z_min = min(points_x + images_x), min(points_y + images_y), min(points_z + images_z)
    x_max, y_max, z_max = max(points_x + images_x), max(points_y + images_y), max(points_z + images_z)

    lower = min(x_min, y_min, z_min)
    upper = max(x_max, y_max, z_max)

    fig = plt.figure()
    ax = plt3d.Axes3D(fig)
    ax.set_xlim3d(lower, upper)
    ax.set_ylim3d(lower, upper)
    ax.set_zlim3d(lower, upper)

    if cycles is not False:
        for cycle in cycles:
            for start, end in cycle:
                ax.plot([points_x[start], points_x[end]], [points_y[start], points_y[end]],
                        [points_z[start], points_z[end]], c='g')

    ax.scatter(points_x, points_y, points_z, c='b')
    ax.scatter(images_x, images_y, images_z, c='r')


if __name__ == '__main__':
    points = points_circles_polar([1, 3, 2], [(-5, 0), (2, 2), (3, 4)], [10, 15, 5], [0.03, 0.1, 0.03])
    images = rotate_circles_polar(points, [2, -1, 3], [10, 15, 5])
    mapped = mapped_points(points, images, distance_circles_polar)
    rips, simp_dict, max_index, epsilon = rips_complex(points, distance_circles_polar, True)
    domain_filt, mapped_simp, domain_dict = domain(rips, simp_dict, mapped)

    # Plotting.
    cycle = {(1, 10): 1, (10, 15): 1, (15, 29): 1}
    visualize_circles_polar(points, images, [cycle])
