import math as m
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as plt3d

def points_circle(r, n, sigma = 0):
    points = []
    for i in range(n):
        points.append((r * m.sin(2 * m.pi * i / n) + np.random.normal(0, sigma),
                       r * m.cos(2 * m.pi * i / n) + np.random.normal(0, sigma)))
    return points

def points_circle_polar(r, n, sigma = 0):
    points = []
    for i in range(n):
        points.append((max(r + np.random.normal(0, sigma), 0),
                       2 * m.pi * i / n + np.random.normal(0, sigma)))
    return points

def power_polar(points, eksp):
    return [(r**eksp, eksp * fi) for (r, fi) in points]

def distance_polar(a, b):
    x = (a[0] * m.sin(a[1]), a[0] * m.cos(a[1]))
    y = (b[0] * m.sin(b[1]), b[0] * m.cos(b[1]))
    return distance(x, y)

def distance(a, b):
    dist = 0
    for i in range(len(a)):
        dist += (a[i] - b[i]) ** 2
    return dist

def points_torus_polar(R, r, N, n, sigma = 0):
    points = []
    for i in range(N):
        for j in range(n):
            big = max(R + np.random.normal(0, sigma), 0)
            small = min(max(r + np.random.normal(0, sigma), 0), big)
            points.append((big, small, 2 * m.pi * i / N + np.random.normal(0, sigma),
                           2 * m.pi * j / n + np.random.normal(0, sigma)))

    return points

def rotate_torus_polar(points, eksp_big, eksp_small):
    return [(R, r, fi * eksp_big, theta * eksp_small) for (R, r, fi, theta) in points]

def distance_torus_polar(a, b):
    x = ((a[0] + a[1] * m.cos(a[3])) * m.cos(a[2]), (a[0] + a[1] * m.cos(a[3])) * m.sin(a[2]), a[1] * m.sin(a[3]))
    y = ((b[0] + b[1] * m.cos(b[3])) * m.cos(b[2]), (b[0] + b[1] * m.cos(b[3])) * m.sin(b[2]), b[1] * m.sin(b[3]))
    return distance(x, y)

def points_circles_polar(radii, centres, numbers, sigmas = None):
    #vsi parametri so seznami enakih dolzin
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
    #number in factor sta seznama enake dolzine, povemo koliko tock zavrtimo za posamezni faktor
    #skupna vsota stevil v seznamu number mora biti ravno dolzina seznama points
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
    (x1, y1), r1, fi1 = a
    (x2, y2), r2, fi2 = b
    first = (x1 + r1 * m.cos(fi1), y1 + r1 * m.sin(fi1))
    second = (x2 + r2 * m.cos(fi2), y2 + r2 * m.sin(fi2))
    return distance(first, second)


def rips_complex(points, distance):
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
                
    simplices.sort()
    rips = []
    current_dist = 0
    filtration_index = 1
    simp_dict = dict() #kljuc: simplex, vrednost: index filtracije
    
    for (d, _, simp) in simplices:
        if d > current_dist:
            current_dist = d
            filtration_index += 1
        rips.append((filtration_index, simp))
        simp_dict[simp] = filtration_index
        
    return rips, simp_dict, filtration_index

def mapped_points(points, images, distance):
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
    mapped_simp = dict() #kljuc: simplex, vrednost: (preslikani simplex, signatura, collapse)
    domain_dict = dict() #kljuc: simplex, vrednost: index filtracije domene
    domain_filtration = []

    signature = {(0,):1, (0,1):1, (1,0):-1, (0,1,2):1, (1,0,2):-1, (0,2,1):-1, (1,2,0):1, (2,1,0):-1, (2,0,1):1}
    
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

        mapped_simp[simp] = (image_simp, signature[permutation], collapse)
        domain_index = max(simp_dict[simp], simp_dict[image_simp])
        domain_dict[simp] = domain_index
        domain_filtration.append((domain_index, simp))
        
    return sorted(domain_filtration, key=lambda x: (x[0], len(x[1]), x[1])), mapped_simp, domain_dict

def visualize_circle_polar(points, images, cycles = False):
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

    if cycles is not False: #cycles je seznam ciklov v prvi homologiji
        for cycle in cycles:
            for start, end in cycle:
                plt.plot([points_x[start], points_x[end]], [points_y[start], points_y[end]], c='g')

    plt.scatter(points_x, points_y, c='b')
    plt.scatter(images_x, images_y, c='r')

def visualize_circles_polar(points, images, cycles = False):
    points_x, points_y = [], []
    images_x, images_y = [], []

    for (a, b), r, fi in points:
        points_x.append(a + r * m.cos(fi))
        points_y.append(b + r * m.sin(fi))

    for (a, b), r, fi in images:
        images_x.append(a + r * m.cos(fi))
        images_y.append(b + r * m.sin(fi))

    plt.axes().set_aspect('equal')

    if cycles is not False: #cycles je seznam ciklov v prvi homologiji
        for cycle in cycles:
            for start, end in cycle: #zacetek in konec 1-simpleksa, ki tvori cikel
                plt.plot([points_x[start], points_x[end]], [points_y[start], points_y[end]], c='g')

    plt.scatter(points_x, points_y, c='b')
    plt.scatter(images_x, images_y, c='r')


def visualize_torus_polar(points, images, cycles = False):
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
    points = points_circles_polar([1, 3, 2], [(-5, 0), (2, 2), (3, 4)], [20, 40, 30], [0.03, 0.1, 0.03])
    images = rotate_circles_polar(points, [2, -1, 3], [20, 40, 30])
    mapped = mapped_points(points, images, distance_circles_polar)
    rips, simp_dict, max_index = rips_complex(points, distance_circles_polar)
    domain_filt, mapped_simp, domain_dict = domain(rips, simp_dict, mapped)

    # plotting
    cycle = {(1, 10): 1, (10, 50): 1, (40, 80): 1}
    visualize_circles_polar(points, images, [cycle])
