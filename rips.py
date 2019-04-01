import math as m
import numpy as np

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
    return (a[0] - b[0])**2 + (a[1] - b[1])**2

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
    
if __name__ == '__main__':
    points = points_circle_polar(1, 140)
    images = power_polar(points, 2)
    mapped = mapped_points(points, images, distance_polar)
    rips, simp_dict, max_index = rips_complex(points, distance_polar)
    domain_filt, mapped_simp, domain_dict = domain(rips, simp_dict, mapped)
