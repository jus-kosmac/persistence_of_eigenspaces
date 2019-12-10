import rips
import towers
import eigenspaces

# Choose a field Z_p for a prime p.
mod = 29

# Generate the points and their images.
points = rips.points_circles_polar([2, 2], [(-4, 0), (4, 0)], [40, 40], [0.01, 0.01])
# points = rips.points_circles_polar([2], [(0, 0)], [80], [0.1])
images = rips.rotate_circles_polar(points, [2, 3], [40, 40])
mapped = rips.mapped_points(points, images, rips.distance_circles_polar)

# Form the filtration of complexes K_i and the domains dom_i.
rips_comp, simp_dict, max_filt_index, epsilon = rips.rips_complex(points, rips.distance_circles_polar)
domain_filt, mapped_simp, domain_dict = rips.domain(rips_comp, simp_dict, mapped)

# Compute persistent homology for the filtration of K_i.
boundary_matrix, simp_to_index, index_to_simp = towers.sparse_boundary_matrix(rips_comp)
reduced_matrix, zero_columns, columns = towers.reduction_sparse_matrix(boundary_matrix, mod)
[hom0, hom1], [hom0_basis, hom1_basis] = \
    towers.persistence_intervals(reduced_matrix, zero_columns, columns, index_to_simp, simp_dict)
hom0_filt_basis = towers.filtration_basis(hom0, max_filt_index)
hom1_filt_basis = towers.filtration_basis(hom1, max_filt_index)

# Compute persistent homology for the filtration of dom_i.
dom_boundary_matrix, dom_simp_to_index, dom_index_to_simp = towers.sparse_boundary_matrix(domain_filt)
dom_reduced_matrix, dom_zero_columns, dom_columns = towers.reduction_sparse_matrix(dom_boundary_matrix, mod)
[dom_hom0, dom_hom1], [dom_hom0_basis, dom_hom1_basis] = \
    towers.persistence_intervals(dom_reduced_matrix, dom_zero_columns, dom_columns, dom_index_to_simp, domain_dict)
dom_hom0_filt_basis = towers.filtration_basis(dom_hom0, max_filt_index)
dom_hom1_filt_basis = towers.filtration_basis(dom_hom1, max_filt_index)

# Map the basis with inclusion and the given map.
inclusion_basis0 = towers.mapped_basis(dom_hom0_basis, dom_index_to_simp, simp_to_index, mod)
inclusion_basis1 = towers.mapped_basis(dom_hom1_basis, dom_index_to_simp, simp_to_index, mod)
mapped_basis0 = towers.mapped_basis(dom_hom0_basis, dom_index_to_simp, simp_to_index, mod, mapped_simp)
mapped_basis1 = towers.mapped_basis(dom_hom1_basis, dom_index_to_simp, simp_to_index, mod, mapped_simp)

# Compute coefficients of the mapped bases with respect to the original basis.
inclusion_coeffs0 = towers.basis_coefficients(inclusion_basis0, hom0_basis, reduced_matrix, mod)
inclusion_coeffs1 = towers.basis_coefficients(inclusion_basis1, hom1_basis, reduced_matrix, mod)
mapped_coeffs0 = towers.basis_coefficients(mapped_basis0, hom0_basis, reduced_matrix, mod)
mapped_coeffs1 = towers.basis_coefficients(mapped_basis1, hom1_basis, reduced_matrix, mod)

# Form a tower of pairs.
inclusion_tower0 = towers.tower_of_pairs(dom_hom0_filt_basis, hom0_filt_basis, inclusion_coeffs0, max_filt_index)
inclusion_tower1 = towers.tower_of_pairs(dom_hom1_filt_basis, hom1_filt_basis, inclusion_coeffs1, max_filt_index)
mapped_tower0 = towers.tower_of_pairs(dom_hom0_filt_basis, hom0_filt_basis, mapped_coeffs0, max_filt_index)
mapped_tower1 = towers.tower_of_pairs(dom_hom1_filt_basis, hom1_filt_basis, mapped_coeffs1, max_filt_index)

# Form the eigenspace tower and transform it to normal form.
eigen_tower0, eigen_tower_basis0 = eigenspaces.eigenspace_tower(inclusion_tower0, mapped_tower0, dom_hom0_filt_basis,
                                                                max_filt_index, 1, mod)
eigen_tower1, eigen_tower_basis1 = eigenspaces.eigenspace_tower(inclusion_tower1, mapped_tower1, dom_hom1_filt_basis,
                                                                max_filt_index, 3, mod)
eigenspaces.tower_normal_form(eigen_tower0, eigen_tower_basis0, max_filt_index, mod)
eigenspaces.tower_normal_form(eigen_tower1, eigen_tower_basis1, max_filt_index, mod)

# Compute the persistent intervals of the eigenspace tower.
intervals0, persist_generators0 = eigenspaces.tower_persistence(eigen_tower0, eigen_tower_basis0, max_filt_index)
intervals1, persist_generators1 = eigenspaces.tower_persistence(eigen_tower1, eigen_tower_basis1, max_filt_index)
intervals0 = eigenspaces.transform_intervals(intervals0, max_filt_index)
intervals1 = eigenspaces.transform_intervals(intervals1, max_filt_index)
persist_cycles0 = eigenspaces.generators_to_cycles(intervals0, persist_generators0, dom_hom0_filt_basis,
                                                   dom_hom0_basis, dom_index_to_simp, mod)
persist_cycles1 = eigenspaces.generators_to_cycles(intervals1, persist_generators1, dom_hom1_filt_basis,
                                                   dom_hom1_basis, dom_index_to_simp, mod)
stretched_intervals0 = eigenspaces.stretched_intervals(intervals0, epsilon)
stretched_intervals1 = eigenspaces.stretched_intervals(intervals1, epsilon)

# Plotting the persistent generators as homology cycles.
rips.visualize_circles_polar(points, images, persist_cycles1)

