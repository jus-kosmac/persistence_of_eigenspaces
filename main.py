import rips
import towers
import eigenspaces

#TODO: v vseh datotekah popravi komentarje

# obseg Z_p
mod = 29

# generiramo tocke in ripsov kompleks
points = rips.points_circle_polar(1, 50, 0.03)
images = rips.power_polar(points, 2)
mapped = rips.mapped_points(points, images, rips.distance_polar)

rips_comp, simp_dict, max_filt_index = rips.rips_complex(points, rips.distance_polar)
domain_filt, mapped_simp, domain_dict = rips.domain(rips_comp, simp_dict, mapped)

# za K_i
boundary_matrix, simp_to_index, index_to_simp = towers.sparse_boundary_matrix(rips_comp)
reduced_matrix, zero_columns, columns = towers.reduction_sparse_matrix(boundary_matrix, mod)
hom0, hom1, hom0_basis, hom1_basis = towers.persistence_intervals(reduced_matrix, zero_columns, columns,
                                                                  index_to_simp, simp_dict)
hom0_filt_basis = towers.filtration_basis(hom0, max_filt_index)
hom1_filt_basis = towers.filtration_basis(hom1, max_filt_index)

# za dom_i
dom_boundary_matrix, dom_simp_to_index, dom_index_to_simp = towers.sparse_boundary_matrix(domain_filt)
dom_reduced_matrix, dom_zero_columns, dom_columns = towers.reduction_sparse_matrix(dom_boundary_matrix, mod)
dom_hom0, dom_hom1, dom_hom0_basis, dom_hom1_basis = \
    towers.persistence_intervals(dom_reduced_matrix, dom_zero_columns, dom_columns, dom_index_to_simp, domain_dict)
dom_hom0_filt_basis = towers.filtration_basis(dom_hom0, max_filt_index)
dom_hom1_filt_basis = towers.filtration_basis(dom_hom1, max_filt_index)

# preslikamo bazo
inclusion_basis0 = towers.mapped_basis(dom_hom0_basis, dom_index_to_simp, simp_to_index, mod)
inclusion_basis1 = towers.mapped_basis(dom_hom1_basis, dom_index_to_simp, simp_to_index, mod)
mapped_basis0 = towers.mapped_basis(dom_hom0_basis, dom_index_to_simp, simp_to_index, mod, mapped_simp)
mapped_basis1 = towers.mapped_basis(dom_hom1_basis, dom_index_to_simp, simp_to_index, mod, mapped_simp)

# razvijemo preslikano bazo po originalni bazi
inclusion_coeffs0 = towers.basis_coefficients(inclusion_basis0, hom0_basis, reduced_matrix, mod)
inclusion_coeffs1 = towers.basis_coefficients(inclusion_basis1, hom1_basis, reduced_matrix, mod)
mapped_coeffs0 = towers.basis_coefficients(mapped_basis0, hom0_basis, reduced_matrix, mod)
mapped_coeffs1 = towers.basis_coefficients(mapped_basis1, hom1_basis, reduced_matrix, mod)

# naredimo osnovna stolpa
inclusion_tower0 = towers.tower_of_pairs(dom_hom0_filt_basis, hom0_filt_basis, inclusion_coeffs0, max_filt_index)
inclusion_tower1 = towers.tower_of_pairs(dom_hom1_filt_basis, hom1_filt_basis, inclusion_coeffs1, max_filt_index)
mapped_tower0 = towers.tower_of_pairs(dom_hom0_filt_basis, hom0_filt_basis, mapped_coeffs0, max_filt_index)
mapped_tower1 = towers.tower_of_pairs(dom_hom1_filt_basis, hom1_filt_basis, mapped_coeffs1, max_filt_index)

# eigenspace tower
eigen_tower0, eigen_tower_basis0 = eigenspaces.eigenspace_tower(inclusion_tower0, mapped_tower0, dom_hom0_filt_basis,
                                                                max_filt_index, 2, mod)
eigen_tower1, eigen_tower_basis1 = eigenspaces.eigenspace_tower(inclusion_tower1, mapped_tower1, dom_hom1_filt_basis,
                                                                max_filt_index, 3, mod)
eigenspaces.tower_normal_form(eigen_tower0, eigen_tower_basis0, max_filt_index, mod)
eigenspaces.tower_normal_form(eigen_tower1, eigen_tower_basis1, max_filt_index, mod)
