

mesh_points = 20;
NUM_TERMS = 2;

coeffs = zeros(NUM_TERMS, 1);

total_sum = 1;
mesh_grid = linspace(0,1,mesh_points);
NUM_POINTS = size(mesh_grid, 1);

for which_term=1:NUM_TERMS,
  for which_point=1:NUM_POINTS,
    new_coeff = mesh_grid(which_point)
    new_sum = total_sum + new_coeff;
    if new_sum > 1, break; end
    coeffs(which_term) = new_coeff;
    total_sum = sum(coeffs);
