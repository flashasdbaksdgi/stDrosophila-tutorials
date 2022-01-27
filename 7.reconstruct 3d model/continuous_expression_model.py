import pyvista
import numpy as np


sphere_a = pyvista.Sphere(center=(1, 0, 0))



sphere_a["values"] = np.array([0] * sphere_a.n_cells)
print(sphere_a["values"])
sphere_b = pyvista.Sphere(center=(0, 1, 0))
sphere_b["values"] = np.array([1] * sphere_b.n_cells)

sphere_c = pyvista.Sphere(center=(0, 0, 1))
sphere_c["values"] = np.array([2] * sphere_c.n_cells)

merged = sphere_a.merge([sphere_b, sphere_c])
print(merged)
print(merged.point_data)
print(merged.cell_data)
p = pyvista.Plotter()
p.add_mesh(merged, scalars="values")
p.show()