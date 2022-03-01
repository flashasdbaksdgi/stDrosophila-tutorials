

import numpy as np
import pyvista as pv
import vtk

mesh = pv.read('./tree.ply')
grid = mesh.cast_to_unstructured_grid()

pl = pv.Plotter()
pl.add_mesh(grid, color='w')


# picking example
picked = []
legend = []


def split_mesh(selected_grid):
    """Adds a new mesh to the plotter each time cells are picked, and
    removes them from the original mesh"""

    # if nothing selected
    if not selected_grid.n_cells:
        return

    # remove the picked cells from main grid
    ghost_cells = np.zeros(grid.n_cells, np.uint8)
    ghost_cells[selected_grid['orig_extract_id']] = 1
    grid.cell_arrays[vtk.vtkDataSetAttributes.GhostArrayName()] = ghost_cells
    grid.RemoveGhostCells()

    # add the selected mesh this to the main plotter
    color = np.random.random(3)
    legend.append(['picked mesh %d' % len(picked), color])
    pl.add_mesh(selected_grid, color=color)
    pl.add_legend(legend)

    # track the picked meshes and label them
    selected_grid['picked_index'] = np.ones(selected_grid.n_points)*len(picked)
    picked.append(selected_grid)


# enable cell picking with our custom callback
pl.enable_cell_picking(mesh=grid, callback=split_mesh, show=False)
pl.show()
# convert these meshes back to surface meshes (PolyData)
meshes = []
for selected_grid in picked:
    meshes.append(selected_grid.extract_surface())


# plot final separated meshes for fun
pv.plot(meshes)

"""
# save these meshes somewhere
for i, mesh in enumerate(meshes):
    mesh.save('mesh_%03d.ply' % i)
"""