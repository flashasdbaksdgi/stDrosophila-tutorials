

import numpy as np
import pyvista as pv

with open("/home/yao/BGIpy37_pytorch113/convolutional_occupancy_networks/out/demo_syn_room/generation/meshes/room_0.off", "r") as file:
    all_content = [i.strip().split(" ") for i in file.readlines()]
    vertices = np.array([i for i in all_content if len(i) == 3]).astype(np.float64)
    faces = np.array([i for i in all_content if len(i) >= 4]).flatten().astype(np.int64)

surf = pv.PolyData(vertices, faces)
surf.plot()