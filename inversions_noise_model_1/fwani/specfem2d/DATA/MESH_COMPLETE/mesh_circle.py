#!/usr/bin/env python3

# FOR SPECFEM2D !!!

import fnmatch
import gmsh
import numpy as np
import pygmsh
import yaml
import sys

# ALL ELEMENTARY ENTITIES ARE DEFINED COUNTER-CLOCKWISE AS REQUIRED BY
# SPECFEM2D

def mesh_cylinder(par):
    # READ PARAMETERS
    # ---------------
    mesh_name = par['MESH_NAME']
    x0 = par['X']
    y0 = par['Y']
    square_side = par['SQUARE_SIDE']
    square_curv = par['SQUARE_CURVATURE']
    nex_theta = par['NEL_ANGULAR']
    nex_r = par['NEL_RADIAL']
    if par["CIRCLE_RADIUS"]:
        circle_radius = par['CIRCLE_RADIUS']
    else:
        circle_radius = None

    file_version = 2.2

    # VERIFY PARAMETERS
    # -----------------
    if (nex_theta % 4) != 0:
        print('NEL_ANGULAR must be a multiple of 4\n')
        sys.exit()

    # SET SOME VARIABLES
    # ------------------
    mesh_algorithm = 8
    mesh_order = 1
    mesh_type = 'Progression'
    coeff = 1.0
    arrangement = 'Left'

    if circle_radius:
        ncircles = len(circle_radius)
    else:
        ncircles = 0

    nex_theta = int(nex_theta / 4)
    nex_theta += 1

    for i in range(0, ncircles):
        nex_r[i] += 1

    point = {}
    line = {}
    surface = {}
    volume = {}

    # INITIALIZE MESH
    # ---------------
    G = pygmsh.geo.Geometry()
    MO = G.__enter__()

    # POINTS
    # ------
    # center
    center_point_coor = np.array([x0, y0, 0.0])
    ep = MO.add_point(center_point_coor)
    point[f'center'] = ep

    # square
    a = square_side / 2.0
    b = a * square_curv

    point_coor = [center_point_coor + np.array([-a, a, 0.0]),
                  center_point_coor + np.array([0.0, b, 0.0]),
                  center_point_coor + np.array([a, a, 0.0]),
                  center_point_coor + np.array([b, 0.0, 0.0]),
                  center_point_coor + np.array([a, -a, 0.0]),
                  center_point_coor + np.array([0.0, -b, 0.0]),
                  center_point_coor + np.array([-a, -a, 0.0]),
                  center_point_coor + np.array([-b, 0.0, 0.0])]

    point['square'] = []

    for coor in point_coor:
        ep = MO.add_point(coor)
        point['square'].append(ep)

    # circles
    if circle_radius:
        for r, radius in enumerate(circle_radius):
            a = radius / np.sqrt(2)

            point_coor = [center_point_coor + np.array([-a, a, 0.0]),
                          center_point_coor + np.array([a, a, 0.0]),
                          center_point_coor + np.array([a, -a, 0.0]),
                          center_point_coor + np.array([-a, -a, 0.0])]

            point[f'circle_{r}'] = []

            for coor in point_coor:
                ep = MO.add_point(coor)
                point[f'circle_{r}'].append(ep)

    # LINES
    # -----
    # square sides
    line['square'] = []

    for i in range(0, 7, 2):
        j = i + 2 if i != 6 else 0

        el = MO.add_spline([point['square'][j],
                            point['square'][i+1],
                            point['square'][i]])

        line['square'].append(el)

    # circles arcs
    if circle_radius:
        for r in range(0, ncircles):
            line[f'circle_{r}'] = []

            for i in range(0, 4):
                j = i + 1 if i != 3 else 0

                el = MO.add_circle_arc(point[f'circle_{r}'][j],
                                       point[f'center'],
                                       point[f'circle_{r}'][i])

                line[f'circle_{r}'].append(el)

        # join first circle with square
        line[f'joint_0'] = []

        for j, i in enumerate(range(0, 7, 2)):
            el = MO.add_line(point[f'circle_0'][j], point['square'][i])
            line[f'joint_0'].append(el)

        # join other circles
        for r in range(1, ncircles):
            line[f'joint_{r}'] = []

            for i in range(0, 4):
                el = MO.add_line(point[f'circle_{r}'][i],
                                 point[f'circle_{r-1}'][i])

                line[f'joint_{r}'].append(el)

    # SURFACES
    # --------
    # square
    ello = MO.add_curve_loop([line['square'][2],
                              line['square'][1],
                              line['square'][0],
                              line['square'][3]])

    es = MO.add_surface(ello)
    surface['square'] = es

    # circles
    if circle_radius:
        for r in range(0, ncircles):
            surface[f'circle_{r}'] = []

            for i in range(0, 4):
                j = i + 1 if i != 3 else 0

                if r == 0:
                    ello = MO.add_curve_loop([-line['square'][i],
                                              -line[f'joint_{r}'][j],
                                              line[f'circle_{r}'][i],
                                              line[f'joint_{r}'][i]])

                else:
                    ello = MO.add_curve_loop([-line[f'circle_{r-1}'][i],
                                              -line[f'joint_{r}'][j],
                                              line[f'circle_{r}'][i],
                                              line[f'joint_{r}'][i]])

                es = MO.add_surface(ello)
                surface[f'circle_{r}'].append(es)

    # TRANSFINITE MESH
    # ----------------
    # lines
    for key, lines in line.items():
        if fnmatch.fnmatch(key, 'joint*'):
            r = int(key[-1])
            num_nodes = nex_r[r]
        else:
            num_nodes = nex_theta

        for el in lines:
            MO.set_transfinite_curve(curve=el, num_nodes=num_nodes,
                                     mesh_type=mesh_type, coeff=coeff)

    # surfaces
    for surfaces in surface.values():
        if not isinstance(surfaces, list):
            surfaces = [surfaces]

        for es in surfaces:
            MO.set_transfinite_surface(surface=es, arrangement=arrangement,
                                       corner_pts=[])

        MO.set_recombined_surfaces(surfaces)

    # DEFINE PHYSICAL LINES
    # -------------------------
    if circle_radius:
        MO.add_physical(line[f'circle_{ncircles-1}'][0], label='Top')
        MO.add_physical(line[f'circle_{ncircles-1}'][1], label='Right')
        MO.add_physical(line[f'circle_{ncircles-1}'][2], label='Bottom')
        MO.add_physical(line[f'circle_{ncircles-1}'][3], label='Left')
    else:
        MO.add_physical(line[f'square'][0], label='Top')
        MO.add_physical(line[f'square'][1], label='Right')
        MO.add_physical(line[f'square'][2], label='Bottom')
        MO.add_physical(line[f'square'][3], label='Left')

    # DEFINE PHYSICAL SURFACES
    # -------------------------
    m1 = []
    m1.append(surface['square'])

    if circle_radius:
        for r in range(0, ncircles):
            m1.extend(surface[f'circle_{r}'])

    MO.add_physical(m1, label='M1')

    # GENERATE AND SAVE MESH
    # ----------------------
    ME = MO.generate_mesh(algorithm=mesh_algorithm,
                          order=mesh_order,
                          verbose=True)

    gmsh.option.setNumber('Mesh.MshFileVersion', file_version)
    gmsh.write(mesh_name)

    return


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('usage: ./mesh_circle.py file\n')
        sys.exit(1)
    else:
        parfile_name = sys.argv[1]

    try:
        with open(parfile_name, 'r') as _file:
            par = yaml.load(_file)
    except IOError:
        print('IOError: {} not found.'.format(parfile_name))

    mesh_cylinder(par)
