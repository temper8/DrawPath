import os
from collections import namedtuple

def read_bounds(f):
    Icms_path = os.path.abspath(f)
    file = open(Icms_path)

    header = file.readline().split()
    print(header)
    lines = file.readlines()

    table = [line.split() for line in lines]
    table = list(filter(None, table))

    R = [float(row[0]) for row in table]
    Z = [float(row[1]) for row in table]
    return R, Z


def read_trajectories(f):
    Icms_path = os.path.abspath(f)
    file = open(Icms_path)

    header = file.readline().split()
    lines = file.readlines()

    table = [line.split() for line in lines]
    table = list(filter(None, table))

    rays = []
    N_traj = 0
    TRay = namedtuple('Ray' , 'R Z N_traj')
    for row in table:
        if N_traj != int(row[12]):
            N_traj = int(row[12])
            ray = TRay ([], [], N_traj)
            rays.append(ray)
        ray.R.append(float(row[0]))
        ray.Z.append(float(row[1]))

    return rays, N_traj

def read_trajectories2(f):
    Icms_path = os.path.abspath(f)
    file = open(Icms_path)

    header = file.readline().replace('=', '_').split()

    print(header)

    empty_ray = dict([ (h, []) for h in header ])

    lines = file.readlines()
    table = [line.split() for line in lines]
    table = list(filter(None, table))

    rays = []
    N_traj = 0
    #header.append('N_traj')
    TRay = namedtuple('Ray' , header)

    for row in table:
        if N_traj != int(row[12]):
            N_traj = int(row[12])
            ray = dict([ (h, []) for h in header ])
            rays.append(ray)
        for index, (p, item) in enumerate(ray.items()):
            item.append(float(row[index]))
        #    print(p, item, index)
        # break


    return rays, N_traj
