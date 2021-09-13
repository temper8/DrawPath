import os
def read_bounds():
    Icms_path = os.path.abspath("data/lcms.dat")
    file = open(Icms_path)

    header = file.readline().split()
    print(header)
    lines = file.readlines()

    table = [line.split() for line in lines]
    table = list(filter(None, table))

    R = [float(row[0]) for row in table]
    Z = [float(row[1]) for row in table]
    return R, Z