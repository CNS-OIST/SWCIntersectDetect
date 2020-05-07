import numpy as np
from typing import Tuple, List

try:
    from neuron import h
    HAS_NEURON = True
except ImportError:
    print("Unable to import NEURON module, NEURON related functions will not be usable.")
    HAS_NEURON = False

SOMA = 1
AXON = 2
DEND = 3
APIC = 4

def h2morph(h):
    """
    Generate morph sectioning data from NEURON h interface with preloaded morphology.
    """
    sections = {}

    for sec in h.allsec():
        sections[sec.name()] = {}
        sections[sec.name()]["name"] = sec.name()
        sr = h.SectionRef(sec=sec)
        if sr.has_parent():
            parent = sr.parent.name()
        else:
            parent = None
        sections[sec.name()]["parent"] = parent

        children = []
        for child in sr.child:
            children.append(child.name())

        sections[sec.name()]["children"] = children
        n3d = int(h.n3d())
        sections[sec.name()]["points"] = []
        for i in range(n3d):
            sections[sec.name()]["points"].append([h.x3d(i),
                                                   h.y3d(i),
                                                   h.z3d(i),
                                                   h.diam3d(i)])
    return sections


def nrn2morph(nrn_file : str):
    """
    Generate morph sectioning data from NEURON .nrn file.
    """

    if HAS_NEURON:
        h.load_file(nrn_file)
        return h2morph(h)
    else:
        print("NEURON module is not available.")


def swc2morph(swc_file : str):
    """
    Generate morph sectioning data from .swc file.
    """
    if HAS_NEURON:
        h.load_file('stdlib.hoc')
        h.load_file('import3d.hoc')

        cell = h.Import3d_SWC_read()
        cell.input(swc_file)

        i3d = h.Import3d_GUI(cell, 0)
        i3d.instantiate(None)
        return neuron2morph(h)

    else:
        print("NEURON module is not available.")

def nrn2swc(nrn_file : str, swc_output : str):
    """
    Convert NEURON .nrn morphology file to .swc file.
    """
    morph_data = nrn2morph(nrn_file)
    root = None
    for branch in morph_data.keys():
        if morph_data[branch]['parent'] == None:
            root = branch
            break

    swc_file = open(swc_output, 'w')
    swc_file.write("# ORIGINAL_FILE %s\n" % nrn_file)
    swc_file.write("# UNIT micrometer\n")
    swc_file.write("\n")
    current_id = 1
    remain_branches = [root]
    point_index = {}
    while len(remain_branches) != 0:
        current_branch = remain_branches[0]
        remain_branches.remove(current_branch)
        remain_branches.extend(morph_data[current_branch]['children'])
        type = 0
        if "soma" in current_branch:
            type = SOMA
        elif "axon" in current_branch:
            type = AXON
        elif "dend" in current_branch:
            type = DEND
        elif "apic" in current_branch:
            type = APIC
        parent = morph_data[current_branch]['parent']
        points = morph_data[current_branch]['points']
        if current_branch == root and len(point_index) == 0:
            swc_file.write("%i %i %f %f %f %f %i\n" % (current_id, type, points[0][0], points[0][1], points[0][2], points[0][3] / 2.0, -1))
            parent_point_id = current_id
            point_index[tuple(points[0])] = current_id
            current_id += 1
        else:
            parent_point_id = point_index[tuple(morph_data[parent]['points'][-1])]
            swc_file.write("%i %i %f %f %f %f %i\n" % (current_id, type, points[0][0], points[0][1], points[0][2], points[0][3] / 2.0, parent_point_id))
            parent_point_id = current_id
            point_index[tuple(points[0])] = current_id
            current_id += 1
        for point in points[1:]:
            swc_file.write("%i %i %f %f %f %f %i\n" % (current_id, type, point[0], point[1], point[2], point[3] / 2.0, parent_point_id))
            point_index[tuple(point)] = current_id
            parent_point_id = current_id
            current_id += 1
    swc_file.close()

def read_swc(file_path : str) -> Tuple[dict, List[str]]:
    """
    Read a .swc morphology file, and return its data.
    """
    file = open(file_path, "r")
    data = {}
    info = []
    for line in file:
        if line[0] == '#':
            info.append(line)
            continue
        secs = line.split()
        if len(secs) == 0:
            continue
        if len(secs) != 7:
            raise Exception("Unable to parser content %s" % line)
        id = int(secs[0])
        data[id] = {"type": int(secs[1]),
                    "point": np.array([float(secs[2]),
                                       float(secs[3]),
                                       float(secs[4]),
                                       float(secs[5])]),
                    "parent": int(secs[6]),
                    "intersected_segs": set([id, int(secs[6])])}
    return data, info

def write_swc(swc_data : dict, file_path : str, info : List[str] =[]):
    """
    Read a .swc morphology file, and return its data
    """
    file = open(file_path, "w")
    for l in info:
        file.write(l)
    for entry in sorted(swc_data):
        file.write("%i %i %f %f %f %f %i\n" % (entry,
                                               swc_data[entry]["type"],
                                               swc_data[entry]["point"][0],
                                               swc_data[entry]["point"][1],
                                               swc_data[entry]["point"][2],
                                               swc_data[entry]["point"][3],
                                               swc_data[entry]["parent"],))
    file.close()
