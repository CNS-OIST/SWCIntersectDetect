import numpy as np
import networkx as nx
from rtree import index

from .morph_io import read_swc, write_swc

import copy

SOMA = 1
AXON = 2
DEND = 3
APIC = 4

def gen_graph_and_spatial_index(swc_data):
    """
    Generate a topological graph and spatial index from swc data.
    """
    graph = nx.Graph()
    p = index.Property()
    p.dimension = 3
    spatial_index = index.Index(properties=p)

    for entry in swc_data:
        parent = swc_data[entry]['parent']
        graph.add_edge(entry, parent)
        if parent != -1:
            bounding_box = get_segment_bounding_box(swc_data[entry]["point"],
                                                    swc_data[parent]["point"])
            spatial_index.insert(entry, bounding_box)
    return graph, spatial_index


def sphere_interp(type, start, end, interp_distance):
    """
    Generate a spherical interpolation of a swc segment.
    """
    distance = np.linalg.norm(end[0:3]-start[0:3])
    n_spheres = int(distance / interp_distance)
    if n_spheres <=2 or type == SOMA:
        return [start, end]
    else:
        interpolation = np.linspace(0.0, 1.0, n_spheres)
        spheres = start + (end - start) * interpolation[:, None]
        return spheres

def sphere_intersect(s1, s2):
    """
    Check if two spheres intersect.
    """
    return np.linalg.norm(s1[0:3]-s2[0:3]) <= s1[3] + s2[3]


def segment_intersect(type0, start0, end0, type1, start1, end1, interp_distance):
    """
    Check if two swc segments intersect.
    """
    spheres_seg0 = sphere_interp(type0, start0, end0, interp_distance)
    spheres_seg1 = sphere_interp(type1, start1, end1, interp_distance)
    for s_seg0 in spheres_seg0:
        for s_seg1 in spheres_seg1:
            if sphere_intersect(s_seg0, s_seg1):
                return True
    return False

def sphere_segment_intersect(sphere, seg_type, seg_start, seg_end, interp_distance):
    """
    Check if a sphere intersects with a swc segment.
    """
    spheres_seg = sphere_interp(seg_type, seg_start, seg_end, interp_distance)
    for s_seg in spheres_seg:
        if sphere_intersect(sphere, s_seg):
            return True
    return False


def get_segment_bounding_box(start, end):
    """
    Return the bounding box of a swc segment.
    """
    start_cord = start[0:3]
    end_cord = end[0:3]
    start_min = start_cord - start[3]
    start_max = start_cord + start[3]
    end_min = end_cord - end[3]
    end_max = end_cord + end[3]
    seg_min = np.minimum(start_min, end_min)
    seg_max = np.maximum(start_max, end_max)
    bounding_box = np.append(seg_min, seg_max)
    return bounding_box

def detect_intersections(swc_data, spatial_index, interp_distance):
    """
    Detect branch intersections using spherical interpolation.
    """
    for node in swc_data:
        node_type = swc_data[node]["type"]
        node_data = swc_data[node]["point"]
        node_parent = swc_data[node]["parent"]
        if node_parent == -1:
            continue
        node_parent_data = swc_data[node_parent]["point"]
        intersected_entries = set(spatial_index.intersection(
                                  get_segment_bounding_box(node_data,
                                                           node_parent_data)))
        for query_node in intersected_entries:
            if query_node == node or query_node == node_parent:
                continue

            query_node_parent = swc_data[query_node]["parent"]
            if node == query_node_parent:
                swc_data[node]["intersected_segs"].add(query_node)
                continue
            query_node_type = swc_data[query_node]["type"]
            query_node_data = swc_data[query_node]["point"]
            query_node_parent_data = swc_data[query_node_parent]["point"]

            if segment_intersect(node_type,
                                 node_data,
                                 node_parent_data,
                                 query_node_type,
                                 query_node_data,
                                 query_node_parent_data,
                                 interp_distance):
                swc_data[node]["intersected_segs"].add(query_node)
                swc_data[query_node]["intersected_segs"].add(node)

def update_tag(swc_data, curation_data, entry, curation_tag):
    """
    Update entry tag for a swc sample entry if it hasn't been changed yet.
    """
    if swc_data[entry]["type"] == curation_data[entry]["type"]:
        curation_data[entry]["type"] = curation_tag - swc_data[entry]["type"]

def classify_intersections(swc_data, graph, interp_distance, curation_tag):
    """
    Perform intersection classification using graph theory.
    """
    curation_data = copy.deepcopy(swc_data)
    for first in swc_data:
        first_intersected = swc_data[first]["intersected_segs"]
        for second in first_intersected:
            if second == -1:
                continue
            path = nx.shortest_path(graph, source=first, target=second)
            if len(path) - 1 == 2:
                if swc_data[first]["parent"] == swc_data[second]["parent"]:
                    continue
                else:
                    mid_point = path[1]
                    mid_point_data = swc_data[mid_point]["point"]
                    first_type = swc_data[first]["type"]
                    first_data = swc_data[first]["point"]
                    first_parent = swc_data[first]["parent"]
                    first_parent_data = swc_data[first_parent]["point"]
                    second_type = second_data = swc_data[second]["type"]
                    second_data = swc_data[second]["point"]
                    second_parent = swc_data[second]["parent"]
                    second_parent_data = swc_data[second_parent]["point"]
                    if sphere_segment_intersect(mid_point_data,
                                                first_type,
                                                first_data,
                                                first_parent_data,
                                                interp_distance) or \
                       sphere_segment_intersect(mid_point_data,
                                                second_type,
                                                second_data,
                                                second_parent_data,
                                                interp_distance):
                        continue
                    else:
                        update_tag(swc_data, curation_data, first, curation_tag)
                        update_tag(swc_data, curation_data, second, curation_tag)
            else:
                second_intersected = swc_data[second]["intersected_segs"]
                for step in path:
                    if step not in first_intersected and step not in second_intersected:
                        update_tag(swc_data, curation_data, first, curation_tag)
                        update_tag(swc_data, curation_data, second, curation_tag)
                        break
    return curation_data

def run(input_swc : str, output_swc : str = None, interp_distance : float = 0.1, curation_tag : int = 999):
    """
    Perform intersection detection on the input_swc file, with interpolation distance
    interp_distance. Label the intersections in the output_swc, with curation_tag.

    If only input_swc is given as "example.swc", 
    the output file will be "example_for_curation.swc".
    """
    print("Read swc morphology from " + input_swc + "\n")
    swc_data, info = read_swc(input_swc)

    if output_swc == None:
        output_swc = input_swc[0:-4] + "_for_curation.swc"
    elif output_swc[0:-4] != ".swc":
        output_swc += ".swc"

    print("Generate topological graph and spatial index from the data\n")
    graph, spatial_index = gen_graph_and_spatial_index(swc_data)

    print("Detect intersections using spherical interpolation with minmum distance %f\n" %(interp_distance))
    detect_intersections(swc_data, spatial_index, interp_distance)

    print("Classify intersections with curation tag [%i - ORIGINAL_TAG]\n" % (curation_tag))
    curation_data = classify_intersections(swc_data, graph, interp_distance, curation_tag)

    print("Write curation data to %s\n" % (output_swc))
    info.append("# This file is for segment intersection curation use only.\n")
    info.append("# Source file: %s\n" % (input_swc))
    info.append("# interp_distance: %f\n" % (interp_distance))
    write_swc(curation_data, output_swc, info)
