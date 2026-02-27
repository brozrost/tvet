#!/usr/bin/env python3

import numpy as np

def load_ele_file(filename, debug=False): 
    '''Load an ELE file.'''
    
    ele = []

    f = open(filename, "r")

    for line in f.readlines()[1:]:
        if len(line) == 0:
            continue
        elif line[0] == '#': 
            continue

        l = line.split()

        ele.append(list(map(int, l[1:])))

    f.close

    if debug:
        print("number of elements = " + str(len(ele)) + ", beginning with element " + str(ele[0]))

    return ele

def load_node_file(filename, debug=False): 
    '''Load a NODE file.'''
    
    node = []

    f = open(filename, "r")

    for line in f.readlines()[1:]:
        if len(line) == 0:
            continue
        elif line[0] == '#': 
            continue

        l = line.split()

        node.append(list(map(float, l[1:])))

    f.close()

    if debug:
        print("number of nodes = " + str(len(node)) + ", beginning with node " + str(node[0]))

    return node

def load_face_file(filename, debug=False): 
    '''Load a FACE file.'''

    face = []
    
    f = open(filename, "r")

    for line in f.readlines()[1:]:
        if len(line) == 0:
            continue
        elif line[0] == '#': 
            continue

        l = line.split()

        face.append(list(map(int, l[1:4])))

    f.close()

    if debug:
       print("number of faces = " + str(len(face)) + ", beginning with face " + str(face[0]))

    return face

def load_obj_file(filename, debug=False):
    '''Load a Wavefront OBJ file.''' 
    
    node = []
    face = []

    f = open(filename, "r")

    for line in f.readlines(): 
        if len(line) == 0:
            continue
        elif line[0] == '#': 
            continue

        l = line.split()
        
        if l[0] == 'v':
            node.append(list(map(float, l[1:])))
        elif l[0] == 'f':
            face.append(list(map(lambda x: int(x) - 1, l[1:])))

    f.close()

    if debug:
        print("number of nodes = " + str(len(node)) + ", beginning with node " + str(node[0]))
        print("number of faces = " + str(len(face)) + ", beginning with face " + str(face[0]))

    return np.array(node), np.array(face)

def check_filetype(filename):
    for index, character in enumerate(filename):
        if character == ".":
            if filename[index+1:] == "obj":
                return True
            else:
                return False

def main():
    ele = load_ele_file("src/sample_files/tri_file_octdecv_1.1.ele", True)
    node = load_node_file("src/sample_files/tri_file_octdecv_1.1.node", True)
    face = load_face_file("src/sample_files/tri_file_octdecv_1.1.face", True)
    obj_node, obj_face = load_obj_file("src/sample_files/tri_file_octdecv_1.obj", True)

if __name__ == "__main__": 
    main()
