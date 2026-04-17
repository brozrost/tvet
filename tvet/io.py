import numpy as np

# MARK: - load_ele_file()

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

# MARK: - load_node_file()

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

# MARK: - load_face_file()

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

# MARK: - load_obj_file()

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

# MARK: - save_obj_file()

def save_obj_file(*, path: str, vertices, faces):
    with open(path, "w", encoding="utf-8") as f:
        for x, y, z in vertices:
            f.write(f"v {x} {y} {z}\n")

        for i, j, k in faces:
            f.write(f"f {i + 1} {j + 1} {k + 1}\n")

# MARK: - save_spin()

def save_spin(*, path: str, l, b, period, epoch, phi0):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{l} {b} {period}\n{epoch} {phi0}")

# MARK: - check_filetype()

def check_filetype(filename):
    for index, character in enumerate(filename):
        if character == ".":
            if filename[index+1:] == "obj":
                return True
            else:
                return False
