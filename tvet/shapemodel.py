import numpy as np

from . import io

class ShapeModel:
    def __init__(self):
        self.vertices = None
        self.faces = None
        self.size = None

        self.centers = None
        self.normals = None

        self.faces_C = None
        self.vertices_C = None
        self.normals_C = None
        self.centers_C = None

        self.nof_faces = 0
        self.nof_nodes = 0

        # Guard against recomputing geometry
        self._geometry_ready = False

    def set_mesh(self, vertices, faces):
        self.vertices = np.asarray(vertices, dtype=np.double)
        self.faces = np.asarray(faces, dtype=np.intc)

        self.size = np.max(self.vertices) - np.min(self.vertices)
        if self.size == 0.0:
            raise ValueError("Invalid mesh: size is zero.")
        
        # This shrinks the model down so the axis and s/o vectors are visible
        self.vertices *= 1.9 / self.size

    def load_obj(self, filename: str):
        vertices, faces = io.load_obj_file(filename)
        self.set_mesh(vertices, faces)
        
    def compute_geometry(self):
        if self.vertices is None or self.faces is None:
            raise ValueError("Geometry not available: no OBJ mesh loaded.")

        if getattr(self, "_geometry_ready", False):
            return

        self.centers = []
        self.normals = []

        for face in self.faces:
            A = self.vertices[face[0]]
            B = self.vertices[face[1]]
            C = self.vertices[face[2]]

            T = 1/3 * (A + B + C)
            a = B - C
            b = C - A
            n = np.cross(a, b)
            n /= np.sqrt(np.dot(n, n))

            self.centers.append(T)
            self.normals.append(n)

        self.centers = np.array(self.centers)
        self.normals = np.array(self.normals)

        self.faces_C = np.ascontiguousarray(self.faces, dtype=np.intc)
        self.vertices_C = np.ascontiguousarray(self.vertices, dtype=np.double)
        self.centers_C = np.ascontiguousarray(self.centers, dtype=np.double)
        self.normals_C = np.ascontiguousarray(self.normals, dtype=np.double)

        self.nof_faces = int(self.faces_C.shape[0])
        self.nof_nodes = int(self.vertices_C.shape[0])

        self._geometry_ready = True
