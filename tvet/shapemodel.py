import numpy as np

from . import io

class ShapeModel:
    """
    Store a triangular mesh together with derived per-face geometric data.

    The class manages:
    - raw mesh data (vertices and triangular faces),
    - basic mesh metadata such as size and element counts,
    - derived geometric quantities such as triangle centers and unit normals,
    - contiguous NumPy buffers prepared for efficient downstream low-level use.

    Geometry is computed lazily via ``compute_geometry()`` and cached until the
    mesh is replaced.
    """

    def __init__(self):
        """
        Initializes an empty shape model with no loaded mesh data.
        """

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
        """
        Store mesh data and normalize its scale for visualization.

        The input arrays are converted to the expected dtypes:
        - vertices -> np.double
        - faces    -> np.intc

        After loading, the mesh is uniformly scaled so that its overall extent
        remains visually manageable when rendered together with axes and
        illumination/observer vectors.

        Any previously computed geometry is invalidated because it no longer
        corresponds to the new mesh.

        Args:
            vertices:
                Array-like object of vertex coordinates with shape (N, 3).
            faces:
                Array-like object of triangle indices with shape (M, 3).

        Raises:
            ValueError:
                If the mesh has zero extent and therefore cannot be scaled.
        """

        self.vertices = np.asarray(vertices, dtype=np.double)
        self.faces = np.asarray(faces, dtype=np.intc)

        self.size = np.max(self.vertices) - np.min(self.vertices)
        if self.size == 0.0:
            raise ValueError("Invalid mesh: size is zero.")

        # Invalidate all geometry derived from the previous mesh.
        self.centers = None
        self.normals = None

        self.faces_C = None
        self.vertices_C = None
        self.normals_C = None
        self.centers_C = None

        self.nof_faces = 0
        self.nof_nodes = 0

        self._geometry_ready = False

    def load_obj(self, filename: str):
        """
        Load a mesh from an OBJ file and store it in this shape model.

        The parsed mesh is passed through ``set_mesh()`` so dtype conversion,
        scaling, and cache invalidation all happen consistently.

        Args:
            filename:
                Path to the OBJ file.
        """
                
        vertices, faces = io.load_obj_file(filename)
        self.set_mesh(vertices, faces)
        
    def compute_geometry(self):
        """
        Compute per-face centers, unit normals, and contiguous working arrays.

        For each triangular face, this method computes:
        - the triangle center as the mean of its three vertices,
        - the face normal as the normalized cross product of two edge vectors.

        It then prepares contiguous arrays suitable for efficient downstream use,
        including calls into lower-level native code.

        The result is cached. If geometry has already been computed for the
        current mesh, the method returns immediately.

        Raises:
            ValueError:
                If no mesh has been loaded before geometry computation is requested.
        """

        if self.vertices is None or self.faces is None:
            raise ValueError("Geometry not available: no OBJ mesh loaded.")

        if self._geometry_ready:
            return

        self.centers = []
        self.normals = []

        for face in self.faces:
            A = self.vertices[face[0]]
            B = self.vertices[face[1]]
            C = self.vertices[face[2]]

            # Triangle centroid.
            T = 1/3 * (A + B + C)

            # Construct two edge vectors and compute the face normal from their
            # cross product. The normal is then normalized to unit length.
            a = B - C
            b = C - A
            n = np.cross(a, b)
            n /= np.sqrt(np.dot(n, n))

            self.centers.append(T)
            self.normals.append(n)

        self.centers = np.array(self.centers)
        self.normals = np.array(self.normals)

        # Create contiguous arrays with stable dtypes for efficient interop with
        # lower-level code that expects raw contiguous memory.
        self.faces_C = np.ascontiguousarray(self.faces, dtype=np.intc)
        self.vertices_C = np.ascontiguousarray(self.vertices, dtype=np.double)
        self.centers_C = np.ascontiguousarray(self.centers, dtype=np.double)
        self.normals_C = np.ascontiguousarray(self.normals, dtype=np.double)

        # Cache mesh sizes used repeatedly by downstream routines.
        self.nof_faces = int(self.faces_C.shape[0])
        self.nof_nodes = int(self.vertices_C.shape[0])

        self._geometry_ready = True
