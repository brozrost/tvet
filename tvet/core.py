import numpy as np
import matplotlib.pyplot as plt
import vispy
import vispy.scene
import vispy.visuals
import vispy.io
import vispy.gloo

from . import scattering
from .io import load_obj
from . import ephemerides as ephems
from . import _tvet

class Asteroid(object):
    def __init__(self, args=None, filename=None):
        self.args = args
        self.filename = filename

        self.vertices, self.faces = load_obj(self.filename)
        self.size = np.max(self.vertices) - np.min(self.vertices)
        self.vertices *= 1.9 / self.size

        # Select scattering function based on CLI argument
        if self.args and hasattr(self.args, "scattering"):
            if self.args.scattering == "lambert":
                self.f_func = scattering.f_lambert
            elif self.args.scattering == "lommel":
                self.f_func = scattering.f_lommel
            elif self.args.scattering == "hapke":
                self.f_func = scattering.f_hapke
            else:
                self.f_func = scattering.f_lambert
        else:
            self.f_func = scattering.f_lambert

    def get_geometry(self):
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

    def get_cosines(self, s=(1, 0, 0), o=(0, 0, 1)):
        self.get_geometry()

        self.s = np.array(s)
        self.o = np.array(o)
        self.alpha = np.arccos(np.dot(s, o))

        mu_i = []
        mu_e = []

        for normal in self.normals:
            mu_i.append(np.dot(s, normal))
            mu_e.append(np.dot(o, normal))

        self.mu_i = np.array(mu_i)
        self.mu_e = np.array(mu_e)

        self.mu_i = np.where(self.mu_i > 0.0, self.mu_i, 0.0)
        self.mu_e = np.where(self.mu_e > 0.0, self.mu_e, 0.0)

        self.nu_i = np.zeros((len(self.faces)))
        self.nu_e = np.zeros((len(self.faces)))

        mu_i_C = np.ascontiguousarray(self.mu_i, dtype=np.double)
        mu_e_C = np.ascontiguousarray(self.mu_e, dtype=np.double)
        nu_i_C = np.ascontiguousarray(self.nu_i, dtype=np.double)
        nu_e_C = np.ascontiguousarray(self.nu_e, dtype=np.double)

        faces_C = np.ascontiguousarray(self.faces, dtype=np.intc)
        vertices_C = np.ascontiguousarray(self.vertices, dtype=np.double)
        normals_C = np.ascontiguousarray(self.normals, dtype=np.double)
        centers_C = np.ascontiguousarray(self.centers, dtype=np.double)

        s_C = np.ascontiguousarray(self.s, dtype=np.double)
        o_C = np.ascontiguousarray(self.o, dtype=np.double)

        nof_faces = len(self.faces)
        nof_nodes = len(self.vertices)

        _tvet.non(mu_i_C, mu_e_C, nof_faces, nu_i_C, nu_e_C)
        _tvet.nu(faces_C, nof_faces, vertices_C, nof_nodes, normals_C, centers_C, s_C, nu_i_C)
        _tvet.nu(faces_C, nof_faces, vertices_C, nof_nodes, normals_C, centers_C, o_C, nu_e_C)

        self.nu_i = nu_i_C
        self.nu_e = nu_e_C

    def get_fluxes(self, s=None, o=None):
        self.get_geometry()
        s = s if s is not None else (self.args.s if self.args and hasattr(self.args, "s") else (1, 0, 0))
        o = o if o is not None else (self.args.o if self.args and hasattr(self.args, "o") else (0, 0, 1))
        self.get_cosines(s=s, o=o)

        phi_s = 1361. # W/m^2
        self.phi_i = phi_s * self.mu_i * self.nu_i

        f = []
        A_w = 0.23
        self.f_L = A_w/(4.0*np.pi)
        for i in range(len(self.mu_e)):
            f.append(self.f_func(self.f_L, self.mu_i[i], self.mu_e[i], self.alpha))
        self.f = np.array(f)

        self.I = self.f * self.phi_i
        self.phi_e = self.I * self.mu_e * self.nu_e

        self.total = np.sum(self.phi_e)

    def get_ephems(
        *,
        body: str,
        start_time: str,
        stop_time: str,
        step_size: str,
        observer_center: str = "500@399",
        sun_center: str = "500@10",
        normalize: bool = True,
        timeout: float = 30.0
    ):
        o_xyz = ephems.fetch_ephems(
            body=body,
            center=observer_center,
            start_time=start_time,
            stop_time=stop_time,
            step_size=step_size,
            timeout=timeout
        )

        s_xyz = ephems.fetch_ephems(
            body=body,
            center=sun_center,
            start_time=start_time,
            stop_time=stop_time,
            step_size=step_size,
            timeout=timeout
        )

        o = np.asarray(o_xyz, dtype=np.float64)
        s = np.asarray(s_xyz, dtype=np.float64)

        if o.shape != s.shape:
            raise ephems.HorizonsError(f"Ephemerides shape mismatch: o={o.shape}, s={s.shape}")

        if not normalize:
            return s, o

        o_norm = np.linalg.norm(o, axis=1, keepdims=True)
        s_norm = np.linalg.norm(s, axis=1, keepdims=True)

        o_norm = np.where(o_norm > 0.0, o_norm, 1.0)
        s_norm = np.where(s_norm > 0.0, s_norm, 1.0)

        o_unit = o / o_norm
        s_unit = s / s_norm

        return s_unit, o_unit


    def light_curve(self, n=100):
        s = self.s if hasattr(self, "s") else np.array(self.args.s if self.args and hasattr(self.args, "s") else (1, 0, 0))
        x, y, z = s
        total = []

        for i in range(n+1):
            gamma = 0 + 2.0*np.pi * i/n

            x_ = x * np.cos(gamma) + y * np.sin(gamma)
            y_ = -x * np.sin(gamma) + y * np.cos(gamma)
            z_ = z

            s_ = np.array([x_, y_, z_])

            self.get_fluxes(s=s_, o=self.o if hasattr(self, "o") else None)
            total.append((gamma, self.total))

        total = np.array(total)

        # Normalize
        total[:, 0] = (total[:, 0] - np.min(total[:, 0])) / (np.max(total[:, 0]) - np.min(total[:, 0]))

        den = np.max(total[:, 1]) - np.min(total[:, 1])
        if den > 0.0:
            total[:, 1] = -(total[:, 1] - np.min(total[:, 1])) / den
        else:
            total[:, 1] = 0.0

        return total

    def plot_light_curve(self, curve_points):
        if curve_points is None:
            curve_points = self.light_curve()

        plt.figure(figsize=(8, 4))
        plt.plot(curve_points[:, 0], curve_points[:, 1], color='blue')
        plt.xlabel('Normalized phase angle')
        plt.ylabel('Normalized flux')
        plt.title('Asteroid Light Curve')
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    def interactive_plot_light_curve(self):
        curve_points = self.light_curve()

        light_curve = vispy.scene.visuals.Line(pos=(curve_points*200) + (40, 560), color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Line(pos=((20, 570), (250, 570)), color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Line(pos=((30, 360), (30, 580)), color='white', parent=self.canvas.scene)

    def interactive_plot(self):
        # Call geometry and flux setup with CLI vectors
        self.get_geometry()
        s = self.args.s if self.args and hasattr(self.args, "s") else (1, 0, 0)
        o = self.args.o if self.args and hasattr(self.args, "o") else (0, 0, 1)
        self.get_cosines(s=s, o=o)
        self.get_fluxes()

        # Provide defaults if args is None
        shininess = self.args.shininess if self.args and hasattr(self.args, "shininess") else 100
        wireframe_width = self.args.wireframe_width if self.args and hasattr(self.args, "wireframe_width") else 1

        self.canvas = vispy.scene.SceneCanvas(keys='interactive')
        self.canvas.size = 1920, 1080
        self.view = self.canvas.central_widget.add_view()

        mesh = vispy.scene.visuals.Mesh(self.vertices, self.faces, color='gray')
        mesh.transform = vispy.scene.transforms.MatrixTransform()
        self.view.add(mesh)

        pos = np.array([self.centers, self.centers + 0.1 * self.normals])
        connect = []
        n = len(self.centers)
        for i in range(n):
            connect.append(np.array([i, n + i]))

        normals = vispy.scene.visuals.Line(pos=pos, connect=np.array(connect), color='white', parent=self.view.scene)
        normals.visible = False
        # text = vispy.scene.visuals.Text(str(i), pos=self.centers[i], font_size=10, color='white')
        # self.view.add(text)

        s = vispy.scene.visuals.Line(pos=np.array([(0, 0.02, 0), self.s+ (0, 0.02, 0)]), color='yellow', parent=self.view.scene)
        o = vispy.scene.visuals.Line(pos=np.array([(0, 0.02, 0), self.o + (0, 0.02, 0)]), color='magenta', parent=self.view.scene)
        
        vispy.scene.visuals.Text("'1' to show phi_i", anchor_x='left', pos=(20, 20), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("'2' to show phi_e", anchor_x='left', pos=(20, 40), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("'3' to show the wireframe", anchor_x='left', pos=(20, 60), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("'4' to show normals", anchor_x='left', pos=(20, 80), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("'5' to show flat model", anchor_x='left', pos=(20, 100), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("'6' to show smooth model", anchor_x='left', pos=(20, 120), font_size=10,
                            color='white', parent=self.canvas.scene)
        
        vispy.scene.visuals.Text("'a' to use Lambert", anchor_x='left', pos=(20, 160), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("'b' to use Lommel", anchor_x='left', pos=(20, 180), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("'c' to use Hapke", anchor_x='left', pos=(20, 200), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("'l' to create a light curve", anchor_x='left', pos=(20, 220), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("'s' to screenshot", anchor_x='left', pos=(20, 240), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("'q' to quit", anchor_x='left', pos=(20, 260), font_size=10,
                            color='white', parent=self.canvas.scene)
        
        vispy.scene.visuals.Text("Number of vertices: %d" %(len(self.vertices)), anchor_x='left', pos=(20, 300), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("Number of faces: %d" %(len(self.faces)), anchor_x='left', pos=(20, 320), font_size=10,
                            color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Text("Asteroid size: %f" %(self.size), anchor_x='left', pos=(20, 340), font_size=10,
                            color='white', parent=self.canvas.scene)
        
        vispy.scene.visuals.XYZAxis(parent=self.view.scene)

        shading_filter = vispy.visuals.filters.ShadingFilter(
            shading='smooth',
            shininess=shininess,
            ambient_coefficient=0.0,
            diffuse_coefficient=1.0,
            specular_coefficient=0.0,
            ambient_light='white',
            diffuse_light='white',
            specular_light='white',
        )
        mesh.attach(shading_filter)

        wireframe_filter = vispy.visuals.filters.WireframeFilter(
            width=wireframe_width,
            color='green',
            wireframe_only=False,
            faces_only=True,
            enabled=False,
        )
        mesh.attach(wireframe_filter)

        self.view.camera = vispy.scene.cameras.TurntableCamera(center=(0, 0, 0))
        self.view.camera.depth_value = 1e3

        light_dir = -self.s
        shading_filter.light_dir = light_dir[:3]

        def plot_fluxes(phi):
            shading_filter.shading = None
            wireframe_filter.enabled = True
            wireframe_filter.wireframe_only = False
            wireframe_filter.faces_only = False
            normals.visible = False
            color = np.array([0.5, 0.5, 0.5]) / np.percentile(phi, 99)
            face_colors = []
            for face in phi:
                face_colors.append(face * color)
            mesh.set_data(self.vertices, self.faces, face_colors=face_colors)
            mesh.update()

        @self.canvas.events.key_press.connect
        def on_key_press(event):
            if event.key in ['q', 'Q']:
                vispy.app.quit()

            elif event.key == 's':
                vispy.io.write_png("output/vispy_screenshot.png", vispy.gloo.util._screenshot())

            elif event.key == '1':
                plot_fluxes(phi=self.phi_i)
                wireframe_filter.enabled = False

            elif event.key == '2':
                plot_fluxes(phi=self.phi_e)

            elif event.key == '3':
                shading_filter.shading = 'flat'
                wireframe_filter.enabled = True
                wireframe_filter.wireframe_only = True
                wireframe_filter.faces_only = False
                normals.visible = False
                mesh.update()

            elif event.key == '4':
                shading_filter.shading = 'flat'
                wireframe_filter.enabled = True
                wireframe_filter.wireframe_only = False
                wireframe_filter.faces_only = False
                normals.visible = True
                face_colors = []
                for i in range(len(self.faces)):
                    face_colors.append(np.array([0.6, 0.6, 0.6]))
                mesh.set_data(self.vertices, self.faces, face_colors=face_colors)
                mesh.update()

            elif event.key == '5':
                shading_filter.shading = 'flat'
                wireframe_filter.enabled = False
                normals.visible = False
                face_colors = []
                for i in range(len(self.faces)):
                    face_colors.append(np.array([0.6, 0.6, 0.6]))
                mesh.set_data(self.vertices, self.faces, face_colors=face_colors)
                mesh.update()

            elif event.key == '6':
                shading_filter.shading = 'smooth'
                wireframe_filter.enabled = False
                normals.visible = False
                face_colors = []
                for i in range(len(self.faces)):
                    face_colors.append(np.array([0.6, 0.6, 0.6]))
                mesh.set_data(self.vertices, self.faces, face_colors=face_colors)
                mesh.update()

            elif event.key == 'a':
                self.f_func = scattering.f_lambert
                self.get_fluxes()

            elif event.key == 'b':
                self.f_func = scattering.f_lommel
                self.get_fluxes()

            elif event.key == 'c':
                scattering.B0 = 1.32
                scattering.minh = 0.20
                scattering.ming = -0.35
                scattering.bartheta = 10.0 * np.pi/180

                scattering.init_hapke(self.alpha)
                self.f_func = scattering.f_hapke
                self.get_fluxes()

            elif event.key == 'l':
                curve_points = self.light_curve()
                self.plot_light_curve(curve_points)

        self.canvas.show()