import numpy as np
import matplotlib.pyplot as plt
import vispy
import vispy.scene
import vispy.visuals
import vispy.io
import vispy.gloo

from . import shapemodel
from . import scattering
from . import horizons
from . import damit
from . import lightcurve
from . import _tvet

class Asteroid:
    def __init__(self, args=None, filename=None):
        self.args = args
        self.filename = filename

        self.s = np.array([1.0, 0.0, 0.0], dtype=np.double)
        self.o = np.array([0.0, 0.0, 1.0], dtype=np.double)

        self.s_array = None
        self.o_array = None

        self.nu_i_C = None
        self.nu_e_C = None

        self.l = 0.0
        self.b = np.pi / 2
        self.period = 1.0
        self.epoch = 0.0
        self.phi0 = 0.0

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

        self.shape = shapemodel.ShapeModel()
        if self.filename is not None:
            self.shape.load_obj(self.filename)
        self.horizons = horizons.HorizonsClient()
        self.damit = damit.DamitClient()
        self.light_curve = lightcurve.LightCurve(self)

    def _rotate_x(self, a, phi):
        x, y, z = a

        x_ = x
        y_ = y * np.cos(phi) - z * np.sin(phi)
        z_ = y * np.sin(phi) + z * np.cos(phi)

        a_ = np.array([x_, y_, z_], dtype=np.double)

        return a_

    def _rotate_y(self, a, phi):
        x, y, z = a

        x_ = x * np.cos(phi) + z * np.sin(phi)
        y_ = y
        z_ = - x * np.sin(phi) + z * np.cos(phi)

        a_ = np.array([x_, y_, z_], dtype=np.double)

        return a_

    def _rotate_z(self, a, phi):
        x, y, z = a

        x_ = x * np.cos(phi) - y * np.sin(phi)
        y_ = x * np.sin(phi) + y * np.cos(phi)
        z_ = z

        a_ = np.array([x_, y_, z_], dtype=np.double)

        return a_

    def _match_vector(self, a, start_time):
        phi1 = 2 * np.pi * (start_time - self.epoch) / self.period + self.phi0
        phi2 = np.pi / 2 - self.b
        phi3 = self.l

        # match damits ecliptic coordinates
        eps = (23.0 + 26.0 / 60.0 + (21.406 / 3600.0)) * np.pi / 180
        a_ = self._rotate_x(a, -eps)

        a_ = self._rotate_z(a_, -phi3)
        a_ = self._rotate_y(a_, -phi2)
        a_ = self._rotate_z(a_, -phi1)

        return a_

    def set_body_frame(self, start_time):
        self.s = self._match_vector(self.s, start_time)
        self.o = self._match_vector(self.o, start_time)

        if self.s_array is not None:
            for vector in self.s_array:
                self._match_vector(vector, start_time)

        if self.o_array is not None:
            for vector in self.o_array:
                self._match_vector(vector, start_time)

    def get_geometry(self):
        if self.shape.vertices is None or self.shape.faces is None:
            raise ValueError("Geometry not available: no mesh loaded.")
        self.shape.compute_geometry()

    def get_cosines(self, s=None, o=None):
        self.shape.compute_geometry()

        if s is None: 
            s = self.s
        if o is None: 
            o = self.o

        s = np.asarray(s, dtype=np.double)
        o = np.asarray(o, dtype=np.double)

        d = np.dot(s, o)
        d = np.clip(d, -1.0, 1.0)
        self.alpha = np.arccos(d)

        self.mu_i = self.shape.normals @ s
        self.mu_e = self.shape.normals @ o

        self.mu_i = np.maximum(self.mu_i, 0.0)
        self.mu_e = np.maximum(self.mu_e, 0.0)

        mu_i_C = np.asarray(self.mu_i, dtype=np.double, order="C")
        mu_e_C = np.asarray(self.mu_e, dtype=np.double, order="C")

        if self.nu_i_C is None or self.nu_i_C.shape[0] != self.shape.nof_faces:
            self.nu_i_C = np.empty(self.shape.nof_faces, dtype=np.double)
            self.nu_e_C = np.empty(self.shape.nof_faces, dtype=np.double)

        self.nu_i_C.fill(0.0)
        self.nu_e_C.fill(0.0)

        _tvet.non(mu_i_C, mu_e_C, self.shape.nof_faces, self.nu_i_C, self.nu_e_C)

        _tvet.nu(
            self.shape.faces_C, self.shape.nof_faces, 
            self.shape.vertices_C, self.shape.nof_nodes, 
            self.shape.normals_C, self.shape.centers_C, 
            s, self.nu_i_C
        )

        _tvet.nu(
            self.shape.faces_C, self.shape.nof_faces, 
            self.shape.vertices_C, self.shape.nof_nodes, 
            self.shape.normals_C, self.shape.centers_C, 
            o, self.nu_e_C
        )

        self.nu_i = self.nu_i_C
        self.nu_e = self.nu_e_C

    def get_fluxes(self, s=None, o=None, f_func=None):
        if s is None: 
            s = self.s
        if o is None: 
            o = self.o
        if f_func is None:
            f_func = self.f_func

        self.get_cosines(s=s, o=o)

        phi_s = 1361. # W/m^2
        self.phi_i = phi_s * self.mu_i * self.nu_i

        f = []
        A_w = 0.23
        self.f_L = A_w/(4.0*np.pi)
        for i in range(len(self.mu_e)):
            f.append(f_func(self.f_L, self.mu_i[i], self.mu_e[i], self.alpha))
        self.f = np.array(f)

        self.I = self.f * self.phi_i
        self.phi_e = self.I * self.mu_e * self.nu_e

        self.total = np.sum(self.phi_e)

    def get_single_ephem(
        self,
        *,
        body: str,
        epoch: str,
        observer_center: str = "500@399",
        sun_center: str = "500@10",
        normalize: bool = True,
        timeout: float = 30.0
    ):
        return self.horizons.fetch_single_so(
            body=body,
            epoch=epoch,
            observer_center=observer_center,
            sun_center=sun_center,
            normalize=normalize,
            timeout=timeout
        )

    def get_ephems(
        self,
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
        return self.horizons.fetch_so(
            body=body,
            start_time=start_time,
            stop_time=stop_time,
            step_size=step_size,
            observer_center=observer_center,
            sun_center=sun_center,
            normalize=normalize,
            timeout=timeout
        )
    
    def get_damit(
        self,
        *,
        model_id: str,
        timeout: float = 30.0
    ):
        vertices, faces = self.damit.fetch_obj(
            model_id=model_id, 
            timeout=timeout
        )
        self.shape.set_mesh(vertices, faces)

        self.l, self.b, self.period, self.epoch, self.phi0 = self.damit.fetch_spin(
            model_id=model_id, 
            timeout=timeout
        )
    
    def get_light_curve(self):
        pass

    def get_light_curve_for_period(self, s=None, o=None, n=100, start=None, period=None, epoch=None, l=None, b=None, phi0=None):
        return self.light_curve.compute_for_period()

    def plot_light_curve(self, curve_points):
        if curve_points is None:
            curve_points = self.get_light_curve_for_period()

        plt.figure(figsize=(8, 4))
        plt.plot(curve_points[:, 0], curve_points[:, 1], marker='+', color='red')
        plt.xlabel('Normalized phase angle')
        plt.ylabel('Normalized flux')
        plt.title('Asteroid Light Curve')
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    def interactive_plot_light_curve(self):
        curve_points = self.get_light_curve()

        light_curve = vispy.scene.visuals.Line(pos=(curve_points*200) + (40, 560), color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Line(pos=((20, 570), (250, 570)), color='white', parent=self.canvas.scene)
        vispy.scene.visuals.Line(pos=((30, 360), (30, 580)), color='white', parent=self.canvas.scene)

    # MARK: - interactive_plot()

    def interactive_plot(self):
        # Call geometry and flux setup with CLI vectors
        self.get_geometry()
        self.get_cosines()
        self.get_fluxes()

        # Provide defaults if args is None
        shininess = self.args.shininess if self.args and hasattr(self.args, "shininess") else 100
        wireframe_width = self.args.wireframe_width if self.args and hasattr(self.args, "wireframe_width") else 1

        self.canvas = vispy.scene.SceneCanvas(keys='interactive')
        self.canvas.size = 1920, 1080
        self.view = self.canvas.central_widget.add_view()

        mesh = vispy.scene.visuals.Mesh(self.shape.vertices, self.shape.faces, color='gray')
        mesh.transform = vispy.scene.transforms.MatrixTransform()
        self.view.add(mesh)

        pos = np.array([self.shape.centers, self.shape.centers + self.shape.normals])
        connect = []
        n = len(self.shape.centers)
        for i in range(n):
            connect.append(np.array([i, n + i]))

        normals = vispy.scene.visuals.Line(pos=pos, connect=np.array(connect), color='white', parent=self.view.scene)
        normals.visible = False
        # text = vispy.scene.visuals.Text(str(i), pos=self.centers[i], font_size=10, color='white')
        # self.view.add(text)

        self.overlays = []

        s_line = vispy.scene.visuals.Line(pos=np.array([(0, 0, 0), self.shape.size * self.s / 1.5]), color='yellow', parent=self.view.scene)
        o_line = vispy.scene.visuals.Line(pos=np.array([(0, 0, 0), self.shape.size * self.o / 1.5]), color='magenta', parent=self.view.scene)

        t1 = vispy.scene.visuals.Text("'1' to show phi_i", anchor_x='left', pos=(20, 20), font_size=10,
                            color='white', parent=self.canvas.scene)
        t2 = vispy.scene.visuals.Text("'2' to show phi_e", anchor_x='left', pos=(20, 40), font_size=10,
                            color='white', parent=self.canvas.scene)
        t3 = vispy.scene.visuals.Text("'3' to show the wireframe", anchor_x='left', pos=(20, 60), font_size=10,
                            color='white', parent=self.canvas.scene)
        t4 = vispy.scene.visuals.Text("'4' to show normals", anchor_x='left', pos=(20, 80), font_size=10,
                            color='white', parent=self.canvas.scene)
        t5 = vispy.scene.visuals.Text("'5' to show flat model", anchor_x='left', pos=(20, 100), font_size=10,
                            color='white', parent=self.canvas.scene)
        t6 = vispy.scene.visuals.Text("'6' to show smooth model", anchor_x='left', pos=(20, 120), font_size=10,
                            color='white', parent=self.canvas.scene)
        
        ta = vispy.scene.visuals.Text("'z' to use Lambert", anchor_x='left', pos=(20, 160), font_size=10,
                            color='white', parent=self.canvas.scene)
        tb = vispy.scene.visuals.Text("'x' to use Lommel", anchor_x='left', pos=(20, 180), font_size=10,
                            color='white', parent=self.canvas.scene)
        tc = vispy.scene.visuals.Text("'c' to use Hapke", anchor_x='left', pos=(20, 200), font_size=10,
                            color='white', parent=self.canvas.scene)
        ts = vispy.scene.visuals.Text("'p' to screenshot", anchor_x='left', pos=(20, 240), font_size=10,
                            color='white', parent=self.canvas.scene)
        th = vispy.scene.visuals.Text("'h' to toggle overlays", anchor_x='left', pos=(20, 260), font_size=10,
                            color='white', parent=self.canvas.scene)
        tq = vispy.scene.visuals.Text("'q' to quit", anchor_x='left', pos=(20, 280), font_size=10,
                            color='white', parent=self.canvas.scene)
        
        tver = vispy.scene.visuals.Text(f"Number of vertices: {len(self.shape.vertices)}", anchor_x='left', pos=(20, 320), font_size=10,
                            color='white', parent=self.canvas.scene)
        tfac = vispy.scene.visuals.Text(f"Number of faces: {len(self.shape.faces)}", anchor_x='left', pos=(20, 340), font_size=10,
                            color='white', parent=self.canvas.scene)
        tsiz = vispy.scene.visuals.Text(f"Asteroid size: {self.shape.size}", anchor_x='left', pos=(20, 360), font_size=10,
                            color='white', parent=self.canvas.scene)
        
        scale = self.shape.size / 1.5
        axis = vispy.scene.visuals.XYZAxis(parent=self.view.scene)
        axis.transform = vispy.visuals.transforms.STTransform(scale=(scale, scale, scale))

        self.overlays.append(s_line)
        self.overlays.append(o_line)
        self.overlays.append(t1)
        self.overlays.append(t2)
        self.overlays.append(t3)
        self.overlays.append(t4)
        self.overlays.append(t5)
        self.overlays.append(t6)
        self.overlays.append(ta)
        self.overlays.append(tb)
        self.overlays.append(tc)
        self.overlays.append(ts)
        self.overlays.append(th)
        self.overlays.append(tq)
        self.overlays.append(tver)
        self.overlays.append(tfac)
        self.overlays.append(tsiz)
        self.overlays.append(axis)

        shading_filter = vispy.visuals.filters.ShadingFilter(
            shading='flat',
            shininess=shininess,
            ambient_coefficient=0.0,
            diffuse_coefficient=1.0,
            specular_coefficient=0.0,
            ambient_light='white',
            diffuse_light='white',
            specular_light='white',
        )
        mesh.attach(shading_filter)
        self._previous_shading = 'flat'

        wireframe_filter = vispy.visuals.filters.WireframeFilter(
            width=wireframe_width,
            color='green',
            wireframe_only=False,
            faces_only=True,
            enabled=False,
        )
        mesh.attach(wireframe_filter)
        
        ox, oy, oz = self.o
        r = np.linalg.norm(self.o)
        azimuth = 180 - np.degrees(np.arctan2(ox, oy))
        elevation = np.degrees(np.arcsin(oz / r))
        self.view.camera = vispy.scene.cameras.TurntableCamera(
            center=(0, 0, 0),
            azimuth=azimuth,
            elevation=elevation,
            fov=0
        )
        self.view.camera.depth_value = 1e3

        light_dir = -self.s
        shading_filter.light_dir = light_dir[:3]

        def plot_fluxes(phi):
            self._previous_shading = shading_filter.shading
            shading_filter.shading = None
            wireframe_filter.enabled = False
            wireframe_filter.wireframe_only = False
            wireframe_filter.faces_only = False
            normals.visible = False
            color = np.array([0.5, 0.5, 0.5]) / np.percentile(phi, 99)
            face_colors = []
            for face in phi:
                face_colors.append(face * color)
            mesh.set_data(self.shape.vertices, self.shape.faces, face_colors=face_colors)
            mesh.update()

        @self.canvas.events.key_press.connect
        def on_key_press(event):
            if event.key in ['q', 'Q']:
                vispy.app.quit()

            elif event.key == 'p':
                vispy.io.write_png("out/vispy_screenshot.png", vispy.gloo.util._screenshot())

            elif event.key == '1':
                plot_fluxes(phi=self.phi_i)

            elif event.key == '2':
                plot_fluxes(phi=self.phi_e)

            elif event.key == '3':
                shading_filter.shading = self._previous_shading
                wireframe_filter.enabled = True
                wireframe_filter.wireframe_only = True
                wireframe_filter.faces_only = False
                normals.visible = False
                mesh.update()

            elif event.key == '4':
                shading_filter.shading = self._previous_shading
                wireframe_filter.enabled = True
                wireframe_filter.wireframe_only = False
                wireframe_filter.faces_only = False
                normals.visible = True
                face_colors = []
                for i in range(len(self.shape.faces)):
                    face_colors.append(np.array([0.6, 0.6, 0.6]))
                mesh.set_data(self.shape.vertices, self.shape.faces, face_colors=face_colors)
                mesh.update()

            elif event.key == '5':
                shading_filter.shading = self._previous_shading
                wireframe_filter.enabled = False
                normals.visible = False
                face_colors = []
                for i in range(len(self.shape.faces)):
                    face_colors.append(np.array([0.6, 0.6, 0.6]))
                mesh.set_data(self.shape.vertices, self.shape.faces, face_colors=face_colors)
                mesh.update()

            elif event.key == '6':
                shading_filter.shading = self._previous_shading
                wireframe_filter.enabled = False
                normals.visible = False
                face_colors = []
                for i in range(len(self.shape.faces)):
                    face_colors.append(np.array([0.6, 0.6, 0.6]))
                mesh.set_data(self.shape.vertices, self.shape.faces, face_colors=face_colors)
                mesh.update()

            elif event.key == 'z':
                self.f_func = scattering.f_lambert
                self.get_fluxes()
                plot_fluxes(self.phi_e)

            elif event.key == 'x':
                self.f_func = scattering.f_lommel
                self.get_fluxes()
                plot_fluxes(self.phi_e)

            elif event.key == 'c':
                scattering.B0 = 1.32
                scattering.minh = 0.20
                scattering.ming = -0.35
                scattering.bartheta = 10.0 * np.pi/180

                scattering.init_hapke(self.alpha)
                self.f_func = scattering.f_hapke
                self.get_fluxes()
                plot_fluxes(self.phi_e)

            elif event.key == 'm':
                if shading_filter.shading == 'smooth':
                    shading_filter.shading = 'flat'
                else: 
                    shading_filter.shading = 'smooth'

                mesh.update()

            elif event.key == 's':
                ox, oy, oz = self.s
                r = np.linalg.norm(self.s)
                azimuth = 180 - np.degrees(np.arctan2(ox, oy))
                elevation = np.degrees(np.arcsin(oz / r))
                self.view.camera = vispy.scene.cameras.TurntableCamera(
                    center=(0, 0, 0),
                    azimuth=azimuth,
                    elevation=elevation,
                    fov=0
                )

            elif event.key == 'o':
                ox, oy, oz = self.o
                r = np.linalg.norm(self.o)
                azimuth = 180 - np.degrees(np.arctan2(ox, oy))
                elevation = np.degrees(np.arcsin(oz / r))
                self.view.camera = vispy.scene.cameras.TurntableCamera(
                    center=(0, 0, 0),
                    azimuth=azimuth,
                    elevation=elevation,
                    fov=0
                )

            elif event.key == 'p':
                vispy.io.write_png("out/vispy_screenshot.png", vispy.gloo.util._screenshot())

            elif event.key == 'h':
                visible = not self.overlays[0].visible
                for v in self.overlays:
                    v.visible = visible

        plot_fluxes(self.phi_e)
        self.canvas.show()