import numpy as np

class LightCurve:
    def __init__(
        self, 
        asteroid, 
        *, 
        f_func=None, 
        start=None, 
        stop=None, 
        step=None
    ):
        if asteroid is None:
            raise ValueError("LightCurve requires an Asteroid instance.")
        self.asteroid = asteroid

        self.f_func = f_func
        self.start = 0 if start is None else start
        self.stop = None if stop is None else stop
        self.step = None if step is None else step

    def compute(self):
        pass

    def compute_for_period(self, s=None, o=None, n=100, start=None, l=None, b=None, period=None, epoch=None, phi0=None):
        asteroid = self.asteroid
        
        if s is None: 
            s = asteroid.s
        if o is None: 
            o = asteroid.o
        if start is None:
            start = self.start
        if l is None:
            l = asteroid.l
        if b is None:
            b = asteroid.b
        if period is None:
            period = asteroid.period
        if epoch is None:
            epoch = asteroid.epoch
        if phi0 is None:
            phi0 = asteroid.phi0

        # Checks
        if n is None or int(n) <= 0:
            raise ValueError("n must be a positive integer.")
        n = int(n)
        
        if period is None or period == 0:
            raise ValueError("period must be non-zero.")

        phi2 = np.pi / 2.0 - b
        phi3 = l
        
        total = np.empty((n, 2), dtype=np.double)
        f_func = self.f_func if self.f_func is not None else asteroid.f_func

        for i in range(n):
            t = start + period * i/n
            phi1 = 2.0 * np.pi * (t - epoch) / period + phi0

            s_ = asteroid._rotate_z(s, phi1)
            s_ = asteroid._rotate_y(s_, phi2)
            s_ = asteroid._rotate_z(s_, phi3)

            o_ = asteroid._rotate_z(o, phi1)
            o_ = asteroid._rotate_y(o_, phi2)
            o_ = asteroid._rotate_z(o_, phi3)

            asteroid.get_fluxes(s=s_, o=o_, f_func=f_func)
            total[i, 0] = t
            total[i, 1] = asteroid.total

        # Normalize
        denx = np.max(total[:, 0]) - np.min(total[:, 0])
        total[:, 0] = (total[:, 0] - np.min(total[:, 0])) / denx if denx > 0 else 0.0

        deny = np.max(total[:, 1]) - np.min(total[:, 1])
        total[:, 1] = -(total[:, 1] - np.min(total[:, 1])) / deny if deny > 0 else 0.0

        return total