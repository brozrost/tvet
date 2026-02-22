#!/usr/bin/env python3

import argparse
import numpy as np
import os
import vispy.app
from .core import Asteroid

def parse_vector(vec_str):
    res = list(map(float, vec_str.split(',')))
    if len(res) != 3:
        raise argparse.ArgumentTypeError("Vector must have exactly 3 components: x,y,z")
    return np.array(res, dtype=np.double)

def main():
    parser = argparse.ArgumentParser(
        description="TVET CLI - Asteroid Visualization & Analysis",
        usage="""
        tvet <filename|body_id> [body_id] [options]

        Example (OBJ, static vectors):
            tvet asteroid.obj --plot-light-curve --s 1,0,0 --o 0,0,1

        Example (OBJ + Horizons ephemerides):
            tvet asteroid.obj 499 --start 2026-02-16T00:00 --stop 2026-02-17T00:00 --step 10m --fluxes

        Example (ephemerides only):
            tvet 499 --start 2026-02-16T00:00 --stop 2026-02-17T00:00 --step 10m
        """
    )

    parser.add_argument("filename", nargs="?", help="Path to OBJ file OR numeric body id (e.g. 499)")
    parser.add_argument("body_id", nargs="?", help="Numeric body id for Horizons (e.g. 499). If provided with an OBJ, s/o are taken from Horizons.")

    parser.add_argument('--s', type=parse_vector, default=None, help="Incident light vector, e.g. '1,0,0'")
    parser.add_argument('--o', type=parse_vector, default=None, help="Observer vector, e.g. '0,0,1'")
    parser.add_argument("--scattering", choices=["lambert", "lommel", "hapke"], default="lambert", help="Scattering law to use: lambert, lommel, or hapke (default: lambert)")
    parser.add_argument("--geometry", action="store_true", help="Returns centers and normals of the asteroid triangle mesh and saves them to out/centers.txt and out/normals.txt")
    parser.add_argument("--cosines", action="store_true", help="Returns mu_i and mu_e and saves them to out/mu_i.txt and out/mu_e.txt")
    parser.add_argument("--fluxes", action="store_true", help="Returns phi_i, phi_e, and total flux and saves them to out/phi_i.txt, out/phi_e.txt, and out/total_flux.txt")
    parser.add_argument("--light-curve", action="store_true", help="Save the asteroid light curve points to out/light_curve.txt")
    parser.add_argument("--plot-light-curve", action="store_true", help="Plot the asteroid light curve using matplotlib")

    parser.add_argument("--interactive-plot", action="store_true", help="Plot the interactive asteroid geometry and light curve")
    parser.add_argument('--shininess', type=float, default=100, help="Shininess factor for the asteroid surface")
    parser.add_argument('--wireframe-width', type=float, default=1, help="Width of the wireframe lines")

    parser.add_argument("--start", dest="start_time", help="Ephemeris start time (e.g. 2026-02-16T00:00)")
    parser.add_argument("--stop", dest="stop_time", help="Ephemeris stop time (e.g. 2026-02-17T00:00)")
    parser.add_argument("--step", dest="step_size", help="Ephemeris step size (e.g. 10m, 1h, 1d)")

    parser.add_argument("--observer-center", default="500@399", help="Horizons center for observer (default: 500@399)")
    parser.add_argument("--sun-center", default="500@10", help="Horizons center for Sun (default: 500@10)")
    parser.add_argument("--timeout", type=float, default=30.0, help="Horizons request timeout in seconds (default: 30)")

    parser.add_argument("-v", "--verbose", type=int, nargs="?", const=3, default=0, help="Show preview of first n elements (defaults to 3)")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress all stdout output (files are still saved)")
    parser.add_argument("--no-save", action="store_true", help="Do not create output files/directories")

    args = parser.parse_args()

    if args.filename is None:
        parser.error("Missing target. Use `tvet asteroid.obj ...` or `tvet 499 --start ... --stop ... --step ...`")

    # Determine mode
    filename_looks_like_obj = isinstance(args.filename, str) and args.filename.lower().endswith(".obj")
    filename_is_int = False
    filename_int = None

    if not filename_looks_like_obj:
        try:
            filename_int = int(str(args.filename))
            filename_is_int = True
        except ValueError:
            filename_is_int = False

    eph_only = filename_is_int and (args.body_id is None)
    obj_plus_eph = filename_looks_like_obj and (args.body_id is not None)
    obj_only = filename_looks_like_obj and (args.body_id is None)

    # Check different mode exceptions
    if obj_plus_eph and ((args.s is not None) or (args.o is not None)):
        parser.error("Do not use --s/--o together with body_id. Ephemerides define s/o.")

    if args.s is None: 
        args.s = np.array([1.0, 0.0, 0.0], dtype=np.double)
    if args.o is None: 
        args.o = np.array([0.0, 0.0, 1.0], dtype=np.double)

    if not (eph_only or obj_plus_eph or obj_only):
        parser.error("Invalid invocation. Use `tvet asteroid.obj ...` or `tvet asteroid.obj 499 --start --stop --step ...` or `tvet 499 --start --stop --step`.")

    if eph_only or obj_plus_eph:
        if not (args.start_time and args.stop_time and args.step_size):
            parser.error("Ephemeris mode requires --start, --stop, and --step.")

    out_dir = "out"
    if not args.no_save:
        os.makedirs(out_dir, exist_ok=True)

    # ===================================
    # MODE 1: eph-only (save eph to file)
    # ===================================
    if eph_only:
        body_id_int = filename_int

        asteroid = Asteroid(args=args, filename=None)
        s_unit, o_unit = asteroid.get_ephems(
            body=str(body_id_int),
            start_time=args.start_time,
            stop_time=args.stop_time,
            step_size=args.step_size,
            observer_center=args.observer_center,
            sun_center=args.sun_center,
            normalize=True,
            timeout=args.timeout
        )

        data = np.hstack([s_unit, o_unit])

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "ephemerides.txt"), data, header="sx sy sz ox oy oz")

        if not args.quiet:
            print(f"\nEphemerides: body [{body_id_int}] points [{data.shape[0]}]\n")

            if not args.no_save:
                print(f"Saved vectors to {out_dir}/ephemerides.txt\n")

            if args.verbose > 0:
                print(f"Vectors [sx sy sz ox oy oz] (first {args.verbose}):\n {data[:args.verbose]}\n ...\n")

        return

    if filename_looks_like_obj and not os.path.isfile(args.filename):
        parser.error(f"OBJ file not found: {args.filename}")

    asteroid = Asteroid(args=args, filename=args.filename)
    asteroid.s = args.s
    asteroid.o = args.o

    out_flags = [
        args.geometry,
        args.cosines,
        args.fluxes,
        args.light_curve,
        args.plot_light_curve
    ]

    # ==================
    # MODE 2: obj + eph
    # ==================
    if obj_plus_eph:
        try:
            body_id_int = int(str(args.body_id))
        except ValueError:
            parser.error(f"body_id must be an integer (e.g. 499), got: {args.body_id}")

        s_unit, o_unit = asteroid.get_ephems(
            body=str(body_id_int),
            start_time=args.start_time,
            stop_time=args.stop_time,
            step_size=args.step_size,
            observer_center=args.observer_center,
            sun_center=args.sun_center,
            normalize=True,
            timeout=args.timeout
        )

        asteroid.s_array = s_unit
        asteroid.o_array = o_unit

        asteroid.s = s_unit[0]
        asteroid.o = o_unit[0]

    # =============================================
    # MODE 3: obj + eph || obj with static vectors
    # =============================================
    if args.geometry:
        asteroid.get_geometry()

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "centers.txt"), asteroid.centers)
            np.savetxt(os.path.join(out_dir, "normals.txt"), asteroid.normals)

        if not args.quiet:
            print(f"\nGeometry: centers [{len(asteroid.centers)}], normals [{len(asteroid.normals)}]\n")

            if not args.no_save:
                print(f"Saved centers to {out_dir}/centers.txt")
                print(f"Saved normals to {out_dir}/normals.txt\n")

            if args.verbose > 0:
                print(f"Centers: (first {args.verbose})\n {asteroid.centers[:args.verbose]}\n")
                print(f"Normals: (first {args.verbose})\n {asteroid.normals[:args.verbose]}\n")

    if args.cosines:
        asteroid.get_cosines()

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "mu_i.txt"), asteroid.mu_i)
            np.savetxt(os.path.join(out_dir, "mu_e.txt"), asteroid.mu_e)

        if not args.quiet:
            print(f"\nCosines: mu_i [{len(asteroid.mu_i)}], mu_e [{len(asteroid.mu_e)}]\n")

            if not args.no_save:
                print(f"Saved mu_i to {out_dir}/mu_i.txt")
                print(f"Saved mu_e to {out_dir}/mu_e.txt\n")

            if args.verbose > 0:
                print(f"mu_i: (first {args.verbose})\n {asteroid.mu_i[:args.verbose]}\n")
                print(f"mu_e: (first {args.verbose})\n {asteroid.mu_e[:args.verbose]}\n")
        
    if args.fluxes:
        asteroid.get_fluxes()

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "phi_i.txt"), asteroid.phi_i)
            np.savetxt(os.path.join(out_dir, "phi_e.txt"), asteroid.phi_e)
            np.savetxt(os.path.join(out_dir, "total_flux.txt"), np.array([asteroid.total]))

        if not args.quiet:
            print(f"\nFluxes: phi_i [{len(asteroid.phi_i)}], phi_e [{len(asteroid.phi_e)}]")
            print(f"Total flux: {asteroid.total}\n")

            if not args.no_save:
                print(f"Saved phi_i to {out_dir}/phi_i.txt")
                print(f"Saved phi_e to {out_dir}/phi_e.txt")
                print(f"Saved total flux to {out_dir}/total_flux.txt\n")

            if args.verbose > 0:
                print(f"phi_i: (first {args.verbose})\n {asteroid.phi_i[:args.verbose]}\n")
                print(f"phi_e: (first {args.verbose})\n {asteroid.phi_e[:args.verbose]}\n")
        

    if args.light_curve:
        curve_points = asteroid.get_light_curve()

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "light_curve.txt"), curve_points)

        if not args.quiet:
            print(f"\nLight curve: points [{len(curve_points)}]\n")

            if not args.no_save:
                print(f"Saved light curve to {out_dir}/light_curve.txt\n")

            if args.verbose > 0:
                print(f"Curve points: (first {args.verbose})\n {curve_points[:args.verbose]}\n")

    if args.plot_light_curve:
        curve_points = asteroid.get_light_curve()

        if not args.quiet:
            asteroid.plot_light_curve(curve_points)

    # Interactive plot is default if no out flags are set, or if explicitly requested
    if not args.quiet and (args.interactive_plot or not any(out_flags)):
        asteroid.interactive_plot()
        vispy.app.run()

if __name__ == "__main__":
    main()
