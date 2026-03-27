#!/usr/bin/env python3

import argparse
import numpy as np
import os
import vispy.app

from .core import Asteroid
from . import io

def parse_vector(vec_str):
    res = list(map(float, vec_str.split(',')))
    if len(res) != 3:
        raise argparse.ArgumentTypeError("Vector must have exactly 3 components: x,y,z")
    return np.array(res, dtype=np.double)

def main():
    parser = argparse.ArgumentParser(
        description="TVET CLI - Asteroid Visualization & Analysis",
        usage="""
        tvet <obj_file> [options]
            Example: tvet asteroid.obj --plot-light-curve

        tvet --jpl <id> --start <t> --stop <t> --step <dt> [options]
            Example: tvet --jpl 103 --start 2026-02-16T00:00 --stop 2026-02-17T00:00 --step 10m

        tvet --damit <id> [options]
            Example: tvet --damit 10809

        tvet <obj_file> --jpl <id> --start <t> --stop <t> --step <dt> [options]
            Example: tvet asteroid.obj --jpl 103 --start 2026-02-16T00:00 --stop 2026-02-17T00:00 --step 10m
        
        tvet --damit <id> --jpl <id> --start <t> --stop <t> --step <dt> [options]
            Example: tvet --damit 10809 --jpl 103 --start 2026-02-16T00:00 --stop 2026-02-17T00:00 --step 10m
        """
    )

    parser.add_argument("filename", nargs="?", help="Path to local OBJ file (optional).")

    parser.add_argument("-j", "--jpl", dest="jpl_id", type=int, help="JPL Horizons body id (integer)")
    parser.add_argument("-d", "--damit", dest="damit_id", type=int, help="DAMIT asteroid/model id (integer)")

    parser.add_argument('--s', type=parse_vector, default=None, help="Incident light vector, e.g. '1,0,0'")
    parser.add_argument('--o', type=parse_vector, default=None, help="Observer vector, e.g. '0,0,1'")
    parser.add_argument("--scattering", choices=["lambert", "lommel", "hapke"], default="lambert", help="Scattering law to use: lambert, lommel, or hapke (default: lambert)")
    parser.add_argument("--geometry", action="store_true", help="Returns centers and normals of the asteroid triangle mesh and saves them to out/centers.txt and out/normals.txt")
    parser.add_argument("--cosines", action="store_true", help="Returns mu_i and mu_e and saves them to out/mu_i.txt and out/mu_e.txt")
    parser.add_argument("--fluxes", action="store_true", help="Returns phi_i, phi_e, and total flux and saves them to out/phi_i.txt, out/phi_e.txt, and out/total_flux.txt")
    parser.add_argument("--spin", type=float, nargs=5, metavar=("PERIOD", "EPOCH", "L", "B", "PHI0"), help="Spin state: period, epoch, l, b, phi0. Example: --spin 4.0 0.0 1.2 1.57 0.0")
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

    has_obj = isinstance(args.filename, str) and args.filename.lower().endswith(".obj")
    has_jpl = args.jpl_id is not None
    has_damit = args.damit_id is not None

    # EXCEPTIONS - coliding/incomplete arguments
    if args.filename is not None and not has_obj:
        parser.error("Positional argument must be a path to a .obj file.")
    if has_obj and not os.path.isfile(args.filename):
        parser.error(f"OBJ file not found: {args.filename}")

    if has_jpl and not args.start_time:
        parser.error("--jpl requires --start.")
    if (args.start_time or args.stop_time or args.step_size) and not has_jpl:
        parser.error("--start/--stop/--step are only valid together with --jpl.")
    
    if args.spin is not None and has_damit:
        parser.error("Can not use --spin and --damit together, they provide coliding information.")

    # Determine mode
    mode_jpl_only = has_jpl and (not has_damit) and (not has_obj)
    mode_damit_only = has_damit and (not has_jpl) and (not has_obj)
    mode_local_obj = has_obj  # may be combined with --jpl and/or --damit
    mode_jpl_damit = has_jpl and has_damit and (not has_obj)

    if not (mode_jpl_only or mode_damit_only or mode_local_obj or mode_jpl_damit):
        parser.error("Invalid invocation. Provide a local .obj file, or use --jpl, or use --damit, or use --jpl + --damit.")

    out_dir = "out"
    if not args.no_save:
        os.makedirs(out_dir, exist_ok=True)

    asteroid = Asteroid(args=args, filename=args.filename)

    # MARK: - jpl only

    if mode_jpl_only:
        if args.stop_time:
            s_unit, o_unit = asteroid.get_ephems(
                body=str(args.jpl_id),
                start_time=args.start_time,
                stop_time=args.stop_time,
                step_size=args.step_size,
                observer_center=args.observer_center,
                sun_center=args.sun_center,
                normalize=True,
                timeout=args.timeout
            )

            data = np.hstack([s_unit, o_unit])
        else:
            s_unit, o_unit = asteroid.get_single_ephem(
                body=str(args.jpl_id),
                epoch=args.start_time,
                observer_center=args.observer_center,
                sun_center=args.sun_center,
                normalize=True,
                timeout=args.timeout
            )

            data = np.hstack([s_unit, o_unit]).reshape(1, 6)

            if args.verbose > 0:
                args.verbose = 1

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "ephemerides.txt"), data, header="sx sy sz ox oy oz")

        if not args.quiet:
            print(f"\nEphemerides: body [{args.jpl_id}] points [{data.shape[0]}]\n")

            if not args.no_save:
                print(f"Saved vectors to {out_dir}/ephemerides.txt\n")

            if args.verbose > 0:
                print(
                    "Vectors [sx sy sz ox oy oz] " +
                    "(first " + (f"{args.verbose}):" if args.verbose <= data.shape[0] else f"{data.shape[0]}):") +
                    f"\n {data[:args.verbose]} \n" +
                    ("...\n" if args.verbose < data.shape[0] else "")
                )

        return
    
    # MARK: - damit only

    if mode_damit_only:
        asteroid.get_damit(model_id=args.damit_id)

        if not args.no_save:
            io.save_obj_file(
                path=os.path.join(out_dir, f"shape_{args.damit_id}.obj"), 
                vertices=asteroid.shape.vertices,
                faces=asteroid.shape.faces
            )

            io.save_spin(
                path=os.path.join(out_dir, f"spin_{args.damit_id}.txt"),
                l=asteroid.light_curve.l, 
                b=asteroid.light_curve.b, 
                period=asteroid.light_curve.period, 
                epoch=asteroid.light_curve.epoch, 
                phi0=asteroid.light_curve.phi0
            )

        if not args.quiet:
            len_v = len(asteroid.shape.vertices)
            len_f = len(asteroid.shape.faces)

            print(f"\nShape model: body {args.damit_id}, " + 
                f"vertices [{len_v}], " +
                f"faces [{len_f}]"
            )
            print(f"Spin:\nl {asteroid.light_curve.l} " + 
                f"b {asteroid.light_curve.b} " +
                f"period {asteroid.light_curve.period}\n" +
                f"epoch {asteroid.light_curve.epoch} " +
                f"phi0 {asteroid.light_curve.phi0}\n"
            )

            if not args.no_save:
                print(f"Saved shape model to {out_dir}/shape_{args.damit_id}.obj")
                print(f"Saved spin to {out_dir}/spin_{args.damit_id}.txt\n")

            if args.verbose > 0:
                print(
                    "Vertices " +
                    "(first " + (f"{args.verbose}" if args.verbose <= len_v else f"{len_v}") +
                    f" of {len_v}):\n {asteroid.shape.vertices[:args.verbose]} \n"
                )

                print(
                    "Faces " +
                    "(first " + (f"{args.verbose}" if args.verbose <= len_f else f"{len_f}") +
                    f" of {len_f}):\n {asteroid.shape.faces[:args.verbose]} \n"
                )

            if args.interactive_plot:
                asteroid.interactive_plot()
                vispy.app.run()

        return

    if args.s is not None:
        asteroid.s = args.s
    if args.o is not None:
        asteroid.o = args.o

    if args.spin is not None:
        period, epoch, l, b, phi0 = args.spin

        if period <= 0.0:
            parser.error(f"--spin PERIOD must be > 0, got: {period}")

        asteroid.light_curve.period = float(period)
        asteroid.light_curve.epoch = float(epoch)
        asteroid.light_curve.l = float(l)
        asteroid.light_curve.b = float(b)
        asteroid.light_curve.phi0 = float(phi0)

    if has_jpl:
        s_unit, o_unit = asteroid.get_ephems(
            body=str(args.jpl_id),
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

    if has_damit:
        asteroid.get_damit(model_id=args.damit_id)

    # MARK: - options

    out_flags = [
        args.geometry,
        args.cosines,
        args.fluxes,
        args.light_curve,
        args.plot_light_curve
    ]

    if args.geometry:
        asteroid.get_geometry()

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "centers.txt"), asteroid.shape.centers)
            np.savetxt(os.path.join(out_dir, "normals.txt"), asteroid.shape.normals)

        if not args.quiet:
            len_c = len(asteroid.shape.centers)
            len_n = len(asteroid.shape.normals)

            print(f"\nGeometry: centers [{len_c}], normals [{len_n}]\n")

            if not args.no_save:
                print(f"Saved centers to {out_dir}/centers.txt")
                print(f"Saved normals to {out_dir}/normals.txt\n")

            if args.verbose > 0:
                print(f"Centers " + 
                    "(first " + (f"{args.verbose}" if args.verbose <= len_c else f"{len_c}") +
                    f" of {len_c}):\n {asteroid.shape.centers[:args.verbose]}\n"
                )
                print(f"Normals " + 
                    "(first " + (f"{args.verbose}" if args.verbose <= len_n else f"{len_n}") +
                    f" of {len_n}):\n {asteroid.shape.normals[:args.verbose]}\n"
                )

    if args.cosines:
        asteroid.get_cosines()

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "mu_i.txt"), asteroid.mu_i)
            np.savetxt(os.path.join(out_dir, "mu_e.txt"), asteroid.mu_e)

        if not args.quiet:
            len_i = len(asteroid.mu_i)
            len_e = len(asteroid.mu_e)

            print(f"\nCosines: mu_i [{len_i}], mu_e [{len_e}]\n")

            if not args.no_save:
                print(f"Saved mu_i to {out_dir}/mu_i.txt")
                print(f"Saved mu_e to {out_dir}/mu_e.txt\n")

            if args.verbose > 0:
                print(f"mu_i " + 
                    "(first " + (f"{args.verbose}" if args.verbose <= len_i else f"{len_i}") +
                    f" of {len_i}):\n {asteroid.mu_i[:args.verbose]}\n"
                )

                print(f"mu_e " + 
                    "(first " + (f"{args.verbose}" if args.verbose <= len_e else f"{len_e}") +
                    f" of {len_e}):\n {asteroid.mu_e[:args.verbose]}\n"
                )
        
    if args.fluxes:
        asteroid.get_fluxes()

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "phi_i.txt"), asteroid.phi_i)
            np.savetxt(os.path.join(out_dir, "phi_e.txt"), asteroid.phi_e)
            np.savetxt(os.path.join(out_dir, "total_flux.txt"), np.array([asteroid.total]))

        if not args.quiet:
            len_i = len(asteroid.phi_i)
            len_e = len(asteroid.phi_e)

            print(f"\nFluxes: phi_i [{len_i}], phi_e [{len_e}]")
            print(f"Total flux: {asteroid.total}\n")

            if not args.no_save:
                print(f"Saved phi_i to {out_dir}/phi_i.txt")
                print(f"Saved phi_e to {out_dir}/phi_e.txt")
                print(f"Saved total flux to {out_dir}/total_flux.txt\n")

            if args.verbose > 0:
                print(f"phi_i " + 
                    "(first " + (f"{args.verbose}" if args.verbose <= len_i else f"{len_i}") +
                    f" of {len_i}):\n {asteroid.phi_i[:args.verbose]}\n"
                )

                print(f"phi_e " + 
                    "(first " + (f"{args.verbose}" if args.verbose <= len_e else f"{len_e}") +
                    f" of {len_e}):\n {asteroid.phi_e[:args.verbose]}\n"
                )

    if args.light_curve:
        curve_points = asteroid.get_light_curve_for_period()

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "light_curve.txt"), curve_points)

        if not args.quiet:
            len_c = len(curve_points)

            print(f"\nLight curve: points [{len_c}]\n")

            if not args.no_save:
                print(f"Saved light curve to {out_dir}/light_curve.txt\n")

            if args.verbose > 0:
                print(f"Curve points " + 
                    "(first " + (f"{args.verbose}" if args.verbose <= len_c else f"{len_c}") +
                    f" of {len_c}):\n {curve_points[:args.verbose]}\n"
                )

    if args.plot_light_curve:
        curve_points = asteroid.get_light_curve_for_period()

        if not args.quiet:
            asteroid.plot_light_curve(curve_points)

    # Interactive plot is default if no out flags are set, or if explicitly requested
    if not args.quiet and (args.interactive_plot or not any(out_flags)):
        asteroid.interactive_plot()
        vispy.app.run()

if __name__ == "__main__":
    main()
