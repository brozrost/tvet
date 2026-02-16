#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import os
import vispy.app
from .core import Asteroid

def parse_vector(vec_str):
    return np.array(list(map(float, vec_str.split(','))))

def main():
    parser = argparse.ArgumentParser(
        description="TVET CLI - Asteroid Visualization & Analysis",
        usage="""
        tvet <filename> [options]

        Example:
            tvet asteroid.obj --plot-light-curve --s 1,0,0 --o 0,0,1
        """
    )

    parser.add_argument("filename", nargs="?", help="Path to OBJ file")
    parser.add_argument('--s', type=parse_vector, default=np.array([1,0,0]), help="Incident light vector, e.g. '1,0,0'")
    parser.add_argument('--o', type=parse_vector, default=np.array([0,0,1]), help="Observer vector, e.g. '0,0,1'")
    parser.add_argument("--scattering", choices=["lambert", "lommel", "hapke"], default="lambert", help="Scattering law to use: lambert, lommel, or hapke (default: lambert)")
    parser.add_argument("--geometry", action="store_true", help="Retruns centers and normals of the asteroid triangle mesh and saves them to out/centers.txt and out/normals.txt")
    parser.add_argument("--cosines", action="store_true", help="Returns mu_i and mu_e and saves them to out/mu_i.txt and out/mu_e.txt")
    parser.add_argument("--fluxes", action="store_true", help="Returns phi_i, phi_e, and total flux and saves them to out/phi_i.txt, out/phi_e.txt, and out/total_flux.txt")
    parser.add_argument("--light-curve", action="store_true", help="Save the asteroid light curve points to out/light_curve.txt")
    parser.add_argument("--plot-light-curve", action="store_true", help="Plot the asteroid light curve using matplotlib")
    parser.add_argument("--interactive-plot", action="store_true", help="Plot the interactive asteroid geometry and light curve")
    parser.add_argument('--shininess', default=100, help="Shininess factor for the asteroid surface")
    parser.add_argument('--wireframe-width', default=1, help="Width of the wireframe lines")

    parser.add_argument(
        "-v", "--verbose", 
        type=int,
        nargs="?",
        const=3,
        default=0,
        help="Show preview of first n elements (defaults to 3)"
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress all stdout output (files are still saved)")
    parser.add_argument("--no-save", action="store_true", help="Do not create output files/directories")

    args = parser.parse_args()

    asteroid = Asteroid(args=args, filename=args.filename)
    asteroid.s = args.s
    asteroid.o = args.o

    out_dir = "out"

    if not args.no_save:
        os.makedirs(out_dir, exist_ok=True)

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
        asteroid.get_cosines(s=args.s, o=args.o)

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
        curve_points = asteroid.light_curve()

        if not args.no_save:
            np.savetxt(os.path.join(out_dir, "light_curve.txt"), curve_points)

        if not args.quiet:
            print(f"\nLight curve: points [{len(curve_points)}]\n")

            if not args.no_save:
                print(f"Saved light curve to {out_dir}/light_curve.txt\n")

            if args.verbose > 0:
                print(f"Curve points: (first {args.verbose})\n {curve_points[:args.verbose]}\n")

    if args.plot_light_curve:
        curve_points = asteroid.light_curve()

        if not args.quiet:
            asteroid.plot_light_curve(curve_points)

    # Interactive plot is default if no out flags are set, or if explicitly requested
    if not args.quiet and (args.interactive_plot or not any(out_flags)):
        asteroid.interactive_plot()
        vispy.app.run()

if __name__ == "__main__":
    main()
