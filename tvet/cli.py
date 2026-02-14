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
    parser.add_argument("--geometry", action="store_true", help="Retruns centers and normals of the asteroid triangle mesh and saves them to output/centers.txt and output/normals.txt")
    parser.add_argument("--cosines", action="store_true", help="Returns mu_i and mu_e and saves them to output/mu_i.txt and output/mu_e.txt")
    parser.add_argument("--fluxes", action="store_true", help="Returns phi_i, phi_e, and total flux and saves them to output/phi_i.txt, output/phi_e.txt, and output/total_flux.txt")
    parser.add_argument("--light-curve", action="store_true", help="Save the asteroid light curve points to output/light_curve.txt")
    parser.add_argument("--plot-light-curve", action="store_true", help="Plot the asteroid light curve using matplotlib")
    parser.add_argument("--interactive-plot", action="store_true", help="Plot the interactive asteroid geometry and light curve")
    parser.add_argument('--shininess', default=100, help="Shininess factor for the asteroid surface")
    parser.add_argument('--wireframe-width', default=1, help="Width of the wireframe lines")

    parser.add_argument("-v", "--verbose", action="store_true", help="Print summaries and small previews to stdout")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress all stdout output (files are still saved)")
    
    args = parser.parse_args()

    asteroid = Asteroid(args=args, filename=args.filename)
    asteroid.s = args.s
    asteroid.o = args.o

    output_dir = "out"
    os.makedirs(output_dir, exist_ok=True)

    output_flags = [
        args.geometry,
        args.cosines,
        args.fluxes,
        args.light_curve,
        args.plot_light_curve
    ]

    if args.geometry:
        asteroid.get_geometry()

        np.savetxt(os.path.join(output_dir, "centers.txt"), asteroid.centers)
        np.savetxt(os.path.join(output_dir, "normals.txt"), asteroid.normals)

        if not args.quiet:
            print(f"\nSaved centers to {output_dir}/centers.txt")
            print(f"Saved normals to {output_dir}/normals.txt\n")

            if args.verbose:
                print("Centers:", len(asteroid.centers))
                print("Normals:", len(asteroid.normals))

                print("\nCenters:\n", asteroid.centers[:3])
                print("\nNormals:\n", asteroid.normals[:3], "\n")

    if args.cosines:
        asteroid.get_cosines(s=args.s, o=args.o)

        np.savetxt(os.path.join(output_dir, "mu_i.txt"), asteroid.mu_i)
        np.savetxt(os.path.join(output_dir, "mu_e.txt"), asteroid.mu_e)

        if not args.quiet:
            print(f"\nSaved mu_i to {output_dir}/mu_i.txt")
            print(f"Saved mu_e to {output_dir}/mu_e.txt\n")

            if args.verbose:
                print("Length of mu_i:", len(asteroid.mu_i))
                print("Length of mu_e:", len(asteroid.mu_e))

                print("\nmu_i:\n", asteroid.mu_i[:3])
                print("mu_e:\n", asteroid.mu_e[:3], "\n")
        
    if args.fluxes:
        asteroid.get_fluxes()

        np.savetxt(os.path.join(output_dir, "phi_i.txt"), asteroid.phi_i)
        np.savetxt(os.path.join(output_dir, "phi_e.txt"), asteroid.phi_e)
        np.savetxt(os.path.join(output_dir, "total_flux.txt"), np.array([asteroid.total]))

        if not args.quiet:
            print(f"\nSaved phi_i to {output_dir}/phi_i.txt")
            print(f"Saved phi_e to {output_dir}/phi_e.txt")
            print(f"Saved total flux to {output_dir}/total_flux.txt\n")

            if args.verbose:
                print("Length of phi_i:", len(asteroid.phi_i))
                print("Length of phi_e:", len(asteroid.phi_e))

                print("\nphi_i:\n", asteroid.phi_i[:3])
                print("phi_e:\n", asteroid.phi_e[:3])
                print("\nTotal flux:", asteroid.total, "\n")
        

    if args.light_curve:
        curve_points = asteroid.light_curve()

        np.savetxt(os.path.join(output_dir, "light_curve.txt"), curve_points)

        if not args.quiet:
            print(f"\nSaved light curve to {output_dir}/light_curve.txt\n")

    if args.plot_light_curve:
        curve_points = asteroid.light_curve()
        asteroid.plot_light_curve(curve_points)

    # Interactive plot is default if no output flags are set, or if explicitly requested
    if args.interactive_plot or not any(output_flags):
        asteroid.interactive_plot()
        vispy.app.run()

if __name__ == "__main__":
    main()
