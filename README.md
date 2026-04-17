Python package for visualization and analysis of asteroid shape models and their light curves.

<div align="center">
  <a href="https://github.com/brozrost/tvet/actions">
    <img src="https://github.com/brozrost/tvet/actions/workflows/python-package.yml/badge.svg">
  </a>
  <a href="https://github.com/brozrost/tvet/graphs/contributors">
    <img src="https://img.shields.io/github/contributors/brozrost/tvet">
  </a>
  <a href="https://github.com/brozrost/tvet/issues">
    <img src="https://img.shields.io/github/issues/brozrost/tvet">
  </a>
  <a href="https://github.com/brozrost/tvet/pulls">
    <img src="https://img.shields.io/github/issues-pr/brozrost/tvet">
  </a>
</div>

## Installation

`tvet` is available on PyPI, download with:

```bash
pip install tvet
```

## CLI usage

After installation, you can use the `tvet` command-line tool:

```bash
tvet --help
```

```bash
tvet asteroid.obj
```

```bash
tvet --jpl 103 --start 2026-02-16T00:00 --stop 2026-02-17T00:00 --step 10m
```

```bash
tvet --damit 10809
```

```bash
tvet asteroid.obj --jpl 103 --start 2026-02-16T00:00 --stop 2026-02-17T00:00 --step 10m
```

```bash
tvet --damit 10809 --jpl 103 --start 2026-02-16T00:00 --stop 2026-02-17T00:00 --step 10m
```

### Options
```bash
  -h, --help            show this help message and exit
  -j, --jpl JPL_ID      JPL Horizons body id (integer)
  -d, --damit DAMIT_ID  DAMIT asteroid/model id (integer)
  --s S                 Incident light vector, e.g. '1,0,0'
  --o O                 Observer vector, e.g. '0,0,1'
  --scattering {lambert,lommel,hapke}
                        Scattering law to use: lambert, lommel, or hapke (default: lambert)
  --geometry            Returns centers and normals of the asteroid triangle mesh and saves them to out/centers.txt and out/normals.txt
  --cosines             Returns mu_i and mu_e and saves them to out/mu_i.txt and out/mu_e.txt
  --fluxes              Returns phi_i, phi_e, and total flux and saves them to out/phi_i.txt, out/phi_e.txt, and out/total_flux.txt
  --spin L B PERIOD EPOCH PHI0
                        Spin state: l, b, period, epoch, phi0. Example: --spin 4.0 0.0 1.2 1.57 0.0
  -l, --light-curve     Save the asteroid light curve points to out/light_curve.txt and plot the light curve
  -i, --interactive-plot
                        Plot the interactive asteroid geometry and light curve
  --shininess SHININESS
                        Shininess factor for the asteroid surface
  --wireframe-width WIREFRAME_WIDTH
                        Width of the wireframe lines
  --start START_TIME    Ephemeris start time (e.g. 2026-02-16T00:00)
  --stop STOP_TIME      Ephemeris stop time (e.g. 2026-02-17T00:00)
  --step STEP_SIZE      Ephemeris step size (e.g. 10m, 1h, 1d)
  --observer-center OBSERVER_CENTER
                        Horizons center for observer (default: 500@399)
  --sun-center SUN_CENTER
                        Horizons center for Sun (default: 500@10)
  --timeout TIMEOUT     Horizons request timeout in seconds (default: 30)
  -v, --verbose [VERBOSE]
                        Show preview of first n elements (defaults to 3)
  -q, --quiet           Suppress all stdout output (files are still saved)
  --no-save             Do not create output files/directories
  --no-normalize        Disable normalization of all vectors before computation
```

For more details, see the [Key concepts](https://github.com/scraptechguy/tvet/blob/main/docs/CONCEPTS.md).

## Resources

+ <b><a href="https://github.com/scraptechguy/tvet/blob/main/docs/CONCEPTS.md">Key concepts</a></b> (explains all key concepts used in tvet)
+ <b><a href="https://sirrah.troja.mff.cuni.cz/~mira/tmp/diplomky/Broz_2024.pdf">Dokumentace v češtině</a></b>
