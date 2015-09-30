# DEUS-EX-Machina

## Dark Energy Universe Simulations Extrema eXploration machina

This toolkit is aimed to explore/analyse extrema statistics ie minima/voids maxima/halos in the DEUS cosmology simulations.

More infos on http://www.deus-consortium.org/

It is composed of the following softwares:

- compute_extrema: MPI Fortran code to detect extrema from particle positions. Calculate density on a grid using CIC, smooth result, then find minima or maxima and write in a FOFextrema file.

- VoidTracer: C++ OpenMP code to trace void/halo from FOFextrema file, and get values: R, delta, avg_delta.

- FOFReaderLib: C++ library to read the different file DEUS formats, include FOFextrema.

- Splotch: C++ raytracer to make images/movies from particle files, modified to read DEUS file formats

- GraphicTools: Python code to display density/mass/velocity profiles




