# GraphicTools

## Abstract

blabla

## DEUSgraphics.py:

Includes a python class DEUSgraphics that can be used to load and plot several quantitites from the .DEUSprofile files. The class can be used in the following way

```
#define an object that will study over-densities from the file 'my_data_file.DEUSprofile'
plotter = DEUSgraphics(True)
plotter.load('my_data_file')

#this last call be be also done in the constructor
plotter = DEUSgraphics(True,'my_data_file')

#plot a mean profile with a comoving r1 in [10.0,10.0 + 1.0] [Mpc/h]
plotter.PlotMeanProfile(10.0,1.0)

#plot the P(r1) histogram
plotter.PlotStatistics()

#plot an individual profile with index 'i' (must verify i >= 0 and i < plotter.getProfileNumber())
plotter.PlotSingleProfile(i) 
```

## File format

header of file:

begining with the simulation parameters

- H0: (float) value of the H0 in [kms/s/Mpc]
- a: (float) experimental value of the scale factor (a = 1 today)
- Wm: (float) value of Omega_matter at z = 0 (i.e. a = 1)
- n: (int) size of the cosmo char[]
- cosmo: (char[]) the cosmology used in the simulation (e.g. 'lcdmw5')
- boxlen: (int) the size of the box of the simulation (e.g. 5184)
- npart: (int) the number of particles in each direction (i.e. the total number of particles in the simulation is npart^3)

properties of this run

- isOverDensity: (bool) 1 if profiles are over-densities, 0 elsewhere
- Nprofiles: (long int) total number of profiles
- Nradius: (int) number of radius for each profile
- R0: (float) value of the first schell radius in Ramses Unit
- Dr: (float) step size for the radius in Ramses Unit

for each profile:

- f[]: (float) array with values of m(r)/m_back(r) ; where m_back(r) = 4Pi/3*r^3*rho_m
- v[]: (float) array with peculiar velocities (Unit ???)

## Ramses Unit

For a simulation with boxlen = BOXLEN, the value of a given size in [Mpc/h] comobile is simply computed by

X(Mpc/h) = X(Ramses)*BOXLEN
