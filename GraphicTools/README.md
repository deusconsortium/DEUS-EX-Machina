# GraphicTools

## Abstract

blabla

## File format

header of file:

- Nprofiles: (long int) total number of profiles
- Nradius: (int) number of radius for each profile
- R0: (float) value of the first schell radius in Ramses Unit
- Dr: (float) step size for the radius in Ramses Unit

for each profile:

- f[]: (float) array with values of m(r)/m_back(r) ; where m_back(r) = 4Pi/3*r^3*rho_m
- v[]: (float) array with peculiar velocities (Unit ???)

## Ramses Unit

For a simulation with BOXLEN and NPART, the value of a given size in [Mpc/h] comobile is simply computed by

X(Mpc/h) = X(Ramses)*BOXLEN/NPART
