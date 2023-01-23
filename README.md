# Molecular_dynamics
MD_2D.py code computates Temperature, Kinetic, Potential and Total energies of 25 particles in a box of length=5. 
These particles experience pair-wise interaction through Lennard-Jones Potential.
mass of particles = 1. 
Potential cutoff = 3*sigma.
epsilon = 1.
PBC => Periodic Boundary conditions.
This code is devoid of Neighbour list algorithm, so it works faster only for bunch of particles like, 25, 50 etc.
This code supports NVE ensemble by default. No other ensemble is implemented yet.
