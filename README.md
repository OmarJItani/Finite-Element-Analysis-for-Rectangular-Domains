# Finite-Element-Analysis-for-Rectangular-Domains

This code models a general rectangular topology doamin under various planar forces/moments. The code solves for the topology deformation and displacements/rotations at all nodes.

The code allows for:

- Domains of variable dimensions, i.e., overall length and width.
- Beam elements of variable lengths, cross-sections, and material.
- Force(s)/Moment(s) of variable magnitudes and locations.

Code Constraints:

- The minimum length of a beam element is 1/5 of the topology's smallest side length.
- The whole structure should be under static equilibrium.
