# VoidFinder
Spherical overdensity void finder

This routine is able to find non-spherical voids as described in https://arxiv.org/abs/1905.12450

# Parameters
parameters.h

## parameters.h is configured to search for spherical voids without any overlaping among them

## if you want to follow https://arxiv.org/abs/1905.12450 you just need to set
define non_spherical_search 1

you can also change the overlap fraction (linking length) among spheres belonging to the same void by setting 'strong_overlap' to another value

## if you want fo improve the small non-spherical voids (improve resolution) you neet to set
const int   no_HD=0

it improves the small non-spherical voids, wich by default will end up being spherical because of the initial grid resolution
 
# INPUTS
dark matter catalog in a cosmological box
ASCII fprmat with 3 columns: x, y, z in Mpc/h

# OUTPUTS
several temporary files you may want to save if the code does not finish in the time limit (if you are running it in a quewe)

esferas_... for spherical voids:
ASCII format x, y, z, r in Mpc/h

voids_... is non-spherical voids are requaired:
ASCII format x, y, z, r in Mpc/h

elipticidad_... is non-spherical voids are requaired:
ASCII format r, a, b ,c whre r is the void radious and abc are the semi-axes of the void
