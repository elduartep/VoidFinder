# VoidFinder
Spherical overdensity void finder

This routine is able to find non-spherical voids as described in https://arxiv.org/abs/1905.12450

# parameters
parameters.h

## parameters.h is configured to search for spherical voids without any overlaping among them

## if you want to follow https://arxiv.org/abs/1905.12450 you just need to set
define non_spherical_search 1

you can also change the overlap fraction (linking length) among spheres belonging to the same void by setting 'strong_overlap' to another value

## if you want fo improve the small non-spherical voids (improve resolution) you neet to set
const int   no_HD=0

it improves the small non-spherical voids, wich by default will end up being spherical because of the initial grid resolution
 
