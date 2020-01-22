// cosmological parameters
const double Om=0.267;
const double Ol=0.733;
const double H0=71.9;


// simulation size
const double Lbox = 256.;		//	box length in [Mpc/h]
const int    np   = 512;		//	number of particles to the power -3

char dark_matter_file[]="512x256_03.txt";			//	dark matter catalog




/***************	WHAT DO YOU WUANT TO DO?  ***************************/

  #define non_spherical_search 0	//	0 if you want spherical voids without overlaping
                                        //	1 if you want to follow https://arxiv.org/abs/1905.12450

#if non_spherical_search == 0
//	spherical voids
  const double strong_overlap=0.0;  // it does not matter when you choose only_spherical=1
  const double weak_overlap=1.0;    // >=1 means no overlap at all among voids
  const int only_spherical=1;      // 1: just spherical voids
#endif

#if non_spherical_search == 1
//	non-spherical voids
  const double strong_overlap=0.6;  // overlap fraction to gather spheres into families (linking length for families)
  const double weak_overlap=0.9;    // overlap fraction to create another sphere (reject candidates with overlap stronger (smaller) than this value)
  const int only_spherical=1;      // 0: gather spheres to make non-spherical voids
#endif






/***************	DO YOU WANT TO GO BEYOND https://arxiv.org/abs/1905.12450?  ***************************/

const int   calcula=1;			//	=1 executes steep 1 and 2 described in section 3.2 https://arxiv.org/abs/1905.12450
					//	=0 read file computed previously
const int   no_HD=1;			//	=1 do not improve beyond https://arxiv.org/abs/1905.12450
					//	=0 improves small non-spherical voids beyond the cell size (before steep 3 in section 3.2 https://arxiv.org/abs/1905.12450)
const int   solo_HD=0;			//	steps 1 and 2 were already computed and you want to improve for small non-spherical voids 
const int   ya_catalogo_randoms=0;	//	=1 random sample for finding spheres around small voids (improving initial cell resolution)
const int   ya_radios_randoms=0;	//	=1 if spheres centered in the random catalog were already computed
const int   solo_depura=0;		//	=1 just gather spheres already founded (nothing is done if only_spherical==1)
const int   solo_vol_elip=0;          //	=1 just estimate volume and elipticity for non-spherical voids previously founded

