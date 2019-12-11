

  //------------------------------------------------------//
  //--  wrapper.h					--//
  //--  Version 2019.01.25				--//
  //--  						--//
  //--  Copyright (C) 2019 - Chieh-An Lin		--//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/	--//
  //------------------------------------------------------//


#ifdef __cplusplus
extern "C" {
#endif

//-- Functions related to a_lm
void readWeight(long long nside, double *weight);
void mapToAlm(long long nside, float *map, int l_maxx, double *alm, double *weight);
void almToMap_spin2(int l_maxx, double *alm, long long nside, float *map1, float *map2);

//-- Functions related to pixel interpolation
void getNgbAndWeight(double pos[2], long long nside, long long neighbor[4], double weight[4]);

#ifdef __cplusplus
}
#endif

