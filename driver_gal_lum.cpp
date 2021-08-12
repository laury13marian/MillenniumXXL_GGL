if(snap==54) red=24;
 else if(snap==63) red=0;
  el#include "galaxy_luminosity.h"


int main()

{

  const double Lmin=-9.;
  const double Lmax=-23.;
  string typeL="LUMINOSITY";
  const double Mmin=1.e+6;
  const double Mmax=1.e+12;
  string typeM="MASS";
  const int NB=30;
  const int iMAG=0;
  const int snap=54;
  ostringstream ossMAG;
  ossMAG<<iMAG; 
  const string path_in = "/ptmp/mpa/reangulo/MXXL_Gals/";
  const string path_out = "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/LUMINOSITY/";
  //string file_out = path_out + "galaxy_counts_mag_" + ossMAG.str() + ".dat";
  string file_out = path_out + "galaxy_counts_stellar_mass_RU.dat";
  int red;
  
  if(snap==54) red=24;
  else if(snap==63) red=0;
  else if(snap==42) red=99;
  else {
    cerr<<"wrong snapshot input for Raul's galaxies"<<endl;
    exit(1);
  }
  
  //galaxy_luminosity OB(iMAG, Lmin, Lmax, NB, path_in, red);
  galaxy_luminosity OB(Mmin, Mmax, NB, path_in, red);
  //OB.read_and_bin_mass(304);
  OB.write(file_out, typeM);

  return 0;

}
