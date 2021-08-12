#include "galaxy_properties.h"

using namespace std;

int main()

{

  const int Nmags=10;
  const int NM=Nmags-1; //number of magnitude/mass bins                                                                  
  vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; mags[5]=-18.;
  mags[6]=-17.; mags[7]=-16.; mags[8]-=15.; mags[9]=-14.;

  /*  const int iMAG=1;
  const int snap=54;
  int red;

  if(snap==54) red=24;
  else if(snap==42) red=99;
  else {
    cerr<<"wrong snapshot input for Raul's galaxies"<<endl;
    exit(1);
  }

  //define mass vector for histogram
  vector<double> MassH;
  const double Mlow=1.e+10;
  const double Msup=3.e+15;
  const int NMass=35;
  const double deltaMass = log(Msup/Mlow)/NMass;
  for(int im=0; im<=NMass; im++)
    MassH.push_back(Mlow*exp(deltaMass*im));
  for(int im=0; im<MassH.size(); im++)
    cerr<<im<<"   "<<MassH[im]<<endl;
  cerr<<endl;

  const string path_in = "/ptmp/mpa/reangulo/MXXL_Gals/";
  const string path_out = "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/GALAXIES/PROPERTIES/";
  //string file_out = path_out + "all_galaxies_per_iMAG_3_properties_type_all.dat";
  //string file_out = path_out + "all_galaxies_per_iMAG_3_properties_type_0.dat";
  const string file_histo = path_out + "histogram_galaxies_FoF_iMAG_1_type_2.dat";
  
  galaxy_properties OB(iMAG, mags, NM, path_in, red);
  vector<vector<double> > Res;
  OB.write_histogram(MassH, file_histo);
  //OB.write(file_out);
  */

  //FOR BRUNO'S CATALOGUE

  const double red=0.24;
  const string path_in = "/galformod/scratch/bmh20/SAM/Henriques2014/snaps/Guo10_bug_fix/MR/";
  const string path_out = "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/GALAXIES/PROPERTIES/"; 
  const int iMAG=17; //corresponding to SDSS red
  string file_out = path_out + "Bruno_galaxies_Mag_per_iMAG_17_properties_type_2.dat"; 
  galaxy_properties_Millennium OB(iMAG, mags, NM, path_in, red);
  OB.write(file_out);

  return 0;

}
