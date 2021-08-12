#ifndef GALALXY_LUMINOSITY_H
#define GALALXY_LUMINOSITY_H

#include "read_data.h"
#include "general.h"

class galaxy_luminosity {

 private:

  static const int NFILE=4096; //member in the class read_galaxy, to read Raul's files
  double Lmin, Lmax, Mmin, Mmax; //the bounds of the luminosity interval where we do the measurements
  const int NB; //nr of bins, luminosity band
  int iBAND;
  vector<double> binsLL, binsLM, binsML, binsMM; //lum and mass bins
  vector<double> LumCount, MassCount, LumCountT;
  double delL, delM; //logarithmic bin
  read_galaxy OB_RG;

 public:

  galaxy_luminosity(const int iBAND_, const double Lmin_, const double Lmax_, const int NB_, string path_in, int red);
  galaxy_luminosity(const double Mmin_, const double Mmax_, const int NB_, string path_in, int red);
  void read_and_bin_lum(int ifile);
  void read_and_bin_mass(int ifile);
  void read_and_bin(string type); 
  void write(string file_out, string type);

};

class galaxy_luminosity_Millennium {

 private:

  static const int NFILE=512; //member in the class read_galaxy, to read Raul's files                          
  double Lmin, Lmax, Mmin, Mmax; //the bounds of the luminosity interval where we do the measurements                      
  const int NB; //nr of bins, luminosity band                                                                         
  int iBAND;
  vector<double> binsLL, binsLM, binsML, binsMM; //lum and mass bins                                               
  vector<double> LumCount, MassCount, LumCountT;
  double delL, delM; //logarithmic bin                                                                              
  read_galaxy_Millennium OB_RG;

 public:

  galaxy_luminosity(const int iBAND_, const double Lmin_, const double Lmax_, const int NB_, 
		    string path_in, const double_ red);
  //galaxy_luminosity(const double Mmin_, const double Mmax_, const int NB_, string path_in, int red);
  void read_and_bin_lum(int ifile);
  void read_and_bin_mass(int ifile);
  void read_and_bin(string type);
  void write(string file_out, string type);

};

#endif

