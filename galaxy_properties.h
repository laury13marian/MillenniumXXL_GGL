#ifndef GALALXY_PROPERTIES_H
#define GALALXY_PROPERTIES_H

#include "read_data.h"
#include "general.h"

class galaxy_properties {

 private:

  static const int NFILE=4096; //member in the class read_galaxy, to read Raul's files
  const int NB; //nr of bins, luminosity band
  int iBAND;
  vector<double> BinsLum; //lum bins; has dimension NB+1
  vector<size_t> COUNT;
  //vector<vector<int> > GType;
  vector<vector<float> > GCenMvir, GMvir, GDistCen; //GColdGas, GBulgeM, GDiskM, GStellarDiskRad, GBlackHoleM;
  vector<double> AveGCenMvir, AveGMvir, AveGDistCen;
  //vector<double> AveGType, AveGColdGas, AveGBulgeM, AveGDiskM, AveGStellarDiskRad, AveGBlackHoleM;
  read_galaxy OB_RG;
  bool FLAG_READ_DATA;

 public:

  galaxy_properties(const int iBAND_, vector<double> BinsLum_, const int NB_, string path_in, int red);
  void read_and_bin_lum(int ifile);
  void read_and_bin();
  void compute_averages(); 
  void write(string file_out);
  vector<vector<double> > mass_histogram(vector<double> VM);
  void write_histogram(vector<double> VM, string file_out);

};

class galaxy_properties_Millennium {

 private:

  static const int NFILE=512; //member in the class read_galaxy, to read Bruno's files                    
  const int NB; //nr of bins, luminosity band                                           
  int iBAND;
  vector<double> BinsLum; //lum bins; has dimension NB+1                                                         
  vector<size_t> COUNT;
  vector<vector<float> > GCenMvir, GMvir, GDistCen; 
  vector<double> AveGCenMvir, AveGMvir, AveGDistCen;
  read_galaxy_Millennium OB_RG;
  bool FLAG_READ_DATA;

 public:

  galaxy_properties_Millennium(const int iBAND_, vector<double> BinsLum_, const int NB_, 
			       string path_in, const double red);
  void read_and_bin_lum(int ifile);
  void read_and_bin();
  void compute_averages();
  void write(string file_out);
  //vector<vector<double> > mass_histogram(vector<double> VM);
  //void write_histogram(vector<double> VM, string file_out);

};

#endif

