#include "galaxy_properties.h"

galaxy_properties::galaxy_properties(const int iBAND_, vector<double> BinsLum_, const int NB_, string path_in, int red):
  BinsLum(BinsLum_), iBAND(iBAND_), NB(NB_), OB_RG(red, path_in)

{
  size_t MAX=250000000;
  GCenMvir.resize(NB); GMvir.resize(NB); GDistCen.resize(NB);  
  COUNT.resize(NB, 0); //number of galaxies in each luminosity bin
  //GType.resize(NB); 
  //GColdGas.resize(NB); GBulgeM.resize(NB); GDiskM.resize(NB);  
  //GStellarDiskRad.resize(NB); GBlackHoleM.resize(NB); 

  for(int ib=2; ib<NB; ib++) {
    //GType[ib].reserve(MAX);
    GCenMvir[ib].reserve(MAX);
    GMvir[ib].reserve(MAX);
    GDistCen[ib].reserve(MAX);
    //GColdGas[ib].reserve(MAX);
    //GBulgeM[ib].reserve(MAX);
    //GDiskM[ib].reserve(MAX);
    //GStellarDiskRad[ib].reserve(MAX);
    //GBlackHoleM[ib].reserve(MAX);
  }

  //AveGType.resize(NB, 0.); 
  AveGCenMvir.resize(NB, 0.); AveGMvir.resize(NB, 0.); AveGDistCen.resize(NB, 0.); 
  //AveGColdGas.resize(NB, 0.); AveGBulgeM.resize(NB, 0.); 
  //AveGDiskM.resize(NB, 0.); AveGStellarDiskRad.resize(NB, 0.); AveGBlackHoleM.resize(NB, 0.);
  FLAG_READ_DATA=false;

}


void galaxy_properties::read_and_bin_lum(int ifile)

{
  //this assumes the bins are negative and going from more negative to less negative
  OB_RG.open_subfile(ifile);
  cerr<<"done reading galaxy file "<<ifile<<endl;
  for(long j=0; j<OB_RG.Ngalaxies[ifile]; j++) {
    for(int ib=0; ib<NB; ib++) {
      if(OB_RG.GalVec[ifile][j].ObsMagDust[iBAND]>BinsLum[ib] && OB_RG.GalVec[ifile][j].ObsMagDust[iBAND]<=BinsLum[ib+1]) {
	if(OB_RG.GalVec[ifile][j].Type==2) { //CAREFUL HERE!!!!!
	COUNT[ib]++;
	GCenMvir[ib].push_back(OB_RG.GalVec[ifile][j].CentralMvir);
	GMvir[ib].push_back(OB_RG.GalVec[ifile][j].Mvir);
	float x[3]; float xdist=0.;
	for(int ipos=0; ipos<3; ipos++) {
	  x[ipos]=OB_RG.GalVec[ifile][j].DistanceToCentralGal[ipos];
	  xdist+=pow(x[ipos], 2);
	}
	GDistCen[ib].push_back(sqrt(xdist));
	//GType[ib].push_back(OB_RG.GalVec[ifile][j].Type);
	//GColdGas[ib].push_back(OB_RG.GalVec[ifile][j].ColdGas);
	//GBulgeM[ib].push_back(OB_RG.GalVec[ifile][j].BulgeMass);
	//GDiskM[ib].push_back(OB_RG.GalVec[ifile][j].DiskMass);
	//GStellarDiskRad[ib].push_back(OB_RG.GalVec[ifile][j].StellarDiskRadius);
	//GBlackHoleM[ib].push_back(OB_RG.GalVec[ifile][j].BlackHoleMass);
	}
      }
    }
  }
  //for(int ib=0; ib<NB; ib++)
  //cerr<<ifile<<' '<<COUNT[ib]<<' '<<GType[ib].size()<<endl;  
  OB_RG.deallocate(ifile);

}


void galaxy_properties::read_and_bin()

{

  for(int ifile=0; ifile<NFILE; ifile++)
      read_and_bin_lum(ifile);
  //check sizes agree with the COUNT vector:
  cerr<<"counted galaxies and the sizes of the vectors: "<<endl;
  for(int ib=0; ib<NB; ib++)
    cerr<<COUNT[ib]<<' '<<GCenMvir[ib].size()<<' '<<GMvir[ib].size()<<' '<<GDistCen[ib].size()<<endl;
  cerr<<endl;
  FLAG_READ_DATA=true;
}


vector<vector<double> > galaxy_properties::mass_histogram(vector<double> VM)

{

  if(!FLAG_READ_DATA) read_and_bin();

  vector<vector<double> > VHist(NB);
  for(int ib=0; ib<NB; ib++) {
    vector<double> TempH(VM.size()-1, 0.);
    for(size_t j=0; j<COUNT[ib]; j++) 
      for(int im=0; im<VM.size()-1; im++) 
	if(VM[im]<GCenMvir[ib][j]*1.e+10 && VM[im+1]>=GCenMvir[ib][j]*1.e+10) 
	  TempH[im]+=1.;
    for(int im=0; im<VM.size()-1; im++) 
      VHist[ib].push_back(TempH[im]);
    TempH.clear();
  }
  
  return VHist;

}


void galaxy_properties::write_histogram(vector<double> VM, string file_out)

{

  vector<vector<double> > Res(NB);
  Res=mass_histogram(VM);
  ofstream out(file_out.c_str(), ios::out);  
  for(int im=0; im<VM.size()-1; im++) {
    out<<VM[im]<<' ';
    for(int ib=0; ib<NB; ib++)
      out<<Res[ib][im]<<' ';
    out<<endl;
  }
  out.close();
  cerr<<"done writing output file " + file_out<<endl;

}


void galaxy_properties::compute_averages()

{

  for(int ib=0; ib<NB; ib++) {
    for(size_t jn=0; jn<COUNT[ib]; jn++) {
      AveGCenMvir[ib]+=GCenMvir[ib][jn];
      AveGMvir[ib]+=GMvir[ib][jn]; 
      AveGDistCen[ib]+=GDistCen[ib][jn];
      /*AveGType[ib]+=GType[ib][jn]; 
      AveGColdGas[ib]+=GColdGas[ib][jn];
      AveGBulgeM[ib]+=GBulgeM[ib][jn];
      AveGDiskM[ib]+=GDiskM[ib][jn];
      AveGStellarDiskRad[ib]+=GStellarDiskRad[ib][jn];
      AveGBlackHoleM[ib]+=GBlackHoleM[ib][jn];*/
    }
  }

  for(int ib=0; ib<NB; ib++) {
    AveGCenMvir[ib]/=COUNT[ib];
    AveGMvir[ib]/=COUNT[ib];
    AveGDistCen[ib]/=COUNT[ib];
    /*AveGType[ib]/=GType[ib].size(); 
    AveGColdGas[ib]/=GType[ib].size();
    AveGBulgeM[ib]/=GType[ib].size();
    AveGDiskM[ib]/=GType[ib].size();
    AveGStellarDiskRad[ib]/=GType[ib].size();
    AveGBlackHoleM[ib]/=GType[ib].size();*/
  }

  /*cerr<<"check that all vectors have the same dimension "<<endl;
  for(int ib=0; ib<NB; ib++)
    cerr<<" bin "<<ib<<' '<<SIZY[ib]<<' '<<GType[ib].size()<<' '<<GCenMvir[ib].size()<<' '<<GMvir[ib].size()
	<<' '<<GDistCen[ib].size()<<' '<<GColdGas[ib].size()<<' '<<GBulgeM[ib].size()<<' '
	<<GDiskM[ib].size()<<' '<<GStellarDiskRad[ib].size()<<' '<<GBlackHoleM[ib].size()<<endl;
	cerr<<endl;*/
  
  //clear the vectors
  
  for(int ib=0; ib<NB; ib++) {
    GCenMvir[ib].clear(); GMvir[ib].clear(); GDistCen[ib].clear();
    //GColdGas[ib].clear(); GBulgeM[ib].clear(); GDiskM[ib].clear(); GStellarDiskRad[ib].clear();
    //GBlackHoleM[ib].clear(); GType[ib].clear(); 
  }
  
}


void galaxy_properties::write(string file_out)

{
  
  read_and_bin();
  compute_averages();
  ofstream out(file_out.c_str(), ios::out);
  for(int ib=0; ib<NB; ib++) {
    out<<BinsLum[ib]<<' '<<BinsLum[ib+1]<<' '<<COUNT[ib]<<' '<<AveGCenMvir[ib]<<' '
       <<AveGMvir[ib]<<' '<<AveGDistCen[ib]<<endl;//' '<<AveGColdGas[ib]<<' '<<AveGBulgeM[ib]<<' '
    //<<AveGDiskM[ib]<<' '<<AveGStellarDiskRad[ib]<<' '<<AveGBlackHoleM[ib]<<endl;

    cerr<<BinsLum[ib]<<' '<<BinsLum[ib+1]<<' '<<COUNT[ib]<<' '<<AveGCenMvir[ib]<<' '
	<<AveGMvir[ib]<<' '<<AveGDistCen[ib]<<endl;//' '<<AveGColdGas[ib]
    //<<' '<<AveGBulgeM[ib]<<' '<<AveGDiskM[ib]<<' '<<AveGStellarDiskRad[ib]<<' '<<AveGBlackHoleM[ib]<<endl;
  }

  out.close();
 
  cerr<<"done writing file " + file_out<<endl;
 
}


//TO DO BRUNO'S CATALOGUE

galaxy_properties_Millennium::galaxy_properties_Millennium(const int iBAND_, vector<double> BinsLum_, const int NB_, 
							   const string path_in, const double red):
  BinsLum(BinsLum_), iBAND(iBAND_), NB(NB_), OB_RG(red, path_in)

{
  size_t MAX=10000000;//250000000;
  GCenMvir.resize(NB); GMvir.resize(NB); GDistCen.resize(NB);
  COUNT.resize(NB, 0); //number of galaxies in each luminosity bin                                  

  for(int ib=2; ib<NB; ib++) {
    GCenMvir[ib].reserve(MAX);
    GMvir[ib].reserve(MAX);
    GDistCen[ib].reserve(MAX);
  }

  AveGCenMvir.resize(NB, 0.); AveGMvir.resize(NB, 0.); AveGDistCen.resize(NB, 0.);
  FLAG_READ_DATA=false;
}


void galaxy_properties_Millennium::read_and_bin_lum(int ifile)

{
  //this assumes the bins are negative and going from more negative to less negative                              
  OB_RG.open_subfile(ifile);
  cerr<<"done reading galaxy file "<<ifile<<endl;
  for(long j=0; j<OB_RG.Ngalaxies[ifile]; j++) {
    for(int ib=0; ib<NB; ib++) {
      //Mag or MagDust, we do both not being sure which one is the right one
      if(OB_RG.GalVec[ifile][j].Mag[iBAND]>BinsLum[ib] && OB_RG.GalVec[ifile][j].Mag[iBAND]<=BinsLum[ib+1]) {
        if(OB_RG.GalVec[ifile][j].Type==2) { //CAREFUL HERE!!!!!                                                  
	  COUNT[ib]++;
	  GCenMvir[ib].push_back(OB_RG.GalVec[ifile][j].CentralMvir);
	  GMvir[ib].push_back(OB_RG.GalVec[ifile][j].Mvir);
	  float x[3]; float xdist=0.;
	  for(int ipos=0; ipos<3; ipos++) {
	    x[ipos]=OB_RG.GalVec[ifile][j].DistanceToCentralGal[ipos];
	    xdist+=pow(x[ipos], 2);
	  }
	  GDistCen[ib].push_back(sqrt(xdist));
	}
      }
    }
  }
  OB_RG.deallocate(ifile);

}


void galaxy_properties_Millennium::read_and_bin()

{

  for(int ifile=0; ifile<NFILE; ifile++)
    read_and_bin_lum(ifile);
  //check sizes agree with the COUNT vector:                                                                      
  cerr<<"counted galaxies and the sizes of the vectors: "<<endl;
  for(int ib=0; ib<NB; ib++)
    cerr<<COUNT[ib]<<' '<<GCenMvir[ib].size()<<' '<<GMvir[ib].size()<<' '<<GDistCen[ib].size()<<endl;
  cerr<<endl;
  FLAG_READ_DATA=true;
}

void galaxy_properties_Millennium::compute_averages()

{

  for(int ib=0; ib<NB; ib++) {
    for(size_t jn=0; jn<COUNT[ib]; jn++) {
      AveGCenMvir[ib]+=GCenMvir[ib][jn];
      AveGMvir[ib]+=GMvir[ib][jn];
      AveGDistCen[ib]+=GDistCen[ib][jn];
    }
  }

  for(int ib=0; ib<NB; ib++) {
    AveGCenMvir[ib]/=COUNT[ib];
    AveGMvir[ib]/=COUNT[ib];
    AveGDistCen[ib]/=COUNT[ib];
  }

  //clear the vectors                                                                                             

  for(int ib=0; ib<NB; ib++) {
    GCenMvir[ib].clear(); GMvir[ib].clear(); GDistCen[ib].clear();
  }

}


void galaxy_properties_Millennium::write(string file_out)

{

  read_and_bin();
  compute_averages();
  ofstream out(file_out.c_str(), ios::out);
  for(int ib=0; ib<NB; ib++) {
    out<<BinsLum[ib]<<' '<<BinsLum[ib+1]<<' '<<COUNT[ib]<<' '<<AveGCenMvir[ib]<<' '
       <<AveGMvir[ib]<<' '<<AveGDistCen[ib]<<endl;

    cerr<<BinsLum[ib]<<' '<<BinsLum[ib+1]<<' '<<COUNT[ib]<<' '<<AveGCenMvir[ib]<<' '
        <<AveGMvir[ib]<<' '<<AveGDistCen[ib]<<endl;
  }
  out.close();
  cerr<<"done writing file " + file_out<<endl;

}
