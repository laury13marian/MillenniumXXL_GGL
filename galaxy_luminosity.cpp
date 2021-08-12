#include "galaxy_luminosity.h"

galaxy_luminosity::galaxy_luminosity(const int iBAND_, const double Lmin_, const double Lmax_, const int NB_, 
				     string path_in, int red):
  Lmin(Lmin_), Lmax(Lmax_), iBAND(iBAND_), NB(NB_), OB_RG(red, path_in)

{

  //first define the luminosity bin vector

  delL=log10(abs(Lmax/Lmin))/(NB-1);
  cerr<<"log bin: "<<delL<<endl<<endl;
  for(int i=0; i<NB; i++) {
    binsLL.push_back(-pow(10, (log10(abs(Lmin))+i*delL)));
    binsLM.push_back(-pow(10, (log10(abs(Lmin))+(i+0.5)*delL)));
  }
  LumCount.resize(NB, 0.);
  //cerr<<"luminosity bins to be used: "<<endl;
  //for(int i=0; i<NB; i++) 
  //cerr<<binsLL[i]<<' '<<binsLM[i]<<endl;
  //cerr<<endl;
  
}


galaxy_luminosity::galaxy_luminosity(const double Mmin_, const double Mmax_, const int NB_,
                                     string path_in, int red):
  Mmin(Mmin_), Mmax(Mmax_), NB(NB_), OB_RG(red, path_in)

{

  //first define the luminosity bin vector

  delM=log10(Mmax/Mmin)/(NB-1);
  cerr<<"stellar mass log bin: "<<delM<<endl<<endl;
  for(int i=0; i<NB; i++) {
    binsML.push_back(pow(10, (log10(Mmin)+i*delM)));
    binsMM.push_back(pow(10, (log10(Mmin)+(i+0.5)*delM)));
  }
  MassCount.resize(NB, 0.);
  cerr<<"stellar mass bins to be used: "<<endl;
  for(int i=0; i<NB; i++)
    cerr<<binsML[i]<<' '<<binsMM[i]<<endl;
  cerr<<endl;

}


void galaxy_luminosity::read_and_bin_lum(int ifile)

{

  double CountMin=0.; 
  double CountMax=0.;
  double CountWrong=0.;
  OB_RG.open_subfile(ifile);
  cerr<<"done reading galaxy file "<<ifile<<endl;
  for(long j=0; j<OB_RG.Ngalaxies[ifile]; j++) {
    if(OB_RG.GalVec[ifile][j].ObsMagDust[iBAND]>0.)
      CountWrong++;
    if(OB_RG.GalVec[ifile][j].ObsMagDust[iBAND]>=Lmin && OB_RG.GalVec[ifile][j].ObsMagDust[iBAND]<0.) {
      CountMin+=1.;
      //cerr<<OB_RG.GalVec[ifile][j].ObsMagDust[iBAND]<<endl;
    }
    else if(OB_RG.GalVec[ifile][j].ObsMagDust[iBAND]<Lmax) {
      CountMax+=1.;
      //cerr<<OB_RG.GalVec[ifile][j].ObsMagDust[iBAND]<<endl;
    }
    else {
      for(int ib=0; ib<NB-1; ib++) {
	if(OB_RG.GalVec[ifile][j].ObsMagDust[iBAND]<binsLL[ib] && OB_RG.GalVec[ifile][j].ObsMagDust[iBAND]>=binsLL[ib+1]) {
	  LumCount[ib]+=1.;
	  //continue;
	}
      }
    }
  }
  LumCount[NB-1]+=CountMax;
  for(int ib=0; ib<NB; ib++)
    cerr<<binsLL[ib]<<' '<<LumCount[ib]<<endl;
  int CHECK=0;
  for(int ib=0; ib<NB; ib++)
    CHECK+=LumCount[ib];
  CHECK+=CountMin;
  CHECK+=CountWrong;
  cerr<<"Total number of galaxies binned: "<<CHECK<<endl;
  cerr<<"wrong gals: "<<CountWrong<<endl;
  OB_RG.deallocate(ifile);

}


void galaxy_luminosity::read_and_bin_mass(int ifile)

{

  double CountMin=0.;
  double CountMax=0.;
  double CountWrong=0.;
  double hfid=0.73;
  double fac=1.e+10/hfid;

  OB_RG.open_subfile(ifile);
  cerr<<"done reading galaxy file "<<ifile<<endl;
  vector<double> TempMass;
  for(long j=0; j<OB_RG.Ngalaxies[ifile]; j++) 
    TempMass.push_back(OB_RG.GalVec[ifile][j].BulgeMass+OB_RG.GalVec[ifile][j].DiskMass);

  for(long j=0; j<TempMass.size(); j++) {
    if(TempMass[j]*fac<Mmin)
      CountMin+=1.;
    else if(TempMass[j]*fac>=Mmax)
      CountMax+=1.;
    else {
      for(int ib=0; ib<NB-1; ib++)
	if(TempMass[j]*fac>=binsML[ib] && TempMass[j]*fac<binsML[ib+1])
	  MassCount[ib]+=1.;
    }
  }
  MassCount[NB-1]+=CountMax;
  TempMass.clear();
  cerr<<"in subfile : "<<ifile<<' ';
  cerr<<"number of low-mass: "<<CountMin<<endl;
  cerr<<"number of high-mass: "<<CountMax<<endl;
  /*for(int ib=0; ib<NB; ib++)
  cerr<<binsML[ib]<<' '<<MassCount[ib]<<endl;
  double CHECK;
  for(int ib=0; ib<NB;ib++)
    CHECK+=MassCount[ib];
  CHECK+=CountMin;

  cerr<<"total number of galaxies binned in mass for subfile: "<<ifile<< "is: "<<CHECK<<endl<<endl;*/

  OB_RG.deallocate(ifile);

}


void galaxy_luminosity::read_and_bin(string type)

{

  if(type=="LUMINOSITY")
    for(int ifile=0; ifile<NFILE; ifile++)
      read_and_bin_lum(ifile);
  else if(type=="MASS")
    for(int ifile=0; ifile<NFILE; ifile++)
      read_and_bin_mass(ifile);

}


void galaxy_luminosity::write(string file_out, string type)

{
  
  read_and_bin(type);
  ofstream out(file_out.c_str(), ios::out);
  for(int ib=0; ib<NB; ib++) {
    if(type=="LUMINOSITY") {
      cerr<<binsLM[ib]<<' '<<LumCount[ib]<<endl;
      out<<binsLL[ib]<<' '<<binsLM[ib]<<' '<<delL<<' '<<LumCount[ib]<<endl;
    }
    else if(type=="MASS") {
      cerr<<binsMM[ib]<<' '<<MassCount[ib]<<endl;
      out<<binsML[ib]<<' '<<binsMM[ib]<<' '<<delM<<' '<<MassCount[ib]<<endl;
    }
  }

  out.close();

}
