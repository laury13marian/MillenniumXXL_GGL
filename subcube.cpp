#include "read_data.h"

extern int dummy;
#define WSKIP out.write((char*) &dummy, sizeof(int));

int main(int argc, char *argv[])

{

  if(argc!=2) {
    cerr<<"introduce the snapshot number."<<endl;
    exit(1);
  }
  int dummy;
  const int snap=atoi(argv[1]);
  ostringstream osSNAP;
  osSNAP<<0; osSNAP<<snap;

  const string path_in = "/ptmp/mpa/reangulo/MXXL_Gals/"; 
  int red;                                                             

  if(snap==54) red=24;
  else if(snap==63) red=0;
  else if(snap==42) red=99;
  else {
    cerr<<"wrong snapshot input for Raul's galaxies"<<endl;
    exit(1);
  }                                                                                                     

  const int NFILE=4096;
  //define subcube parameters

  const float Lsub=500.;
  const int Nsub=6;
  const int NT = Nsub*Nsub*Nsub;

  //things to be computed. How many gals for each subcube?
  vector<long long> SS(NT, 0); 
    
  //first determine how many galaxies per subcube we have
  
  read_galaxy Ob(red, path_in); 
  for(int ifile=0; ifile<NFILE; ifile++) {
    Ob.open_subfile(ifile);
    cerr<<"done reading galaxy file "<<ifile<<endl;
    for(long j=0; j<Ob.Ngalaxies[ifile]; j++) {
      int ip[3];
      for(int k=0; k<3; k++) {
	ip[k]=floor(Ob.GalVec[ifile][j].Pos[k]/Lsub);
	//if(ip[k]==Nsub) 
	//ip[k]=0;
      }
      int isub=ip[0]+ip[1]*Nsub+ip[2]*Nsub*Nsub; 
      SS[isub]++;
    }
    Ob.deallocate(ifile);
  }

  string path_out = "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/";
  string file_n = path_out + "subcube_population_snapshot" + osSNAP.str() + ".dat";
  ofstream out(file_n.c_str(), ios::out);
  for(int i=0; i<NT; i++) {
    cerr<<i<<' '<<SS[i]<<endl;
    out<<i<<' '<<SS[i]<<endl;
  }

  out.close();
  cerr<<"done writing file " + file_n<<endl<<endl;*/

  string file_n = "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/subcube_population_snapshot_" 
    + osSNAP.str() + ".dat";

  ifstream in0(file_n.c_str(), ios::in);
  check_stream(in0, file_n);
  for(int i=0; i<NT; i++) {
    int s; long long ns;
    in0>>s>>ns;
    SS[i]=ns;
  }
  in0.close();

  long long NGALTOT=0;
  for(int i=0; i<NT; i++)
    NGALTOT+=SS[i]; //tottal number of galaxies in this snapshot

  galaxy_rec* VecSub[NT];
  vector<long long> Count(NT, 0);
  SScube_header Header[NT];
  string subcube_file[NT];
  vector<string> path_gen(2);                               
  string pcon="GalsLaura";
  path_gen[0] =
    "/ptmp/mpa/rosm/Millennium/Millennium-XXL/subcubes/subcubes_" + osSNAP.str() + '/' + pcon;
  path_gen[1] = "/snapshot_" + osSNAP.str() + "_SubCube_" + pcon + '_';

  for(int i=0; i<NT; i++) {
    ostringstream ooo;
    ooo<<i;
    subcube_file[i]=path_gen[0] + path_gen[1] + ooo.str() + ".0";
    //cerr<<subcube_file[i]<<endl;
  }
  
  //write the header file
  for(int ix=0; ix<Nsub; ix++) 
    for(int iy=0; iy<Nsub; iy++)
      for(int iz=0; iz<Nsub; iz++) {
	int i=ix+iy*Nsub+iz*Nsub*Nsub;
	Header[i].NSubCubeTOT=NT; //nr subcubes in sim
	Header[i].NumberOfFilesPerSubCube=1; //nr output files per subcube
	Header[i].NSubCubePartTOT=NGALTOT; //nr gal per snapshot
	Header[i].NSubCubePart=SS[i]; //nr gal per sub-subcube
	Header[i].time=0.;
	Header[i].redshift = red/100.;
	Header[i].BoxSize=3000.;
	Header[i].Omega0=0.; //no time to fiddle with the next items
	Header[i].OmegaLambda=0.;
	Header[i].HubbleParam=0.;
	Header[i].mass=0.;
	Header[i].LSubCube=Lsub;
	Header[i].Min[0]=ix*Lsub; Header[i].Min[1]=iy*Lsub; Header[i].Min[2]=iz*Lsub;
	Header[i].Max[0]=(ix+1)*Lsub; Header[i].Max[1]=(iy+1)*Lsub; Header[i].Max[2]=(iz+1)*Lsub;
      }
 
  vector<int> split;
  int unit=45;
  int N_unit = floor(NT/unit);
  for(int i=0; i<=N_unit; i++)
    split.push_back(unit*i);
  split.push_back(NT);
  
  read_galaxy Ob(red, path_in);

  for(int k=0; k<split.size()-1; k++) {
    for(int i=split[k]; i<split[k+1]; i++) {
      try {
	VecSub[i]=new galaxy_rec [SS[i]];
      } catch(bad_alloc xa) {
	cerr<<"allocation failure for subcube : "<<i<<endl;
	exit(1);
      }
    }
    for(int ifile=0; ifile<NFILE; ifile++) {
      Ob.open_subfile(ifile);
      cerr<<"done reading galaxy file "<<ifile<<endl;
      for(long j=0; j<Ob.Ngalaxies[ifile]; j++) {
	int ip[3];
	for(int jk=0; jk<3; jk++) {
	  ip[jk]=floor(Ob.GalVec[ifile][j].Pos[jk]/Lsub);
	  if(ip[jk]==Nsub)
	    ip[jk]=0;
	}
	int isub=ip[0]+ip[1]*Nsub+ip[2]*Nsub*Nsub;
	for(int cc=split[k]; cc<split[k+1]; cc++) 
	  if(cc==isub) {
	    VecSub[cc][Count[cc]]=Ob.GalVec[ifile][j];
	    Count[cc]++;
	  }
      }
      Ob.deallocate(ifile);
    }
    cerr<<"compare precomputed sizes of subcube vectors to the currently computed ones: "<<endl;
    for(int cc=split[k]; cc<split[k+1]; cc++)
      cerr<<"subcube "<<cc<<" precomputed size: "<<SS[cc]<<" current size: "<<Count[cc]<<endl;
    //now write subcube files
    for(int cc=split[k]; cc<split[k+1]; cc++) {
      ofstream out(subcube_file[cc].c_str(), ios::out | ios::binary);
      dummy=sizeof(Header[cc]);
      WSKIP;
      out.write((char*) &Header[cc], sizeof(SScube_header));
      dummy=sizeof(Header[cc]);
      WSKIP;
      for(long nk=0; nk<SS[cc]; nk++) {
	dummy=sizeof(galaxy_rec);
	WSKIP;
	out.write((char*) &VecSub[cc][nk], sizeof(galaxy_rec));
	dummy=sizeof(galaxy_rec);
	WSKIP;
      }
      out.close();
      cerr<<"done writing subcube output file for subcube: "<<cc<<endl;
    }
    for(int i=split[k]; i<split[k+1]; i++) {
      delete [] VecSub[i];
      cerr<<"deleting subcube vector: "<<i<<endl;
    }
  }
  

  return 0;

} 
