//to read the MillenniumXXL data

#ifndef READ_DATA_H
#define READ_DATA_H

#include <vector>
#include "general.h"
#include "RandomNumbers.h"
#include "tree_code.h"
#include <new>

using namespace std;

extern int dummy;
#define SKIP in.read((char*)&dummy, sizeof(int));

class group_rec {

 public:

  int Len;
  int Offset;
  long long Nr;
  float CM[3];
  float Vel[3];
  float Pos[3];
  float M_Mean200;
  float M_Crit200;
  float M_TopHat200;
  float Vel_Disp;
  float Vel_Disp_Mean200;
  float Vel_Disp_Crit200;
  float Vel_Disp_TopHat200;
  int Nsubs;
  int FirstSub;

};


class subgroup_rec {

 public:
  
  int Len;
  int Offset;
  long long GrNr;
  long long SubNr;
  float Pos[3];
  float Vel[3];
  float CM[3];
  float Spin[3]; 
  float Vel_Disp;
  float Vmax;
  float VmaxRad;
  float HalfmassRad;
  float Shape[6];
  float BindingEnergy;
  float Potential;
  float Profile[9];
  
};


class particle_rec {

 public:

  float Pos[3];
  float Vel[3];
  unsigned long long ID; //specifically for XXL; they could otherwise be unsigned int
};


class particle_data { //for subcubes

 public:

  float pos[3];
};


// Header for the standard file format

class io_header {

 public:

  unsigned int npart[6]; //number of particles of each type in this file 
  double mass[6]; //mass of particles of each type. If 0, then the masses are explicitly
		  //stored in the mass-block of the snapshot file, otherwise they are omitted 
  double time; //time of snapshot file
  double redshift; //redshift of snapshot file
  int flag_sfr;  //flag for star formation included or not
  int flag_feedback;  //flag for feedback (obsolete)
  unsigned int npart_Total[6]; //total number of particles of each type in this snapshot. 
                               //This candiffer from npart for a multi-file snapshot
  int flag_cooling; //flag for cooling  
  int num_subsnapshot;  //number of files in multi-file snapshot 
  double BoxSize;
  double Om;
  double OmLambda; //Hubble parameter in units of 100 km/sec/Mpc
  double Hubble;
  int flag_stellarage; // flag for formation times of star particles 
  int flag_metals; // flags for metallicity values for gas and star
  unsigned int npartTotalHighWord[6]; //High word of the total number of particles of each type
  int flag_entropy_instead_u;   // flag that IC-file contains entropy instead of u 
  int flag_doubleprecision; //flag for double-precision instead of single precision 
  
  int flag_ic_info; 
  float lpt_scalingfactor;  //scaling factor for 2lpt initial conditions 
  char fill[48];  // fills to 256 Bytes 

};



class read_data {

 private:

  static const int Nsnapshots=63; //not sure i need this actually
  static const int Nsubsnapshot=3072;
  vector<group_rec> Group_Vec[Nsubsnapshot];
  vector<subgroup_rec> Subgroup_Vec[Nsubsnapshot];
  particle_rec* Particle_Vec[Nsubsnapshot];
  int snapshot;
  string path_snap, num_snap; //snapshot path
  string convert_snapshot();

 public:

  read_data(string path_snap_, int snapshot_);
  vector<particle> PV[Nsubsnapshot];
  vector<particle> stream_particles(signed long seed, const int nsamp);
  vector<particle> return_formatted_data(signed long seed, const int nsamp);
  void destroy(int subsnap);
  void read_groups(int subsnap);
  void read_ids(int subsnap);
  unsigned int read_particles(int subsnap);
  void group_output(group_rec &T);
  void subgroup_output(subgroup_rec &T);
  void header_output(io_header &T);
  void particle_output(particle_rec &T);

};


#define NMAG 5
class galaxy_rec {

 public:

  int Type; //Galaxy type: 0 for central galaxies of a main halo, 1 for
            //central galaxies in sub-halos, 2 for satellites without halo.
  int Ghost; //the snapshot number where this galaxy was identified 
  float CentralMvir; // 10^10/h Msun virial mass of background (FOF) halo containing this galaxy   
  float Pos[3]; // Mpc/h - Galaxy Positions 
  float Vel[3]; // Mpc/h - Galaxy Positions 
  float Mvir; // 10^10/h Msun - Virial mass of the subhalo the galaxy is/was the center of. 
  float Rvir; // Mpc/h - Virial radius of the subhalo the galaxy is/was the center of. 
  float Vvir; // km/s - Virial velocity of the subhalo the galaxy is/was the center of. 
  float DistanceToCentralGal[3];

  // baryonic reservoirs
  float ColdGas; // 10^10/h Msun - Mass in cold gas.
  float BulgeMass; // 10^10/h Msun - Mass in the bulge
  float DiskMass; // 10^10/h Msun - Mass in hot gas   
  float StellarDiskRadius;
  float BlackHoleMass; // 10^10/h Msun - Mass in black hole
  float Sfr;
  float ObsMagDust[NMAG]; // rest-frame absolute mags //is, rs, us, gs, zs  

};

class galaxy_rec_Millennium {

 public:

  int Type;
  int HaloIndex;
  int SnapNum;
  float LookBackTimeToSnap;
  float CentralMvir;
  float CentralRvir;
  float Pos[3];
  float Vel[3];
  int Len;
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float GasSpin[3];
  float StellarSpin[3];
  float InfallVmax;
  int InfallSnap;
  float InfallHotGas;
  float HotRadius;
  float OriMergTime;
  float MergTime;
  float DistanceToCentralGal[3];
  float ColdGas;
  float BulgeMass;
  float DiskMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
  float BlackHoleGas;
  float ICM;
  float MetalsColdGas;
  float MetalsBulgeMass;
  float MetalsDiskMass;
  float MetalsHotGas;
  float MetalsEjectedMass;
  float MetalsICM;
  float PrimordialAccretionRate;
  float CoolingRate;
  float CoolingRate_beforeAGN;
  float Sfr;
  float SfrBulge;
  float XrayLum;
  float BulgeSize;
  float StellarDiskRadius;
  float GasDiskRadius;
  float CosInclination;
  int DisruptOn;
  int MergeOn;
  float CoolingRadius;
  float QuasarAccretionRate;
  float RadioAccretionRate;
  float Mag[40];
  float MagBulge[40];
  float MagDust[40];
  float MassWeightAge;
  float rbandWeightAge;
  int sfh_ibin;
  float sfh_time[20];
  float sfh_dt[20];
  float sfh_DiskMass[20];
  float sfh_BulgeMass[20];
  float sfh_ICM[20];
  float sfh_MetalsDiskMass[20];
  float sfh_MetalsBulgeMass[20];
  float sfh_MetalsICM[20];
};


//to read the subcube data created by Robert

class SScube_header { //each sub-subcube has in fact this header

 public:

  int NSubCubeTOT; //nr subcubes in sim
  int NumberOfFilesPerSubCube; //nr subsubcubes
  long long NSubCubePartTOT; //nr particles in sim
  long long NSubCubePart; //nr particles in sub-subcube

  double time;         //  time of snapshot file 
  double redshift;     //  redshift of snapshot file
  double BoxSize;      //  box-size of simulation in case periodic boundaries were used 
  double Omega0;       //  matter density in units of critical density 
  double OmegaLambda;  //  cosmological constant density parameter 
  double HubbleParam;  //  Hubble parameter in units of 100 km/sec/Mpc
  double mass;         //  mass of particles 
  double LSubCube; //size subcube
  double Min[3];  //particle with lowest position, actually subcube lower boundaries
  double Max[3]; //particle with highest position, actually subcube upper boundaries

}; 


class load_SScube { //make this for each subcube

 private:


  static const int NSUBCUBE=216; //nr subcubes per sim//6 per cube side
  static const int NFILE=10; //nr files (sub-subcubes) per subcube
  static const int NFILEG=1; //nr of sub-subcubes for the galaxy species (special case)
  const int iSC; // iSSC; //subcube number
  const string SPECIES;
  const signed long rng_seed;
  int nsamp, NFILE0;
  int iMAG; //iMAG determines in which magnitude band we bin the galaxies 
  int NBINS; //number of mass bins for groups or magnitude bins for luminosities
  size_t NTOT; //number of particles in subcube //use this to estimate the shot noise piece
  long long NP_SAMPLED[NFILE]; //number of sampled particles per sub-subcube
  long long NPTOT_SAMPLED; //total number of smapled particles per subcube
  double massPart;
  SScube_header SSH[NFILE]; //CAREFUL, THIS IS ONLY FOR A SUB-SUBCUBE
  string file_SScube[NFILE]; 
  bool FLAG_SAMPLED[NFILE];
  bool FLAG_BINNED;
  vector<particle> FINAL;
  vector<vector<particle> > FINALB;
  vector<int> nsampB;
  long long SampleSize[NFILE]; 
  vector<long long> SampleSizeB[NFILE];

  int dummy;

 public:

  vector<double> binsM;
  particle_data* VecPart[NFILE]; //can also try pointers here and memory allocation   
  group_rec* VecGroup[NFILE];
  vector<group_rec*> VecGroupB[NFILE];
  subgroup_rec* VecSubGroup[NFILE];
  galaxy_rec* VecGal[NFILEG];
  vector<galaxy_rec*> VecGalB[NFILEG];

  load_SScube(const int iSC_, vector<string> pg_, const string SP_, 
	      const int nsamp_, const signed long rng_seed_);
  //define second constructor for mass-binned groups
  load_SScube(const int iSC_, vector<string> pg_, const string SP_, 
	      const signed long rng_seed_, vector<double> binsM_, const int iMAG_);
  void read_species(int sub);
  void read_species(); //read all of the subfiles
  void format_data();
  void deallocate(int sub); 
  vector<particle> return_formatted_data();
  vector<particle> return_formatted_data(int ibin);
  vector<double> lower_bound();
  long long info_output() { 
    cerr<<"total number of particles in subcube "<<iSC<<' '<<NTOT<<endl;
    if(SPECIES=="PARTICLE") 
      cerr<<"total number of sampled particles: "<<NPTOT_SAMPLED<<endl; 
    //works for the particle data, not galaxies or groups
    cerr<<"rate of sampling: "<<nsamp<<endl;
    cerr<<"random seed used: "<<rng_seed<<endl<<endl; 
    return NPTOT_SAMPLED; //same, just particle data, not galaxies or groups
  }

  ~load_SScube();

  void output(SScube_header &T);
  void output(group_rec &T); //group output  
  void output(subgroup_rec &T); //subgroup output 
  void output(particle_data &T);  //particle output
  void output(galaxy_rec &T); //galaxy output

};


//reading galaxy catalogues


class read_galaxy {

 private:
  //path_in: /ptmp/mpa/reangulo/MXXL_Gals
  //SA_zXX_[0-4095] SA_z0.00_999 SA_z0.24_999
  static const int NFILE=4096;
  const int snapshot; 
  string path_in, file_gen;
  bool FLAG_READ[NFILE];

 public:

  read_galaxy(const int snap_, const string path_in_);
  void open_subfile(int sub);
  void output(galaxy_rec &T);
  galaxy_rec* GalVec[NFILE];
  void deallocate(int sub);
  long Ngalaxies[NFILE]; //number of galaxies per file
};

class read_galaxy_Millennium {

 private:
  
  static const int NFILE=512;
  const double red;
  string path_in, file_gen;
  bool FLAG_READ[NFILE];

 public:

  read_galaxy_Millennium(const double red_, const string path_in_);
  void open_subfile(int sub);
  galaxy_rec_Millennium* GalVec[NFILE];
  void deallocate(int sub);
  void output(galaxy_rec_Millennium &T);
  long Ngalaxies[NFILE];
};


#endif
