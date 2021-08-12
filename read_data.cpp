#include "read_data.h"


read_data::read_data(string path_snap_, int snapshot_):path_snap(path_snap_), snapshot(snapshot_)

{

  num_snap=convert_snapshot();


}



string read_data::convert_snapshot()

{
  //snapshot defined in the constructor
  string num;
  ostringstream oss0, oss1;
  oss0<<0;
  oss1<<snapshot;
  if(snapshot >= 0 && snapshot < 10 ) 
    num = oss0.str() + oss0.str() + oss1.str();
  else if(snapshot >= 10 && snapshot < 100)
    num = oss0.str() + oss1.str();
  else
    num = oss1.str();
  
  return num;

}


void read_data::read_groups(int subsnap) //groups and subgroups //start at 0

{
  
  string file_groups;
  int Ngroups, Nids, NTask, Nsubgroups;
  long long TotNgroups, TotNids, TotNsubgroups;
  ostringstream oss;
  oss<<subsnap;

  file_groups = path_snap + 
    "groups_tab/groups_" + num_snap + "/subhalo_tab_" + num_snap + '.' + oss.str();

  ifstream ifs(file_groups.c_str(), ios::binary);
  check_stream(ifs, file_groups);
  
  //read header

  ifs.read((char*) &Ngroups, sizeof(int));
  ifs.read((char*) &TotNgroups, sizeof(long long));
  ifs.read((char*) &Nids, sizeof(int));
  ifs.read((char*) &TotNids, sizeof(long long));
  ifs.read((char*) &NTask, sizeof(int));
  ifs.read((char*) &Nsubgroups, sizeof(int));
  ifs.read((char*) &TotNsubgroups, sizeof(long long));

  cerr<<"Ngroups: "<<Ngroups<<" TotNgroups: "<<TotNgroups<<" Nids: "<<Nids
      <<" TotNids: "<<TotNids<<" NTask: "<<NTask<<" Nsubgroups: "<<Nsubgroups
      <<" TotNsubgroups: "<<TotNsubgroups<<endl<<endl; 

  Group_Vec[subsnap].resize(Ngroups);
  Subgroup_Vec[subsnap].resize(Nsubgroups);

  for(int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].Len, sizeof(int));
  for(int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].Offset, sizeof(int));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].Nr, sizeof(long long));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].CM, 3*sizeof(float));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].Vel, 3*sizeof(float));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].Pos, 3*sizeof(float));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].M_Mean200, sizeof(float));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].M_Crit200, sizeof(float));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].M_TopHat200, sizeof(float));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].Vel_Disp, sizeof(float));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].Vel_Disp_Mean200, sizeof(float));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].Vel_Disp_Crit200, sizeof(float));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].Vel_Disp_TopHat200, sizeof(float));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].Nsubs, sizeof(int));
  for (int i=0; i<Ngroups; i++) 
    ifs.read((char*) &Group_Vec[subsnap][i].FirstSub, sizeof(int));
 
  cerr<<"done reading groups."<<endl<<endl;
  group_output(Group_Vec[subsnap][0]);
  cerr<<endl;
  group_output(Group_Vec[subsnap][99]);

  cerr<<"start reading subgroups"<<endl;

  for(int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].Len, sizeof(int));
  for(int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].Offset, sizeof(int));
  for (int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].GrNr, sizeof(long long));
  for (int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].SubNr, sizeof(long long));
  for (int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].Pos, 3*sizeof(float));
  for (int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].Vel, 3*sizeof(float));
  for (int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].CM, 3*sizeof(float));
  for (int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].Spin, 3*sizeof(float));
  for(int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].Vel_Disp, sizeof(float));
  for(int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].Vmax, sizeof(float));
  for(int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].VmaxRad, sizeof(float));
  for(int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].HalfmassRad, sizeof(float));
  for (int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].Shape, 6*sizeof(float));
  for(int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].BindingEnergy, sizeof(float));
  for(int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].Potential, sizeof(float));
  for (int i=0; i<Nsubgroups; i++) 
    ifs.read((char*) &Subgroup_Vec[subsnap][i].Profile, 9*sizeof(float));
  
  subgroup_output(Subgroup_Vec[subsnap][0]);
  subgroup_output(Subgroup_Vec[subsnap][99]);
 
  ifs.close();
  cerr<<"done reading group and subgroup  file " + file_groups<<endl;
  
}


void read_data::group_output(group_rec &T)

{

  cerr<<"outputting group "<<T.Nr<<endl;
  cerr<<"length: "<<T.Len<<"; offset: "<<T.Offset<<"; Nr: "<<T.Nr<<endl;
  cerr<<"CM position: "<<T.CM[0]<<' '<<T.CM[1]<<' '<<T.CM[2]<<endl;
  cerr<<"velocity: "<<T.Vel[0]<<' '<<T.Vel[1]<<' '<<T.Vel[2]<<endl;
  cerr<<"potential center: "<<T.Pos[0]<<' '<<T.Pos[1]<<' '<<T.Pos[2]<<endl;
  cerr<<"M_mean200: "<<T.M_Mean200<<"; M_crit200: "<<T.M_Crit200<<"; M_tophat200: "
      <<T.M_TopHat200<<endl;
  cerr<<"Vel: "<<T.Vel_Disp<<"; Vel_mean: "<<T.Vel_Disp_Mean200<<"; Vel_crit: "
      <<T.Vel_Disp_Crit200<<"; Vel_TH: "<<T.Vel_Disp_TopHat200<<endl;
  cerr<<" Nsubs: "<<T.Nsubs<<"; first sub: "<<T.FirstSub<<endl;

  cerr<<"done outputting group "<<T.Nr<<endl<<endl;

}


void read_data::subgroup_output(subgroup_rec &T)

{

  cerr<<"outputting subgroup "<<T.GrNr<<endl;
  cerr<<"length (nr particles): "<<T.Len<<"; offset: "<<T.Offset<<"; Nr: "<<T.GrNr<<endl;
  cerr<<"potential center: "<<T.Pos[0]<<' '<<T.Pos[1]<<' '<<T.Pos[2]<<endl;
  cerr<<"velocity: "<<T.Vel[0]<<' '<<T.Vel[1]<<' '<<T.Vel[2]<<endl;
  cerr<<"CM position: "<<T.CM[0]<<' '<<T.CM[1]<<' '<<T.CM[2]<<endl;
  cerr<<"spin: "<<T.Spin[0]<<' '<<T.Spin[1]<<' '<<T.Spin[2]<<endl;
  cerr<<"Vel: "<<T.Vel_Disp<<"; Vmax "<<T.Vmax<<"; VmaxRad: "<<T.VmaxRad
      <<"; HalfmassRad: "<<T.HalfmassRad<<endl;

  cerr<<"done outputting subgroup "<<T.GrNr<<endl<<endl;

}


void read_data::particle_output(particle_rec &T)

{

  cerr<<"position: "<<T.Pos[0]<<' '<<T.Pos[1]<<' '<<T.Pos[2]<<endl;
  cerr<<"velocity: "<<T.Vel[0]<<' '<<T.Vel[1]<<' '<<T.Vel[2]<<endl;
  cerr<<"ID: "<<T.ID<<endl<<endl;

}


void read_data::header_output(io_header &T)

{

  cerr<<"outputting header file from the particle data."<<endl;
  cerr<<"redshift: "<<T.redshift<<"; time: "<<T.time<<"; BoxSize: "<<T.BoxSize<<endl;
  cerr<<"number of particles: "<<T.npart[0]<<' '<<T.npart[1]<<' '<<T.npart[2]<<' '
      <<T.npart[3]<<' '<<T.npart[4]<<' '<<T.npart[5]<<endl;
  cerr<<"mass: "<<T.mass[0]<<' '<<T.mass[1]<<' '<<T.mass[2]<<' '<<T.mass[3]<<' '
      <<T.mass[4]<<' '<<T.mass[5]<<endl;
  cerr<<"Nsnapshots: "<<T.num_subsnapshot<<endl;
  cerr<<"Om: "<<T.Om<<"; OmLambda: "<<T.OmLambda<<"; Hubble: "<<T.Hubble<<endl;
  cerr<<"done outputting header file."<<endl<<endl;

}


void read_data::destroy(int subsnap)

{

  if(Particle_Vec[subsnap]!=0) {
    delete [] Particle_Vec[subsnap];
    cerr<<"deleting Particle_Vec for subsnapshot: "<<subsnap<<endl;
  }
  if(Group_Vec[subsnap].size()>0) {
    Group_Vec[subsnap].clear();
    cerr<<"cleared Group_Vec for subsnapshot: "<<subsnap<<endl;
  }
  if(Subgroup_Vec[subsnap].size()>0) {
    Subgroup_Vec[subsnap].clear();
    cerr<<"cleared Subgroup_Vec for subsnapshot: "<<subsnap<<endl;
  }
}


unsigned int read_data::read_particles(int subsnap)

{

  //these files available 27, 41, 54 63
  ostringstream oss;
  oss<<subsnap;
  io_header header;
  particle_rec temp;
  int dummy;
  string file_snapshot = path_snap 
    + "/snapshots/snapdir_" + num_snap + "/snapshot_" + num_snap + '.' + oss.str(); 

  ifstream in(file_snapshot.c_str(), ios::binary);
  check_stream(in, file_snapshot);
  SKIP;
  in.read((char*) &header, sizeof(io_header));
  SKIP;
  header_output(header);
  unsigned int NumPart=header.npart[1];//WHY NOT 0?
  try {
    Particle_Vec[subsnap] = new particle_rec [NumPart];
  } catch(bad_alloc xa) {
    cerr<<"allocation failure for particles in subsnapshot: "<<subsnap<<endl;
   return 1;
  }
  
  SKIP;
  for(unsigned int i=0; i<NumPart; i++) 
    in.read((char*) &Particle_Vec[subsnap][i].Pos, 3*sizeof(float));
  SKIP;
  SKIP;
  for(unsigned int i=0; i<NumPart; i++) 
    in.read((char*) &Particle_Vec[subsnap][i].Vel, 3*sizeof(float));
  SKIP;
  SKIP;
  for(unsigned int i=0; i<NumPart; i++) 
    in.read((char*) &Particle_Vec[subsnap][i].ID, sizeof(unsigned long long));
  SKIP;
  cerr<<"done reading particle data from " + file_snapshot<<endl;
  particle_output(Particle_Vec[subsnap][0]);  

  in.close();

  return NumPart;

}

vector<particle> read_data::stream_particles(const signed long seed, const int nsamp)

{
  random2 OBran(seed);
  vector<particle> VV;
  //VV.reserve(3000000);
  for(int s=0; s<Nsubsnapshot; s++) {
  //for(int s=1500; s<1600; s++) {
    int unsigned Num;
    if(Particle_Vec[s]==0) 
      Num=read_particles(s);
    for(register size_t i=0; i<Num; i++) {
      if(OBran() <= 1./nsamp) {   
	vector<float> temp(3);
	for(int j=0; j<3; j++)
	  temp[j]=Particle_Vec[s][i].Pos[j];
	particle OT(temp);
	PV[s].push_back(OT);
      }
    }
    cerr<<"size of and PV for subsnap: "<<s<<" are: "<<PV[s].size()<<endl;
    cerr<<PV[s][0].pos[0]<<' '<<PV[s][0].pos[1]<<' '<<PV[s][0].pos[2]<<endl; 
    destroy(s);
    for(long j=0; j<PV[s].size(); j++)
      VV.push_back(PV[s][j]);
    PV[s].clear();
  }
  
  return VV;
  
}


/*vector<particle> read_data::return_formatted_data(const signed long seed, const int nsamp)

{
 
  vector<particle> VV;
  random2 OBran(seed);
  for(int i=0; i<2; i++) {
    stream_particles(i);
    for(register size_t k=0; k<PV[i].size(); k++) {
      if(OBran() <= 1./nsamp)
	VV.push_back(PV[i][k]);
    }
    PV[i].clear();
  }
  
  cerr<<"size of sampled particle data: "<<VV.size()<<endl;
  
  return VV;
  }*/


/*void read_data::read_ids(int subsnap)

{

  //the only snapshots with these files are: 0, 1, 2, 3, 9, 10, 11, 27, 54, 63
  ostringstream oss;
  oss<<subsnap;
  string filegroup_ids = path_snap + 
    "groups_ids/groups_ids_" + num_snap + "/subhalo_ids_" + num_snap + '.' + oss.str();
  long Ngroups, Nids, NTask, Need;
  long long TotNgroups, TotNids, NidsPrevious;
  vector<unsigned long> SubhaloIDs_x(Need), SubhaloIDs(Nids);

  ifstream in(filegroup_ids.c_str(), ios::binary);
  check_stream(in, filegroup_ids);
  SKIP;
  in.read((char*) &Ngroups, sizeof(long));
  SKIP;
  in.read((char*) &TotNgroups, sizeof(long long));
  SKIP;
  in.read((char*) &Nids, sizeof(long));
  SKIP;
  in.read((char*) &TotNids, sizeof(long long));
  SKIP;
  in.read((char*) &NTask, sizeof(long));
  SKIP;
  in.read((char*) &NidsPrevious, sizeof(long long));
  SKIP;

  if(FLAG_COMPRESS) {
    in.read((char*) &Need, sizeof(long));
    SKIP;
    in.read((char*) &SubhaloIDs_x, sizeof(Need*unsigned long));
    SKIP;
    //something else here, uncompress
  }
  else {
    in.read((char*) &SubhaloIDs, sizeof(Nids*unsigned long));
    SKIP;
  }
}*/



 //SUBCUBES CODE



load_SScube::load_SScube(const int iSC_, vector<string> pg, const string SP_, 
			 const int nsamp_, const signed long rng_seed_):
  iSC(iSC_), SPECIES(SP_), nsamp(nsamp_), rng_seed(rng_seed_)

{

  if( iSC<0 || iSC>=NSUBCUBE ) {
    cerr<<"subcube number out of bounds."<<endl;
    exit(10);
  }
  FLAG_BINNED=false;
  cerr<<"random seed number: "<<rng_seed<<endl<<endl; 
  string path_spec;
  NTOT=0.;
  if(SPECIES=="PARTICLE") 
    path_spec = "Part";
  else if(SPECIES=="GROUP")
    path_spec = "Grp";
  else if(SPECIES=="SUBGROUP")
    path_spec = "SubGrp";
  else if(SPECIES=="GALAXY")
    path_spec = "GalsLaura";
  else {
    cerr<<"SPECIES can be only 'PARTICLE', 'GROUP', 'SUBGROUP', 'GALAXY' "<<endl;
    exit(10);
  }

  ostringstream oss0;
  oss0<<iSC;
  if(SPECIES=="GALAXY")
    NFILE0=NFILEG;
  else
    NFILE0=NFILE;
  for(int iSSC=0; iSSC<NFILE0; iSSC++) {
    ostringstream oss1;
    oss1<<iSSC;
    file_SScube[iSSC] = pg[0] + path_spec + pg[1] + path_spec + '_'+ oss0.str() + '.' + oss1.str(); 
    FLAG_SAMPLED[iSSC]=false;
    //to read the general header, whatever the species: particles, groups, subgroups
    ifstream in(file_SScube[iSSC].c_str(), ios::binary);
    check_stream(in, file_SScube[iSSC]);
    //read header
  
    SKIP;
    in.read((char*) &SSH[iSSC], sizeof(SScube_header));
    SKIP;
    NTOT+=SSH[iSSC].NSubCubePart;
    in.close();
    //output(SSH[iSSC]);
  }
  output(SSH[0]);
  cerr<<endl;
  //cerr<<"critical density of the Universe: " //valid only for particles
  //<<pow(SSH[0].BoxSize, -3)* SSH[0].NSubCubePartTOT * SSH[0].mass *1.e+10/SSH[0].Omega0<<endl;

}


load_SScube::load_SScube(const int iSC_, vector<string> pg, const string SP_, 
			 const signed long rng_seed_, vector<double> binsM_, const int iMAG_):
  iSC(iSC_), SPECIES(SP_), rng_seed(rng_seed_), binsM(binsM_), iMAG(iMAG_) 
 
{

  FLAG_BINNED=true;
  NBINS=binsM.size()-1; //the assumption is that binsM[binsM.size()-1] is an upper limit on the mass or magnitude 

  if(SPECIES=="GALAXY") {
    NFILE0=NFILEG;
    for(int i=0; i<NFILEG; i++) {
      VecGalB[i].resize(NBINS);
      SampleSizeB[i].resize(NBINS);
    }
  }
  else if(SPECIES=="GROUP") {
    NFILE0=NFILE;
    for(int i=0; i<NFILE; i++) {
      VecGroupB[i].resize(NBINS);
      SampleSizeB[i].resize(NBINS);
    }
  }
  
  FINALB.resize(NBINS);

  if( iSC<0 || iSC>=NSUBCUBE ) {
    cerr<<"subcube number out of bounds."<<endl;
    exit(10);
  }
  cerr<<"random seed number: "<<rng_seed<<endl<<endl; 
  string path_spec;
  NTOT=0.;
  if(SPECIES=="GROUP")
    path_spec = "Grp";
  else if(SPECIES=="SUBGROUP")
    path_spec = "SubGrp";
  else if(SPECIES=="GALAXY")
    path_spec = "GalsLaura";
  else {
    cerr<<"SPECIES can be only 'GROUP', 'SUBGROUP'"<<endl;
    exit(10);
  }
  
  ostringstream oss0;
  oss0<<iSC;
  
  for(int iSSC=0; iSSC<NFILE0; iSSC++) {
    ostringstream oss1;
    oss1<<iSSC;
    file_SScube[iSSC] = pg[0] + path_spec + pg[1] + path_spec + '_'+ oss0.str() + '.' + oss1.str(); 
    FLAG_SAMPLED[iSSC]=false;
    //to read the general header, whatever the species: particles, groups, subgroups
    ifstream in(file_SScube[iSSC].c_str(), ios::binary);
    check_stream(in, file_SScube[iSSC]);
    //read header
  
    SKIP;
    in.read((char*) &SSH[iSSC], sizeof(SScube_header));
    SKIP;
    NTOT+=SSH[iSSC].NSubCubePart;
    in.close();
    //output(SSH[iSSC]);
  }
  output(SSH[0]);
  cerr<<endl;

}


vector<double> load_SScube::lower_bound() 
  
{
  
  vector<double> Bound(3);
  for(int i=0; i<3; i++)
    Bound[i]=SSH[0].Min[i];
  
  return Bound;
  
}



void load_SScube::read_species(int sub) //0<=sub<NFILE; runs through sub-subcubes

{

  int dummy;
  ifstream in(file_SScube[sub].c_str(), ios::binary);
  check_stream(in, file_SScube[sub]);
  //read header
  SKIP;
  cerr<<"DUMMY AFTER FIRST SKIP IN HEADER "<<dummy<<endl;
  in.read((char*) &SSH[sub], sizeof(SScube_header));
  SKIP;
  cerr<<"DUMMY AFTER SECOND SKIP IN HEADER "<<dummy<<endl;
  random2 OBran(rng_seed);
  //read particle species

  long long NPart=SSH[sub].NSubCubePart;
  cerr<<"number of objects of the type "<<SPECIES<<" in sub-subcube "<<sub<<" of subcube "
      <<iSC<<" is: "<<NPart<<endl;
  if(SPECIES=="PARTICLE") {
    particle_data* Temp;
    particle_data* TempS;
    try {
      Temp=new particle_data [NPart];
    } catch(bad_alloc xa) {
	cerr<<"STEP 1: allocation failure for temp particles in subcube: "
	    <<iSC<<" and subsubcube: "<<sub<<endl;
	exit(1);
    }
    SKIP;
    for(long long i=0; i<NPart; i++)
      in.read((char*) &Temp[i], sizeof(class particle_data));
    SKIP;  

    long long NPartT=static_cast<long long>(NPart/nsamp *(1+0.1)); //guess how big the random sample is
    try {
      TempS=new particle_data [NPartT];
    } catch(bad_alloc xa) {
      cerr<<"STEP 2: allocation failure for temp particles in subcube: "
	  <<iSC<<" and subsubcube: "<<sub<<endl;
      exit(1);
    }
    long long NPartS=0;  //actual number of sampled particles
    for(long long j=0; j<NPart; j++)
      if(OBran() <= 1./nsamp) {
	TempS[NPartS]=Temp[j]; 
	NPartS++;
      }
    delete [] Temp;
    cerr<<"size of sampled particle vector: "<<NPartS<<endl;

    try {
      VecPart[sub]=new particle_data [NPartS];
    } catch(bad_alloc xa) {
      cerr<<"STEP 3: allocation failure for temp particles in subcube: "
	  <<iSC<<" and subsubcube: "<<sub<<endl;
      exit(1);
    }
    for(long long j=0; j<NPartS; j++)
      VecPart[sub][j]=TempS[j];
    NP_SAMPLED[sub]=NPartS;
    delete [] TempS;
    cerr<<"first particle: "; output(VecPart[sub][0]); 
    cerr<<"last particle: "; output(VecPart[sub][NPartS-1]);
    SampleSize[sub]=NPartS;
    FLAG_SAMPLED[sub]=true;
   
  }


  else if(SPECIES=="GALAXY") { 
    if(!FLAG_BINNED) {
      try {
	VecGal[sub]=new galaxy_rec [NPart];
      } catch(bad_alloc xa) {
	cerr<<"allocation failure for temp galaxies in subcube: "<<iSC<<" and subsubcube: "<<sub<<endl;
	exit(1);
      }
    for(long long i=0; i<NPart; i++) {
      SKIP;
      in.read((char*) &VecGal[sub][i], sizeof(galaxy_rec));
      SKIP;
    }
    for(int i=0; i<5  ; i++) {
      output(VecGal[sub][i]);    
      cerr<<endl;                    
    }            
    SampleSize[sub]=NPart;
    FLAG_SAMPLED[sub]=true;
    }

    else {
      galaxy_rec* Temp;
      galaxy_rec* TempB[NBINS];
      galaxy_rec* TempB2[NBINS];

      try {
        Temp=new galaxy_rec [NPart];
      } catch(bad_alloc xa) {
        cerr<<"STEP 1:allocation failure for temp galaxies in subcube: "<<iSC<<" and subsubcube: "<<sub<<endl;
        exit(1);
      }
      
      for(long long i=0; i<NPart; i++) {
        SKIP;
        in.read((char*) &Temp[i], sizeof(class galaxy_rec));
        SKIP;
      }
      
      vector<long long> NPartT(NBINS, 0); //see how many objects per bin
      for(long long j=0; j<NPart; j++)
        for(int ib=0; ib<NBINS; ib++)
          if(binsM[ib]<=Temp[j].ObsMagDust[iMAG] && binsM[ib+1]>Temp[j].ObsMagDust[iMAG])
            NPartT[ib]++;

      long long check_bins=0;
      for(int ib=0; ib<NBINS; ib++)
        check_bins+=NPartT[ib];
      if(check_bins != NPart) {
	cerr<<"THE GALAXY BINNING DOES NOT ADD UP for sub-subcube "<<sub<<" of subcube "<<iSC
            <<" . Original size: "<<NPart<<" binned size: "<<check_bins<<endl<<endl;
      }

      cerr<<"original number of bins: "<<endl;
      for(int ib=0; ib<NBINS; ib++)
	cerr<<binsM[ib]<<' '<<binsM[ib+1]<<' '<<NPartT[ib]<<endl;

      for(int ib=0; ib<NBINS; ib++) { //define binned arrays
	try {
	  TempB[ib]=new galaxy_rec [NPartT[ib]];
	} catch(bad_alloc xa) {
	  cerr<<"STEP 2:allocation failure for binned galaxies in subcube: "<<iSC
	      <<" and subsubcube: "<<sub<<endl;
	}
      }
      
      vector<long long> county(NBINS, 0);
      for(register long long j=0; j<NPart; j++)
	for(int ib=0; ib<NBINS; ib++) 
          if(binsM[ib]<=Temp[j].ObsMagDust[iMAG] && binsM[ib+1]>Temp[j].ObsMagDust[iMAG]) {
	    TempB[ib][county[ib]]=Temp[j];
	    county[ib]++;
	  }
      delete [] Temp;
      county.clear();
      
      const long Nref=150000; //we want to keep the galaxy vectors at this level
      
      for(int ib=0; ib<NBINS; ib++) {
	long Neff;
	if(NPartT[ib]<Nref) 
	  Neff=NPartT[ib];
       	else 
	  Neff=(1+0.1)*Nref;
	try {
	  TempB2[ib]=new galaxy_rec [Neff];
	} catch(bad_alloc xa) {
	  cerr<<"STEP 3 :allocation failure for binned galaxies in subcube: "<<iSC
	      <<" , subsubcube: "<<sub<<" and magnitude bin: "<<ib<<endl;
	  exit(1);
	}
      }
      
      vector<long long> county2(NBINS, 0);
      for(int ib=0; ib<NBINS; ib++) {
	if(NPartT[ib]<Nref) {
	  for(long j=0; j<NPartT[ib]; j++) {
	    TempB2[ib][j]=TempB[ib][j];
	    county2[ib]++;
	  }
	}
	else {
	  float rest = (Nref+0.)/NPartT[ib];
	  const signed long seed_bin = rng_seed + ib*31;
	  random2 OBranGal(seed_bin);
	  for(register long j=0; j<NPartT[ib]; j++) 
	    if(OBranGal() <= rest) {
	      TempB2[ib][county2[ib]]=TempB[ib][j];
	      county2[ib]++;
	    }
	}
	delete [] TempB[ib];
      }

      for(int ib=0; ib<NBINS; ib++) {
	try {
	  VecGalB[sub][ib] = new galaxy_rec [county2[ib]];
	} catch(bad_alloc xa) {
          cerr<<"STEP 4 :allocation failure for binned galaxies in subcube: "<<iSC
              <<" , subsubcube: "<<sub<<" and magnitude bin: "<<ib<<endl;
          exit(1);
        }
      }
      
      for(int ib=0; ib<NBINS; ib++) {
	for(long long j=0; j<county2[ib]; j++)
	  VecGalB[sub][ib][j] = TempB2[ib][j]; 
	delete [] TempB2[ib];
	SampleSizeB[sub][ib]=county2[ib];      
      }
      
      cerr<<"size of binned galaxy vector for subcube "<<iSC<<" and subsubcube: "<<sub<<endl;
      for(int ib=0; ib<NBINS; ib++)
	cerr<<ib<<' '<<binsM[ib]<<' '<<binsM[ib+1]<<' '<<SampleSizeB[sub][ib]<<endl;
      FLAG_SAMPLED[sub]=true;
    }
  }
  
  else if(SPECIES=="GROUP") {
    if(!FLAG_BINNED) {
      group_rec* Temp;
      group_rec* TempS;
      
      try {
	Temp=new group_rec [NPart];
      } catch(bad_alloc xa) {
	cerr<<"allocation failure for temp groups in subcube: "<<iSC<<" and subsubcube: "<<sub<<endl;
	exit(1);
      }
      
      for(long long i=0; i<NPart; i++) {
	SKIP;
	in.read((char*) &Temp[i], sizeof(class group_rec));
	SKIP;
      }
              
      //alocate memory to sampled array;

      long long NPartT = static_cast<long long>(NPart/nsamp *(1+0.1)); //temporary size, just a guess
      cerr<<NPart<<' '<<NPartT<<endl;
      
      try {
	TempS=new group_rec [NPartT];
      } catch(bad_alloc xa) {
	cerr<<"allocation failure for temporray sampled groups in subcube: "
	    <<iSC<<" and subsubcube: "<<sub<<endl;
	exit(1);
      }
      
      long long NPartS = 0; //real size of sample vector
      for(long long j=0; j<NPart; j++) {
	if(OBran() <= 1./nsamp) {//the non binning part	
	  TempS[NPartS]=Temp[j];
	  NPartS++;
	}
	if(NPartS>=NPartT) {
	  cerr<<"temporary sampled group array for subcube "<<iSC
	      <<" and sub-subcube "<<sub<<" is too large"<<endl;
	  exit(1);
	}
      }
      delete [] Temp;
      cerr<<"size of sampled vector: "<<NPartS<<endl;
      
      try {
	VecGroup[sub]=new group_rec [NPartS];
      } catch(bad_alloc xa) {
	cerr<<"allocation failure for sampled groups in subcube: "<<iSC<<" and subsubcube: "
	    <<sub<<endl;
	exit(1);
      }
      
      for(long long j=0; j<NPartS; j++)
	VecGroup[sub][j]=TempS[j];
      
      delete [] TempS;
    
      SampleSize[sub] = NPartS; //write down size
      FLAG_SAMPLED[sub]=true;
    }
    
    else  {
      group_rec* Temp;
      group_rec* TempB[NBINS]; 
      group_rec* TempB2[NBINS];
      try {
	Temp=new group_rec [NPart];
      } catch(bad_alloc xa) {
	cerr<<"allocation failure for temp groups in subcube: "<<iSC<<" and subsubcube: "<<sub<<endl;
	exit(1);
      }
      
      for(long long i=0; i<NPart; i++) {
	SKIP;
	in.read((char*) &Temp[i], sizeof(class group_rec));
	SKIP;
      }
    
      vector<long long> NPartT(NBINS, 0), NPartT2(NBINS), NPartS(NBINS, 0);
      for(long long j=0; j<NPart; j++) 
	for(int ib=0; ib<NBINS; ib++) 
	  if(binsM[ib]<=Temp[j].Len*massPart && binsM[ib+1]>Temp[j].Len*massPart)
	    NPartT[ib]++;
	
      long long check_bins=0;
      for(int ib=0; ib<NBINS; ib++)
	check_bins+=NPartT[ib];
      if(check_bins != NPart) {
	cerr<<"THE BINNING DOES NOT ADD UP for sub-subcube "<<sub<<" of subcube "<<iSC
	    <<" . Original size: "<<NPart<<" binned size: "<<check_bins<<endl<<endl;
	//exit(1);
      }
      
      for(int ib=0; ib<NBINS; ib++) {
	try {
	  TempB[ib]=new group_rec [NPartT[ib]];
	} catch(bad_alloc xa) {
	  cerr<<"STEP 2:allocation failure for sampled binned groups in subcube: "<<iSC
	      <<" and subsubcube: "<<sub<<endl;
	  exit(1);
	}
      }
      vector<long long> county(NBINS, 0);
      for(long long j=0; j<NPart; j++) 
	for(int ib=0; ib<NBINS; ib++) 
	  if(binsM[ib]<=Temp[j].Len*massPart && binsM[ib+1]>Temp[j].Len*massPart) {
	    TempB[ib][county[ib]]=Temp[j];
	    county[ib]++;
	  }
      delete [] Temp;

      for(int ib=0; ib<NBINS; ib++) {
	NPartT2[ib]=static_cast<long long>(NPartT[ib]/nsampB[ib]*(1+0.1));
      	try {
	  TempB2[ib]=new group_rec [NPartT2[ib]];
	} catch(bad_alloc xa) {
	  cerr<<"STEP 3: allocation failure for sampled binned groups in subcube: "<<iSC
	      <<" and subsubcube: "<<sub<<endl;
	  exit(1);
	}

	for(long long j=0; j<NPartT[ib]; j++)  
	  if(OBran() <= 1./nsampB[ib]) {
	    TempB2[ib][NPartS[ib]]=TempB[ib][j];
	    NPartS[ib]++; 
	  }
	delete [] TempB[ib];
	
	try {
	  VecGroupB[sub][ib]=new group_rec [NPartS[ib]];
	} catch(bad_alloc xa) {
	  cerr<<"STEP 4: allocation failure for sampled binned groups in subcube: "<<iSC
	      <<" and subsubcube: "<<sub<<endl;
	  exit(1);
	}
	for(long long j=0; j<NPartS[ib]; j++) 
	  VecGroupB[sub][ib][j]=TempB2[ib][j];
	delete [] TempB2[ib];
	SampleSizeB[sub][ib]=NPartS[ib];

      }
      cerr<<"size of binned and sampled group vector for subsubcube: "<<sub<<endl;
      for(int ib=0; ib<NBINS; ib++) 
	cerr<<ib<<' '<<binsM[ib]<<' '<<SampleSizeB[sub][ib]<<endl;
      FLAG_SAMPLED[sub]=true;
    }       
  }
  
  else if(SPECIES=="SUBGROUP") {
    subgroup_rec* Temp;
    subgroup_rec* TempS;
    
    try {
      Temp=new subgroup_rec [NPart];
    } catch(bad_alloc xa) {
      cerr<<"allocation failure for temp groups in subcube: "<<iSC<<" and subsubcube: "<<sub<<endl;
      exit(1);
    }
    
    for(long long i=0; i<NPart; i++) {
      SKIP;
      in.read((char*) &Temp[i], sizeof(class subgroup_rec));
      SKIP;
    }
    
    long long NPartT = static_cast<long long>(NPart/nsamp *(1+0.1)); //temporary size, just a guess
    cerr<<NPart<<' '<<NPartT<<endl;
    
    try {
      TempS=new subgroup_rec [NPartT];
    } catch(bad_alloc xa) {
	cerr<<"allocation failure for temporray sampled groups in subcube: "
	    <<iSC<<" and subsubcube: "<<sub<<endl;
	exit(1);
      }
    
    long long NPartS = 0; //real size of sample vector
    for(long long j=0; j<NPart; j++) {
      if(OBran() <= 1./nsamp) {//the non binning part	
	TempS[NPartS]=Temp[j];
	NPartS++;
      }
      if(NPartS>=NPartT) {
	cerr<<"temporary sampled subgroup array for subcube "<<iSC
	    <<" and sub-subcube "<<sub<<" is too large"<<endl;
	exit(1);
      }
    }
    delete [] Temp;
    cerr<<"size of sampled subgroup vector: "<<NPartS<<endl;
    
    try {
      VecSubGroup[sub]=new subgroup_rec [NPartS];
    } catch(bad_alloc xa) {
      cerr<<"allocation failure for sampled groups in subcube: "<<iSC<<" and subsubcube: "
	  <<sub<<endl;
      exit(1);
    }
    
    for(long long j=0; j<NPartS; j++)
      VecSubGroup[sub][j]=TempS[j];
      
    delete [] TempS;
    
    SampleSize[sub] = NPartS; //write down size
    FLAG_SAMPLED[sub]=true;
  }
  
  else {
    cerr<<"SPECIES can only be PARTICLE, GALAXY, GROUP, SUBGROUP"<<endl;
    return;
  }
  
  FLAG_SAMPLED[sub]=true;
}


void load_SScube::read_species() 
  
{
  
  for(int sub=0; sub<NFILE0; sub++) 
    if(!FLAG_SAMPLED[sub])
      read_species(sub);
}


void load_SScube::format_data()

{

  for(int sub=0; sub<NFILE0; sub++) {
    if(!FLAG_SAMPLED[sub])
      read_species(sub);
    if(SPECIES=="PARTICLE") {
      NPTOT_SAMPLED=0;
      for(int s=0; s<NFILE; s++)
	NPTOT_SAMPLED+=NP_SAMPLED[s];
      for(long long j=0; j<SampleSize[sub]; j++) {
	vector<float> temp(3);
	for(int k=0; k<3; k++) {
	  if( VecPart[sub][j].pos[k] < SSH[sub].Min[k] || 
	      VecPart[sub][j].pos[k] > SSH[sub].Max[k] ) {
	    cerr<<"anomaly in read particles for subcube "<<iSC
		<<" sub-subcube "<<sub<<' '; 
	    output(VecPart[sub][j]); 
	    cerr<<endl;
	    exit(10);
	  }
	  temp[k]=VecPart[sub][j].pos[k];
	}
	particle Ob(temp);
	FINAL.push_back(Ob);
      }
    }

    else if(SPECIES=="GALAXY") {
      if(FLAG_BINNED) {
        for(int ib=0; ib<NBINS; ib++)
          for(long long j=0; j<SampleSizeB[sub][ib]; j++) {
            vector<float> temp(3);
            for(int k=0; k<3; k++) {
              if( VecGalB[sub][ib][j].Pos[k] < SSH[sub].Min[k] ) {
                  //VecGalB[sub][ib][j].Pos[k] > SSH[sub].Max[k] ) {
                cerr<<"anomaly in read galaxies for subcube "<<iSC
                    <<" sub-subcube "<<sub<<' ';
                output(VecGalB[sub][ib][j]);
		cerr<<endl;
                exit(10);
              }
              temp[k]=VecGalB[sub][ib][j].Pos[k];
            }
            particle Ob(temp);
            FINALB[ib].push_back(Ob);
	  }
      }

      else {
        for(long j=0; j<SampleSize[sub]; j++) {
          vector<float> temp(3);
          for(int k=0; k<3; k++) {
            if( VecGal[sub][j].Pos[k] < SSH[sub].Min[k] ||
                VecGal[sub][j].Pos[k] > SSH[sub].Max[k] ) {
              cerr<<"anomaly in read galaxies for subcube "<<iSC
                  <<" sub-subcube "<<sub<<' ';
              output(VecGal[sub][j]);
              cerr<<endl;
              exit(10);
            }
            temp[k]=VecGal[sub][j].Pos[k];
          }
          particle Ob(temp);
	  FINAL.push_back(Ob);
        }
      }
    }

    else if(SPECIES=="GROUP") {
      if(FLAG_BINNED) {
	for(int ib=0; ib<NBINS; ib++)
	  for(long long j=0; j<SampleSizeB[sub][ib]; j++) {
	    vector<float> temp(3);
	    for(int k=0; k<3; k++) {
	      if( VecGroupB[sub][ib][j].Pos[k] < SSH[sub].Min[k] || 
		  VecGroupB[sub][ib][j].Pos[k] > SSH[sub].Max[k] ) {
		cerr<<"anomaly in read groups for subcube "<<iSC
		    <<" sub-subcube "<<sub<<' '; 
		output(VecGroupB[sub][ib][j]); 
		cerr<<endl;
		exit(10);
	      }
	      temp[k]=VecGroupB[sub][ib][j].Pos[k];
	    }
	    particle Ob(temp);
	    FINALB[ib].push_back(Ob);
	  }
      }

      else {
	for(long j=0; j<SampleSize[sub]; j++) {
	  vector<float> temp(3);
	  for(int k=0; k<3; k++) {
	    if( VecGroup[sub][j].Pos[k] < SSH[sub].Min[k] || 
		VecGroup[sub][j].Pos[k] > SSH[sub].Max[k] ) {
	      cerr<<"anomaly in read groups for subcube "<<iSC
		  <<" sub-subcube "<<sub<<' '; 
	      output(VecGroup[sub][j]); 
	      cerr<<endl;
	      exit(10);
	    }
	    temp[k]=VecGroup[sub][j].Pos[k];
	  }
	  particle Ob(temp);
	  FINAL.push_back(Ob);
	}
      }
    }
    else if(SPECIES=="SUBGROUP") {
      for(long long j=0; j<SampleSize[sub]; j++) {
	vector<float> temp(3);
	for(int k=0; k<3; k++) {
	  if( VecSubGroup[sub][j].Pos[k] < SSH[sub].Min[k] || 
	      VecSubGroup[sub][j].Pos[k] > SSH[sub].Max[k] ) {
	    cerr<<"anomaly in read groups for subcube "<<iSC
		<<" sub-subcube "<<sub<<' '; 
	    output(VecSubGroup[sub][j]); 
	    cerr<<endl;
	    exit(10);
	  }
	  temp[k]=VecSubGroup[sub][j].Pos[k];
	}
	particle Ob(temp);
	FINAL.push_back(Ob);
      }
    }
    
    deallocate(sub);
  }

}



vector<particle> load_SScube::return_formatted_data()
  
{
  
  if(FINAL.size()==0)
    format_data();
  return FINAL;

}



vector<particle> load_SScube::return_formatted_data(int ibin)
  
{
  if(!FLAG_BINNED) {
    cerr<<"wrong call to function return_formatted_data(int ibin). The data is not binned."<<endl;
    exit(1);
  }
  if(ibin<0 || ibin>NBINS) {
    cerr<<"mass bin number out-of-bounds in 'return_formatted_data(ibin) "<<ibin<<endl;
    exit(1);
  }
  
  if(FINALB[ibin].size()==0) //what if the bin is naturally empty because of poor limits?
    format_data();
  
  return FINALB[ibin];
  
}



void load_SScube::deallocate(int sub)
  
{
  
  if(FLAG_SAMPLED[sub]) {
    if(SPECIES=="PARTICLE")
      delete [] VecPart[sub];
    else if(SPECIES=="GALAXY") {
      if(!FLAG_BINNED)
	delete [] VecGal[sub];
      else {
        for(int ib=0; ib<NBINS; ib++)
          delete [] VecGalB[sub][ib];
      }
    }
    else if(SPECIES=="GROUP") {
      if(!FLAG_BINNED)
	delete [] VecGroup[sub];
      else {
	for(int ib=0; ib<NBINS; ib++)
	  delete [] VecGroupB[sub][ib];
      }
    }
    else if(SPECIES=="SUBGROUP")
      delete [] VecSubGroup[sub];
    FLAG_SAMPLED[sub] = false; 
  }
  cerr<<"deallocating vector objects of type "<<SPECIES<<" for sub-subcube "
      <<sub<<" of subcube "<<iSC<<endl; 
  
}



load_SScube::~load_SScube()
  
{

  for(int sub=0; sub<NFILE0; sub++) {
    if(FLAG_SAMPLED[sub]) 
      deallocate(sub);
  }
  
  if( FINAL.size() !=0 )
    FINAL.clear();
  
  cerr<<"DESTROYING load_SScube object for subcube "<<iSC<<endl;
  
}



void load_SScube::output(group_rec &T)
  
{

  cerr<<"outputting group "<<T.Nr<<endl;
  cerr<<"length: "<<T.Len<<"; offset: "<<T.Offset<<"; Nr: "<<T.Nr<<endl;
  cerr<<"CM position: "<<T.CM[0]<<' '<<T.CM[1]<<' '<<T.CM[2]<<endl;
  cerr<<"velocity: "<<T.Vel[0]<<' '<<T.Vel[1]<<' '<<T.Vel[2]<<endl;
  cerr<<"potential center: "<<T.Pos[0]<<' '<<T.Pos[1]<<' '<<T.Pos[2]<<endl;
  cerr<<"M_mean200: "<<T.M_Mean200<<"; M_crit200: "<<T.M_Crit200<<"; M_tophat200: "
      <<T.M_TopHat200<<endl;
  cerr<<"Vel: "<<T.Vel_Disp<<"; Vel_mean: "<<T.Vel_Disp_Mean200<<"; Vel_crit: "
      <<T.Vel_Disp_Crit200<<"; Vel_TH: "<<T.Vel_Disp_TopHat200<<endl;
  cerr<<" Nsubs: "<<T.Nsubs<<"; first sub: "<<T.FirstSub<<endl;

  cerr<<"done outputting group "<<T.Nr<<endl<<endl;

}



void load_SScube::output(subgroup_rec &T)

{

  cerr<<"outputting subgroup "<<T.GrNr<<endl;
  cerr<<"length (nr particles): "<<T.Len<<"; GrNr: "<<T.GrNr<<"; SubNr: "<<T.SubNr<<endl;
  cerr<<"potential center: "<<T.Pos[0]<<' '<<T.Pos[1]<<' '<<T.Pos[2]<<endl;
  cerr<<"velocity: "<<T.Vel[0]<<' '<<T.Vel[1]<<' '<<T.Vel[2]<<endl;
  cerr<<"CM position: "<<T.CM[0]<<' '<<T.CM[1]<<' '<<T.CM[2]<<endl;
  cerr<<"spin: "<<T.Spin[0]<<' '<<T.Spin[1]<<' '<<T.Spin[2]<<endl;
  cerr<<"Vel: "<<T.Vel_Disp<<"; Vmax "<<T.Vmax<<"; VmaxRad: "<<T.VmaxRad
      <<"; HalfmassRad: "<<T.HalfmassRad<<endl;

  cerr<<"done outputting subgroup "<<T.GrNr<<endl<<endl;

}


void load_SScube::output(galaxy_rec &T)

{

  cerr<<"Galaxy type: "<<T.Type<<endl;
  cerr<<"Original snapshot: "<<T.Ghost<<endl;
  cerr<<"Mass of host halo: "<<T.CentralMvir<<endl;
  cerr<<"Galaxy position: "<<T.Pos[0]<<' '<<T.Pos[1]<<' '<<T.Pos[2]<<endl;
  cerr<<"Galaxy velocity: "<<T.Vel[0]<<' '<<T.Vel[1]<<' '<<T.Vel[2]<<endl;
  cerr<<"Mass of subhalo: "<<T.Mvir<<endl;
  cerr<<"Virial radius of subhalo: "<<T.Rvir<<endl;
  cerr<<"Virial velocity of subhalo: "<<T.Vvir<<endl;
  cerr<<"Distance to central galaxy: "<<T.DistanceToCentralGal[0]<<' '<<T.DistanceToCentralGal[1]
      <<' '<<T.DistanceToCentralGal[2]<<endl;
  cerr<<"Cold Gas: "<<T.ColdGas<<endl;
  cerr<<"Bulge mass: "<<T.BulgeMass<<endl;
  cerr<<"Disk mass: "<<T.DiskMass<<endl;
  cerr<<"Black Hole Mass: "<<T.BlackHoleMass<<endl;
  cerr<<"Rest-frame magnitudes: "<<T.ObsMagDust[0]<<' '<<T.ObsMagDust[1]<<' '<<T.ObsMagDust[2]
      <<' '<<T.ObsMagDust[3]<<' '<<T.ObsMagDust[4]<<endl;
  cerr<<endl;

}


void load_SScube::output(particle_data &T)

{

  cerr<<" position: "<<T.pos[0]<<' '<<T.pos[1]<<' '<<T.pos[2]<<endl;
  
}


void load_SScube::output(SScube_header &T)

{

  cerr<<"outputting sub-subcube header file for type: "<<SPECIES<<"s; subcube: "<<iSC<<endl;

  cerr<<"redshift: "<<T.redshift<<"; time: "<<T.time<<"; BoxSize: "<<T.BoxSize
      <<"; subcube size: "<<T.LSubCube<<endl;

  cerr<<"total number of subcubes in sim: "<<T.NSubCubeTOT
      <<"; number of sub-subcubes per subcube: "<<T.NumberOfFilesPerSubCube<<endl;

  cerr<<"total number of "<<SPECIES<<"s in the sim.:"<<T.NSubCubePartTOT
      <<"  ; and in the sub-subcube: "<<T.NSubCubePart<<endl;

  cerr<<"mass: "<<T.mass<<"; min "<<SPECIES<<" position: "<<T.Min[0]<<' '<<T.Min[1]<<' '<<T.Min[2]
      <<" ; and max "<<SPECIES<<" position: "<<T.Max[0]<<' '<<T.Max[1]<<' '<<T.Max[2]<<endl;

  cerr<<"cosmology: Om: "<<T.Omega0<<"; OmLambda: "<<T.OmegaLambda
      <<"; Hubble: "<<T.HubbleParam<<endl;

  cerr<<"done outputting header file for type: "<<SPECIES<<"s; subcube: "<<iSC<<endl;

}



//read the galaxy data


read_galaxy::read_galaxy(const int snap_, const string path_in_):snapshot(snap_), path_in(path_in_)

{

  ostringstream oss;
  oss<<snapshot;
  if(snapshot==0)
    oss<<0; 
  file_gen = path_in + "SA_z0." + oss.str() + '_';
  for(int i=0; i<NFILE; i++)
    FLAG_READ[i]=false;
}


void read_galaxy::open_subfile(int sub)

{

  if(sub>=NFILE) {
    cerr<<"wrong file number for Raul's galaxies. Maximum number is: "<<NFILE-1<<endl;
    exit(1);
  }
  ostringstream oss;
  int Ntrees, Ngals;
  oss<<sub;
  string file_in=file_gen + oss.str();

  ifstream in(file_in.c_str(), ios::in);
  check_stream(in, file_in);
  in.read((char*) &Ntrees, sizeof(int)); 
  in.read((char*) &Ngals, sizeof(int));
  cerr<<"number of trees: "<<Ntrees<<" and number of galaxies: "<<Ngals<<endl<<endl;

  int O[Ntrees];
  try {
    GalVec[sub] = new galaxy_rec [Ngals]; 
  } catch(bad_alloc xa) {
    cerr<<"allocation failure for galaxies in file : "<<sub<<endl;
    exit(1);
  }

  for(int i=0; i<Ntrees; i++)
    in.read((char*) &O[i], sizeof(int));
  //or in.read((char*) &O[0], Ntrees*sizeof(int));
  for(int i=0; i<Ngals; i++) 
    in.read((char*) &GalVec[sub][i], sizeof(galaxy_rec));    
  //or in.read((char*) &V[0], Ngals*sizeof(galaxy_rec));
  in.close();

  //output(GalVec[sub][0]); cerr<<endl;
  //output(GalVec[sub][1]); cerr<<endl;
  //output(GalVec[sub][Ngals-1]); cerr<<endl;

  Ngalaxies[sub]=Ngals;
  FLAG_READ[sub]=true;
  
}


void read_galaxy::deallocate(int sub)

{
  if(sub>=NFILE) {
    cerr<<"wrong file number for Raul's galaxies. Maximum number is: "<<NFILE-1<<endl;
    exit(1);
  }
  
  if(FLAG_READ[sub]) {
    delete [] GalVec[sub];
    FLAG_READ[sub]=false;
    cerr<<"deallocating galaxy vector from file "<<sub<<endl;
  }

}

void read_galaxy::output(galaxy_rec &T)

{

  cerr<<"Galaxy type: "<<T.Type<<endl;
  cerr<<"Original snapshot: "<<T.Ghost<<endl;
  cerr<<"Mass of host halo: "<<T.CentralMvir<<endl;
  cerr<<"Galaxy position: "<<T.Pos[0]<<' '<<T.Pos[1]<<' '<<T.Pos[2]<<endl;
  cerr<<"Galaxy velocity: "<<T.Vel[0]<<' '<<T.Vel[1]<<' '<<T.Vel[2]<<endl;
  cerr<<"Mass of subhalo: "<<T.Mvir<<endl;
  cerr<<"Virial radius of subhalo: "<<T.Rvir<<endl;
  cerr<<"Virial velocity of subhalo: "<<T.Vvir<<endl;
  cerr<<"Distance to central galaxy: "<<T.DistanceToCentralGal[0]<<' '<<T.DistanceToCentralGal[1]
      <<' '<<T.DistanceToCentralGal[2]<<endl;
  cerr<<"Cold Gas: "<<T.ColdGas<<endl;
  cerr<<"Bulge mass: "<<T.BulgeMass<<endl;
  cerr<<"Disk mass: "<<T.DiskMass<<endl;
  cerr<<"Black Hole Mass: "<<T.BlackHoleMass<<endl;
  cerr<<"Rest-frame magnitudes: "<<T.ObsMagDust[0]<<' '<<T.ObsMagDust[1]<<' '<<T.ObsMagDust[2]
      <<' '<<T.ObsMagDust[3]<<' '<<T.ObsMagDust[4]<<endl;
  cerr<<endl;

}


read_galaxy_Millennium::read_galaxy_Millennium(const double red_, const string path_in_):
  red(red_), path_in(path_in_)
{

  ostringstream oss;
  oss<<red;
  file_gen = path_in + "SA_z" + oss.str() + '_';
  for(int i=0; i<NFILE; i++)
    FLAG_READ[i]=false;
}


void read_galaxy_Millennium::open_subfile(int sub)

{

  if(sub>=NFILE) {
    cerr<<"wrong file number for Bruno's galaxies. Maximum number is: "<<NFILE-1<<endl;
    exit(1);
  }
  ostringstream oss;
  oss<<sub;
  string file_in=file_gen + oss.str();
  int Ntrees, Ngals;  
  ifstream in(file_in.c_str(), ios::in);
  check_stream(in, file_in);
  in.read((char*) &Ntrees, sizeof(int));
  in.read((char*) &Ngals, sizeof(int));
  cerr<<"number of trees: "<<Ntrees<<" and number of galaxies: "<<Ngals<<endl<<endl;

  int O[Ntrees];
  try {
    GalVec[sub] = new galaxy_rec_Millennium [Ngals];
  } catch(bad_alloc xa) {
    cerr<<"allocation failure for galaxies in file : "<<sub<<endl;
    exit(1);
  }

  for(int i=0; i<Ntrees; i++)
    in.read((char*) &O[i], sizeof(int));
  //or in.read((char*) &O[0], Ntrees*sizeof(int));                                              
  for(int i=0; i<Ngals; i++)
    in.read((char*) &GalVec[sub][i], sizeof(galaxy_rec_Millennium));
  //or in.read((char*) &V[0], Ngals*sizeof(galaxy_rec));                        
  in.close();
  //output(GalVec[sub][0]); cerr<<endl;                                                            
  //output(GalVec[sub][1]); cerr<<endl;                                                            
  //output(GalVec[sub][Ngals-1]); cerr<<endl;                                                    
  Ngalaxies[sub]=Ngals;
  FLAG_READ[sub]=true;

}

void read_galaxy_Millennium::deallocate(int sub)

{
  if(sub>=NFILE) {
    cerr<<"wrong file number for Bruno's galaxies. Maximum number is: "<<NFILE-1<<endl;
    exit(1);
  }

  if(FLAG_READ[sub]) {
    delete [] GalVec[sub];
    FLAG_READ[sub]=false;
    cerr<<"deallocating galaxy vector from file "<<sub<<endl;
  }
}

void read_galaxy_Millennium::output(galaxy_rec_Millennium &T)

{

  cerr<<"Galaxy type: "<<T.Type<<endl;
  cerr<<"Mass of host halo: "<<T.CentralMvir<<endl;
  cerr<<"Galaxy position: "<<T.Pos[0]<<' '<<T.Pos[1]<<' '<<T.Pos[2]<<endl;
  cerr<<"Galaxy velocity: "<<T.Vel[0]<<' '<<T.Vel[1]<<' '<<T.Vel[2]<<endl;
  cerr<<"Mass of subhalo: "<<T.Mvir<<endl;
  cerr<<"Virial radius of subhalo: "<<T.Rvir<<endl;
  cerr<<"Virial velocity of subhalo: "<<T.Vvir<<endl;
  cerr<<"Distance to central galaxy: "<<T.DistanceToCentralGal[0]<<' '<<T.DistanceToCentralGal[1]
      <<' '<<T.DistanceToCentralGal[2]<<endl;
  cerr<<"Cold Gas: "<<T.ColdGas<<endl;
  cerr<<"Bulge mass: "<<T.BulgeMass<<endl;
  cerr<<"Disk mass: "<<T.DiskMass<<endl;
  cerr<<"Black Hole Mass: "<<T.BlackHoleMass<<endl;
  cerr<<"SDSS Dust magnitudes: "<<T.MagDust[15]<<' '<<T.MagDust[16]<<' '<<T.MagDust[17]
      <<' '<<T.MagDust[18]<<' '<<T.MagDust[19]<<endl;
  cerr<<"SDSS magnitudes: "<<T.Mag[15]<<' '<<T.Mag[16]<<' '<<T.Mag[17]
      <<' '<<T.Mag[18]<<' '<<T.Mag[19]<<endl;
  cerr<<"SDSS Bulge magnitudes: "<<T.MagBulge[15]<<' '<<T.MagBulge[16]<<' '<<T.MagBulge[17]
      <<' '<<T.MagBulge[18]<<' '<<T.MagBulge[19]<<endl;
  cerr<<endl;

}
