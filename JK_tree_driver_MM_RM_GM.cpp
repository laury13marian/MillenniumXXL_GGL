#include "mpi.h"
#include "read_data.h"
#include "tree_processing_projected.h"

using namespace std;

int main(int argc, char *argv[])

{
  
  if(argc != 4) {
    cerr<<"must input subcube number, MM sampling rate, JKsample."<<endl;
    exit(10);
  }

  MPI::Init(argc, argv);
  int ntasks, id, master, tag;
  double wtime, InitTime, FinTime;
  id=MPI::COMM_WORLD.Get_rank(); //get this processor's id 
  ntasks=MPI::COMM_WORLD.Get_size(); //get number of processors 
  wtime=MPI::Wtime();
  InitTime=MPI::Wtime();
  MPI::Status status;
  master=0;

  const int snap = 54;
  const size_t NRAN = 1000000;
  const float zMAX = 100.;
  const double rbox=500.;
  const double DeltaJK=62.5; //JKno=8 for this value
  int JKno=static_cast<int>(rbox/DeltaJK); //number of JK samples per side
  cerr<<"number of JK samples: "<<JKno*JKno<<endl;
  const int ns = atoi(argv[1])-1; //because the script starts subcube at 1, not 0!
  const int nsampM=atoi(argv[2]);
  const int JKsample = atoi(argv[3])-1;  //the JK sample for which the calculation is done
  int Jix, Jiy;  
  //JKsample=Jix*NJK+Jiy
  Jiy=JKsample % JKno;
  Jix=(JKsample-Jiy)/JKno;
  
  if(id==master) {
    cerr<<"Run specifics: snap: "<<snap<<" subcube: "<<ns<<" nsampM: "
	<<nsampM<<" NRAN: "<<NRAN<<" NTASKS: "<<ntasks<<" zmax: "<<zMAX<<endl;
    cerr<<"CALCULATIONS FOR JK SAMPLE: "<<JKsample<<' '<<Jix<<' '<<Jiy<<endl<<endl;	
  }

  const float rMAX=100.;//110.;//30.;                                              
  const float rMIN=0.01;
  const float zMIN=0.;
  const int NR=30; //20; //10;
  const int NZ=10;
  float rb[3][NR]; //radial bin vector
  float zb[3][NZ]; //z-bin vector      
  const int dim=3;

  const int Nmags=5;
  const int NM=Nmags-1; //number of magnitude/mass bins   
  vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; 
  const int iMAG=1;
  
  double npairsMM[ntasks][NR][NZ], npairsGM[ntasks][NM][NR][NZ], npairsRM[ntasks][NR][NZ];
  string file_out_MM_RM[ntasks], file_out_GM[ntasks];

  ostringstream osNR, osNZ, osNM, osSNAP, osSC, osNTASK, osNRAN, osZMAX, osSAMPM, osJK;
  osNR<<NR; osNZ<<NZ; osNM<<NM; osSNAP<<0;  osSNAP<<snap; osJK<<JKsample;
  osSC<<ns; osNTASK<<ntasks; osNRAN<<NRAN; osZMAX<<zMAX; osSAMPM<<nsampM;  

  const string path_out = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/JACK_KNIFE/";
  //to register the MM and RM results
  const string file_MM_RM = path_out + "PAIR_COUNTS_MM_RM/PC_snap_" + osSNAP.str() + "_sub_" 
    + osSC.str() + "_JKsample_" + osJK.str() + "_binsR_" + osNR.str() + "_NRAN_" + osNRAN.str() 
    + "_nsamp_" + osSAMPM.str()+  "_NT_" + osNTASK.str() + "_task_" ;

  const string file_GM = path_out + "PAIR_COUNTS_GM/PC_snap_" + osSNAP.str() + "_sub_"
    + osSC.str() + "_JKsample_" + osJK.str() + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + 
    + "_nsamp_" + osSAMPM.str() + "_NT_" + osNTASK.str() + "_task_" ;

  if(id==0) 
    cerr<<"beginning the master process for snapshot: "<<snap<<" subcube: "<<ns
	<<" JK sample: "<<Jix<<' '<<Jiy<<' '<<JKsample<<" MM sampling rate 1/"<<nsampM
	<<" NRAN: "<<NRAN<<endl;
  
  tag=id;
  cerr<<"beginning slave process :"<<id<<endl;
  double wtime1, wtime2;
  vector<string> path_data(4);
  path_data[0] = 
    "/ptmp/mpa/rosm/Millennium/Millennium-XXL/subcubes/subcubes_" + osSNAP.str() +'/';
  path_data[1] = "/snapshot_" + osSNAP.str() + "_SubCube_"; 
  
  const string SPM = "PARTICLE"; 
  const string SPG = "GALAXY"; //"GROUP";
  const signed long rng_seed_gen = -1381718569;
  //const signed long rng_seed_gen = -static_cast<signed long>(std::time(NULL));
  const signed long rng_seed = rng_seed_gen + id*id*id*id;
  const signed long rng_seed_GAL = rng_seed_gen/4 - 123*id*id*id*id;
  cerr<<"random seed used for process: "<<id<<' '<<rng_seed_gen<<' '<<rng_seed_GAL<<' '<<rng_seed<<endl; 
  

  cerr<<"READING PARTICLE DATA"<<endl;
  wtime1=MPI::Wtime();
  load_SScube OBM(ns, path_data, SPM, nsampM, rng_seed);
  vector<double> Bound=OBM.lower_bound();
  if(id==master) {
    cerr<<"lower subcube bound: ";
    for(int i=0; i<Bound.size(); i++)
      cerr<<Bound[i]<<' ';
    cerr<<endl<<endl;
  }

  vector<particle> VpartM0, VpartM;
  VpartM0=OBM.return_formatted_data();

  cerr<<"INITIAL size of formatted particles: "<<VpartM0.size()<<endl;
  for(int j=0; j<VpartM0.size(); j++) {
    if( (VpartM0[j].pos[0] >= Bound[0] + Jix*DeltaJK) && (VpartM0[j].pos[0] < Bound[0] + (Jix+1)*DeltaJK)
	&& (VpartM0[j].pos[1] >= Bound[1] + Jiy*DeltaJK) && (VpartM0[j].pos[1] < Bound[1] + (Jiy+1)*DeltaJK) )
      continue; 
    else
      VpartM.push_back(VpartM0[j]);
    }
  cerr<<"FINAL size for matter particle vector: "<<VpartM.size()<<endl;
  //now delete first vector
  vector<particle>::iterator p;
  p=VpartM0.begin();
  VpartM0.erase(p, p+VpartM0.size());
  cerr<<"erased vector new size: "<<VpartM0.size()<<endl<<endl;
 
  wtime2=MPI::Wtime();
  cerr<<"time to read and format the particle data for slave process: "
      <<id<<" is: "<<wtime2-wtime1<<endl;
  
  
  cerr<<"READING GALAXY DATA"<<endl;
  wtime1=MPI::Wtime();
  load_SScube OBG(ns, path_data, SPG, rng_seed_GAL, mags, iMAG);
  vector<particle> VpartG0[NM], VpartG[NM];
  for(int im=0; im<NM; im++)
    VpartG0[im]=OBG.return_formatted_data(im);

  for(int im=0; im<NM; im++) {
    cerr<<"initial size of galaxy vector for im="<<im<<" is: "<<VpartG0[im].size()<<endl;
    for(int j=0; j<VpartG0[im].size(); j++) {
      if( (VpartG0[im][j].pos[0] >= Bound[0] + Jix*DeltaJK) && (VpartG0[im][j].pos[0] < Bound[0] + (Jix+1)*DeltaJK)
	  && (VpartG0[im][j].pos[1] >= Bound[1] + Jiy*DeltaJK) && (VpartG0[im][j].pos[1] < Bound[1] + (Jiy+1)*DeltaJK) ) 
	continue; 
      else
	VpartG[im].push_back(VpartG0[im][j]);
    }
    cerr<<"final size for galaxy vector for im="<<im<<" is: "<<VpartG[im].size()<<endl;
    //now delete first vector
    vector<particle>::iterator p;
    p=VpartG0[im].begin();
    VpartG0[im].erase(p, p+VpartG0[im].size());
    cerr<<"erased vector new size: "<<VpartG0[im].size()<<endl<<endl;
  }

  OBG.info_output();
  wtime2=MPI::Wtime();
  cerr<<"time to read and format the group data for slave process: "
      <<id<<" is: "<<wtime2-wtime1<<endl;

  const signed long rng_seed1 = rng_seed_gen + id*id;
  const signed long rng_seed2 = rng_seed_gen/3 + id*id*id*id;
  const signed long rng_seed3 = 2*rng_seed_gen + id*id;
  cerr<<"seeds for random catalogue for process "<<id<<' '
      <<rng_seed1<<' '<<rng_seed2<<' '<<rng_seed3<<endl;
  vector<particle> VpartR;
  
  random2 OBran1(rng_seed1), OBran2(rng_seed2), OBran3(rng_seed3);
  wtime1=MPI::Wtime();
  for(size_t i=0; i<NRAN; i++) {
    vector<float> temp(3);
    temp[0]=Bound[0] + rbox*OBran1(); temp[1]=Bound[1] + rbox*OBran2(); temp[2]=Bound[2] + rbox*OBran3();
    if( (temp[0] >= Bound[0] + Jix*DeltaJK ) && (temp[0] < Bound[0] + (Jix+1) * DeltaJK ) && 
	(temp[1] >= Bound[1] + Jiy*DeltaJK ) && (temp[1] < Bound[1] + (Jiy+1) * DeltaJK ) )
      continue;
    else {
      particle Otemp(temp);
      VpartR.push_back(Otemp);
    }
  } 
  cerr<<"size of random groups: "<<VpartR.size()<<endl;
  wtime2=MPI::Wtime();
  cerr<<"time to generate random data for slave process:"<<id<<" is: "<<wtime2-wtime1<<endl;

  //now build tree

  wtime1=MPI::Wtime();
  build_tree OtreeM(rbox, VpartM);
  wtime2=MPI::Wtime();
  cerr<<"time taken to build the M tree by process: "<<id<<" = "<<wtime2-wtime1<<endl;
  build_tree *OM1, *OM2;
  OM1=&OtreeM; OM2=&OtreeM; 
  
  wtime1=MPI::Wtime();
  build_tree OtreeR(rbox, VpartR);
  wtime2=MPI::Wtime();
  cerr<<"time taken to build the R tree by process: "<<id<<" = "<<wtime2-wtime1<<endl;
  build_tree *OR;
  OR=&OtreeR; 
  
  const vector<float> vecPER(dim, 0.);

  projected_correlations::BinType bt = projected_correlations::LOG10;
  bool flag_per=false; //true;
  cerr<<"box size: "<<rbox<<endl;
  
  projected_correlations OcorrMM(OM1, OM2, bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER);
  projected_correlations OcorrRM(OM1, OR, bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(id==master) {
    vector<float> bin_R=OcorrMM.return_bins_rperp();
    vector<float> bin_Z=OcorrMM.return_bins_z();
    for(int i=0; i<NR; i++) {
      rb[0][i]=bin_R[i]; 
      rb[1][i]=bin_R[2*NR+i];
      rb[2][i]=bin_R[NR+i];
    }

    for(int i=0; i<NZ; i++) {
      zb[0][i]=bin_Z[i]; 
	zb[1][i]=bin_Z[2*NZ+i];
	zb[2][i]=bin_Z[NZ+i];
    }
  }
  
  const string dummy = "dummy";
  
  wtime1=MPI::Wtime();
  vector<double> tcorrMM=OcorrMM.dual_tree_2p(dummy);
  for(int j=0; j<NR; j++)
    for(int jz=0; jz<NZ; jz++)
      npairsMM[id][j][jz]=tcorrMM[j*NZ+jz];
    wtime2=MPI::Wtime();
    MPI::COMM_WORLD.Barrier();
    if(id==master) 
      cerr<<"time taken to compute MM correlations:= "<<wtime2-wtime1<<endl<<endl;
    
    wtime1=MPI::Wtime();
    vector<double> tcorrRM=OcorrRM.dual_tree_2p(dummy);
    for(int j=0; j<NR; j++)
      for(int jz=0; jz<NZ; jz++)
        npairsRM[id][j][jz]=tcorrRM[j*NZ+jz];
    wtime2=MPI::Wtime();
    MPI::COMM_WORLD.Barrier();
    if(id==master)
      cerr<<"time taken to compute RM correlations:= "<<wtime2-wtime1<<endl<<endl;

    if(id!=master)
      {
	MPI::COMM_WORLD.Send(&npairsMM[id], NR*NZ, MPI::DOUBLE, master, tag);
	MPI::COMM_WORLD.Send(&npairsRM[id], NR*NZ, MPI::DOUBLE, master, tag);
      }

    else {
      //now the master part
      for(int np=1; np<ntasks; np++)
	{
	  tag=np;
	  MPI::COMM_WORLD.Recv(&npairsMM[np], NR*NZ, MPI::DOUBLE, np, tag, status);
	  MPI::COMM_WORLD.Recv(&npairsRM[np], NR*NZ, MPI::DOUBLE, np, tag, status);
	}
      
      for(int np=0; np<ntasks; np++) {
	ostringstream oss;
	oss<<np;
        file_out_MM_RM[np] = file_MM_RM + oss.str() + ".dat";
	ofstream out(file_out_MM_RM[np].c_str(), ios::out);
	for(int jr=0; jr<NR; jr++)
	  for(int iz=0; iz<NZ; iz++)
	    out<<rb[1][jr]<<' '<<zb[1][iz]<<' '<<npairsMM[np][jr][iz]<<' '<<npairsRM[np][jr][iz]<<endl;
	out.close();
      }
      cerr<<"the master has written the MM and RM output files: " + file_MM_RM<<endl;
    }

    //compute the GM correlations

    wtime1=MPI::Wtime();
    build_tree OtreeG[NM] = { build_tree(rbox, VpartG[0]), build_tree(rbox, VpartG[1]),
			      build_tree(rbox, VpartG[2]), build_tree(rbox, VpartG[3])} ;
    wtime2=MPI::Wtime();
    cerr<<"time taken to build the G tree by process: "<<id<<" = "<<wtime2-wtime1<<endl;
    build_tree *OG[NM];
    for(int ib=0; ib<NM; ib++)
      OG[ib]=&OtreeG[ib];
    
    
    projected_correlations OcorrGM[NM] = { projected_correlations(OM1, OG[0], bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER),
					   projected_correlations(OM1, OG[1], bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER),
					   projected_correlations(OM1, OG[2], bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER),
					   projected_correlations(OM1, OG[3], bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER) };

    vector<double> tcorrGM[NM];
    const int NMPAR=NM;
    for(int im=0; im<NMPAR; im++) {
      wtime1=MPI::Wtime();
      tcorrGM[im]=OcorrGM[im].dual_tree_2p(dummy);
      for(int jr=0; jr<NR; jr++)
        for(int jz=0; jz<NZ; jz++)
          npairsGM[id][im][jr][jz]=tcorrGM[im][jr*NZ+jz];
      wtime2=MPI::Wtime();
      if(id==master)
        cerr<<"time taken to compute GM correlations for mag bin: "<<im<<" is: = "<<wtime2-wtime1<<endl<<endl;
    }
    
    MPI::COMM_WORLD.Barrier();
    
    if(id!=master) 
      {
	for(int im=0; im<NMPAR; im++) 
	  MPI::COMM_WORLD.Send(&npairsGM[id][im], NR*NZ, MPI::DOUBLE, master, tag);
      }
    
    else {      
      //now the master part
      for(int np=1; np<ntasks; np++) 
	{
	  tag=np;
	  for(int im=0; im<NMPAR; im++) 
	    MPI::COMM_WORLD.Recv(&npairsGM[np][im], NR*NZ, MPI::DOUBLE, np, tag, status);
	}
      //output files
      
      for(int np=0; np<ntasks; np++) {
        ostringstream oss;
        oss<<np;
	file_out_GM[np] = file_GM + oss.str() + ".dat";
	ofstream outGM(file_out_GM[np].c_str(), ios::out);
        for(int jr=0; jr<NR; jr++)
          for(int iz=0; iz<NZ; iz++) {
            outGM<<rb[1][jr]<<' '<<zb[1][iz];
	    for(int im=0; im<NMPAR; im++)
	      outGM<<' '<<mags[im]<<' '<<mags[im+1]<<' '<<npairsGM[np][im][jr][iz];
	    outGM<<endl;
	  }
	outGM.close();
	cerr<<"done writing file "<<file_out_GM[np]<<endl;      
      }
    }
    
    MPI::COMM_WORLD.Barrier();
    FinTime=MPI::Wtime();
    cerr<<"total time taken: "<<FinTime-InitTime<<endl<<endl;
    MPI::Finalize();
    
    return 0;
    
}

 
