//this code will only do RG and GG. Later maybe G1G2 (i.e. cross correlation of different bins)
//the GM, MM and RM are done separately.

#include "mpi.h"
#include "read_data.h"
#include "tree_processing_projected.h"

using namespace std;

int main(int argc, char *argv[])

{
  
  if(argc != 3) {
    cerr<<"must input subcube number, JKsample number"<<endl; //varying from 1 to (Lsub/DeltaJK)^2
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
  const int JKsample = atoi(argv[2])-1;  //the JK sample for which the calculation is done
  int Jix, Jiy;
  //JKsample=Jix*NJK+Jiy
  Jiy=JKsample % JKno;
  Jix=(JKsample-Jiy)/JKno;
  
  if(id==master) {
    cerr<<"Run specifics: snap: "<<snap<<" subcube: "<<ns<<" NRAN: "<<NRAN
	<<" NTASKS: "<<ntasks<<" zmax: "<<zMAX<<endl<<endl;
    cerr<<"CALCULATIONS FOR JK SAMPLE: "<<JKsample<<' '<<Jix<<' '<<Jiy<<endl<<endl;
  }

  const float rMAX=100.;//110.;//30.;                                              
  const float rMIN=0.01;
  const float zMIN=0.;
  int NR=30; //temporarily
  const int NZ=10;
  float rb[3][NR]; //radial bin vector
  float zb[3][NZ]; //z-bin vector      
  const int dim=3;
  //define magnitude bins, from the brightest (most negative) to the faintest (least negative)

  const int Nmags=5;
  const int NM=Nmags-1; //number of magnitude/mass bins
  vector<double> mags(Nmags, 0.);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; 
  const int iMAG=1; //red band

  double npairsGG[ntasks][NM][NR][NZ], npairsRG[ntasks][NM][NR][NZ]; //array containing the pair counts
  string file_out_GG[ntasks], file_out_RG1[ntasks];

  ostringstream osNR, osNZ, osNM, osSNAP, osSC, osNTASK, osNRAN, osZMAX, osJK;
  osNR<<NR; osNZ<<NZ; osNM<<NM; osSNAP<<0;  osSNAP<<snap; osJK<<JKsample;
  osSC<<ns; osNTASK<<ntasks; osNRAN<<NRAN; osZMAX<<zMAX;
  const string path_out = 
    "/u/lmarian/LauraTreeCode/XXL_CALCULATIONS/RESULTS/PROJECTED_CORR_FUNC/NEW_BINS/JACK_KNIFE/";
  //to register the GG and RG results
  const string file_GG = path_out + "PAIR_COUNTS_GG_RG/GG_snap_" + osSNAP.str() + "_sub_" + osSC.str() + 
    "_JKsample_" + osJK.str() + "_binsM_" + osNM.str() + "_binsR_" + osNR.str()  + "_NT_" + osNTASK.str() 
    + "_task_" ;
  
  const string file_RG = path_out + "PAIR_COUNTS_GG_RG/RG_snap_" + osSNAP.str() + "_sub_" + osSC.str()  + 
    "_JKsample_" + osJK.str() + "_binsM_" + osNM.str() + "_binsR_" + osNR.str() + "_NRAN_" + osNRAN.str() 
    + "_NT_" + osNTASK.str() + "_task_" ;
  
  if(id==master) 
    cerr<<"beginning the GG_RG master process for snapshot: "<<snap<<" subcube: "<<ns<<" NRAN: "<<NRAN<<endl;

  tag=id;
  cerr<<"beginning slave process :"<<id<<endl;
  double wtime1, wtime2;
  vector<string> path_data(4);
  path_data[0] = 
    "/ptmp/mpa/rosm/Millennium/Millennium-XXL/subcubes/subcubes_" + osSNAP.str() +'/';
  path_data[1] = "/snapshot_" + osSNAP.str() + "_SubCube_"; 
  
  const string SPG = "GALAXY"; //"GROUP";
  const signed long rng_seed_gen = -1381718569;
  //const signed long rng_seed_gen = -static_cast<signed long>(std::time(NULL));
  const signed long rng_seed_GAL = rng_seed_gen/4 - 123*id*id*id*id;
  cerr<<"random seed used for process: "<<id<<' '<<rng_seed_gen<<" and for galaxies: "<<rng_seed_GAL<<endl; 
  
  cerr<<"READING GALAXY DATA"<<endl;
  wtime1=MPI::Wtime();
  load_SScube OBG(ns, path_data, SPG, rng_seed_GAL, mags, iMAG);
  vector<particle> VpartG0[NM], VpartG[NM]; 
  for(int im=0; im<NM; im++)
    VpartG0[im]=OBG.return_formatted_data(im);  
  vector<double> Bound=OBG.lower_bound();
  if(id==master) {
    cerr<<"lower subcube bound: ";
    for(int i=0; i<Bound.size(); i++)
      cerr<<Bound[i]<<' ';
    cerr<<endl<<endl;
  }

  //now exclude the paralleliped of interest, defined by Jix, Jiy

  for(int im=0; im<NM; im++) {
    cerr<<"old size of galaxy vector for im="<<im<<" is: "<<VpartG0[im].size()<<endl;
    for(int j=0; j<VpartG0[im].size(); j++) {
      if( (VpartG0[im][j].pos[0] >= Bound[0] + Jix*DeltaJK) && (VpartG0[im][j].pos[0] < Bound[0] + (Jix+1)*DeltaJK)
	  && (VpartG0[im][j].pos[1] >= Bound[1] + Jiy*DeltaJK) && (VpartG0[im][j].pos[1] < Bound[1] + (Jiy+1)*DeltaJK) ) {
	//cerr<<VpartG0[im][j].pos[0]<<' '<<VpartG0[im][j].pos[1]<<' '<<VpartG0[im][j].pos[2]<<endl;
	continue; }
      else
	VpartG[im].push_back(VpartG0[im][j]);
    }
    cerr<<"new size if galaxy vector for im="<<im<<" is: "<<VpartG[im].size()<<endl;
    //now delete first vector
    vector<particle>::iterator p;
    p=VpartG0[im].begin();
    VpartG0[im].erase(p, p+VpartG0[im].size()-1);
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
  
  //now build tree for binned galaxies
  
  wtime1=MPI::Wtime(); 
  build_tree OtreeG[NM] = { build_tree(rbox, VpartG[0]), build_tree(rbox, VpartG[1]),
			    build_tree(rbox, VpartG[2]), build_tree(rbox, VpartG[3])} ;
  wtime2=MPI::Wtime();
  cerr<<"time taken to build the G tree by process: "<<id<<" = "<<wtime2-wtime1<<endl;
  build_tree *OG1[NM], *OG2[NM];
  for(int ib=0; ib<NM; ib++) {
    OG1[ib]=&OtreeG[ib]; OG2[ib]=&OtreeG[ib];
  }
  
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
  //NR=NR+2; //CAREFUL HERE!!!!
  projected_correlations OcorrGG[NM] = { projected_correlations(OG1[0], OG2[0], bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER),
					 projected_correlations(OG1[1], OG2[1], bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER),
					 projected_correlations(OG1[2], OG2[2], bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER),
					 projected_correlations(OG1[3], OG2[3], bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER)};
  //NR=NR-2; //CAREFUL AGAIN
  if(id==master) {
    cerr<<"THE RADIAL BINS FROM THE MASTER: "<<endl;
    vector<float> bin_R=OcorrGG[0].return_bins_rperp();
    cerr<<"size of returned bins: "<<bin_R.size()<<endl;
    vector<float> bin_Z=OcorrGG[0].return_bins_z();
    for(int i=0; i<NR; i++) {
      rb[0][i]=bin_R[i];
      rb[1][i]=bin_R[2*NR+i];
      rb[2][i]=bin_R[NR+i];
      cerr<<rb[0][i]<<' '<<rb[1][i]<<' '<<rb[2][i]<<endl;
    }
    cerr<<endl;
    for(int i=0; i<NZ; i++) {
      zb[0][i]=bin_Z[i];
      zb[1][i]=bin_Z[2*NZ+i];
      zb[2][i]=bin_Z[NZ+i];
    }
  }

  const string dummy = "dummy";

  vector<double> tcorrGG[NM];
  for(int im=0; im<NM; im++) {
    wtime1=MPI::Wtime();
    tcorrGG[im]=OcorrGG[im].dual_tree_2p(dummy);
    //if(id==master)
      //cerr<<"computing corr for mag bin "<<im<<endl;
    for(int jr=0; jr<NR; jr++)
      for(int jz=0; jz<NZ; jz++) {
	npairsGG[id][im][jr][jz]=tcorrGG[im][jr*NZ+jz];
      }
    wtime2=MPI::Wtime();
    if(id==master)
      cerr<<"time taken to compute GG correlations for mag bin: "<<im<<" is: = "<<wtime2-wtime1<<endl<<endl;
  }
  MPI::COMM_WORLD.Barrier();
  
  if(id!=master)
    {
      for(int im=0; im<NM; im++)
	MPI::COMM_WORLD.Send(&npairsGG[id][im], NR*NZ, MPI::DOUBLE, master, tag);
    }

  else {
    //now the master part                                                                                                                    
    for(int np=1; np<ntasks; np++)
      {
	tag=np;
	for(int im=0; im<NM; im++)
	  MPI::COMM_WORLD.Recv(&npairsGG[np][im], NR*NZ, MPI::DOUBLE, np, tag, status);
      }
    
    for(int np=0; np<ntasks; np++) {
      ostringstream oss;
      oss<<np;
      file_out_GG[np] = file_GG + oss.str() + ".dat";
    }
    
    for(int np=0; np<ntasks; np++) {
      ofstream out(file_out_GG[np].c_str(), ios::out);
      for(int jr=0; jr<NR; jr++)
	for(int iz=0; iz<NZ; iz++) {
	  out<<rb[1][jr]<<' '<<zb[1][iz];
	  for(int im=0; im<NM; im++)
	    out<<' '<<mags[im]<<' '<<mags[im+1]<<' '<<npairsGG[np][im][jr][iz];                         
	  out<<endl;
	}
      out.close();
      
      cerr<<"done writing file "<<file_out_GG[np]<<endl;
    }
  }
  MPI::COMM_WORLD.Barrier();
  //NR=NR+2;  //CAREFUL HERE!!!!				 
  projected_correlations OcorrRG[NM] = { projected_correlations(OG1[0], OR, bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER),
                                         projected_correlations(OG1[1], OR, bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER),
                                         projected_correlations(OG1[2], OR, bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER),
                                         projected_correlations(OG1[3], OR, bt, rMIN, rMAX, NR, zMIN, zMAX, NZ, flag_per, vecPER)};
					 
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  //NR=NR-2; //CAREFUL HERE!!!
    //wtime1=MPI::Wtime();
    vector<double> tcorrRG[NM];
    const int NMPAR=NM;
    for(int im=0; im<NMPAR; im++) {
      //for(int im=NMPAR; im<NM; im++) {
      wtime1=MPI::Wtime();
      tcorrRG[im]=OcorrRG[im].dual_tree_2p(dummy);
      for(int jr=0; jr<NR; jr++)
	for(int jz=0; jz<NZ; jz++)
	  npairsRG[id][im][jr][jz]=tcorrRG[im][jr*NZ+jz];
      wtime2=MPI::Wtime();
      if(id==master)
	cerr<<"time taken to compute RG correlations for mag bin: "<<im<<" is: = "<<wtime2-wtime1<<endl<<endl;
    }

    MPI::COMM_WORLD.Barrier();
        
    if(id!=master) 
      {
	for(int im=0; im<NMPAR; im++) 
	  //for(int im=NMPAR; im<NM; im++) 
	  MPI::COMM_WORLD.Send(&npairsRG[id][im], NR*NZ, MPI::DOUBLE, master, tag);
      }

    else {
      cerr<<"receiving RG from slaves: "<<endl;      
      //now the master part
      for(int np=1; np<ntasks; np++) 
	{
	  tag=np;
	  for(int im=0; im<NMPAR; im++)
	    //for(int im=NMPAR; im<NM; im++) 
	    MPI::COMM_WORLD.Recv(&npairsRG[np][im], NR*NZ, MPI::DOUBLE, np, tag, status);
	}
      cerr<<"done receiving RG"<<endl;
      //output RG files
      
      for(int np=0; np<ntasks; np++) {
	ostringstream oss;
	oss<<np;
	file_out_RG1[np] = file_RG + oss.str() + ".dat";
      }
      cerr<<"done creating output names"<<endl;

      
      for(int np=0; np<ntasks; np++) {
	cerr<<"opening file " + file_out_RG1[np]<<endl;
        ofstream out(file_out_RG1[np].c_str(), ios::out);
	cerr<<"opened file " + file_out_RG1[np]<<endl;
        for(int jr=0; jr<NR; jr++)
          for(int iz=0; iz<NZ; iz++) {
            out<<rb[1][jr]<<' '<<zb[1][iz];
            for(int im=0; im<NMPAR; im++)
              out<<' '<<mags[im]<<' '<<mags[im+1]<<' '<<npairsRG[np][im][jr][iz];
            out<<endl;
          }
        out.close();
	cerr<<"done writing file "<<file_out_RG1[np]<<endl;
      }
      
      }
    
    MPI::COMM_WORLD.Barrier();
    FinTime=MPI::Wtime();
    cerr<<"total time taken: "<<FinTime-InitTime<<endl<<endl;
    MPI::Finalize();
    
    return 0;
    
}


