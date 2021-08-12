#include "read_data.h"

int main(int argc, char *argv[])

{

  /*const string path_in = "/ptmp/mpa/reangulo/MXXL_Gals/";
  const int snap=24;
  read_galaxy Ob(snap, path_in);
  Ob.open_subfile(0);*/
  

  if(argc!=2) {
    cerr<<"introduce the snapshot number"<<endl;
    exit(1);
  }
  int snap=atoi(argv[1]);
  ostringstream osSNAP;
  osSNAP<<0; osSNAP<<snap;
  vector<string> path_gen(2);
  path_gen[0] = 
    "/ptmp/mpa/rosm/Millennium/Millennium-XXL/subcubes/subcubes_" + osSNAP.str() + '/';
  path_gen[1] = "/snapshot_" + osSNAP.str() + "_SubCube_"; //Grp_";
  const string SP1 = "GROUP";//"PARTICLE"; //"GROUP";
  const string SP2 = "PARTICLE";
  const string SP = "GALAXY";
  //const signed long rng_seed = -static_cast<signed long>(std::time(NULL));
  const signed long rng_seed = -1387584831;
  cerr<<"seed is "<<rng_seed<<endl<<endl;
  const int Nmags=10;
  vector<double> mags(Nmags);
  mags[0]=-30.; mags[1]=-22.; mags[2]=-21.; mags[3]=-20.; mags[4]=-19.; mags[5]=-18.;
  mags[6]=-17.; mags[7]=-16.; mags[8]-=15.; mags[9]=-14.;
  const int iMAG=1;
  const int nsamp=1;
  //int NBM=11;
  //vector<double> MASSB(NBM), VGM[NBM];
  //MASSB[6]=1.e+17; MASSB[5]=5.e+13; MASSB[4]=1.e+13; MASSB[3]=5.e+12; 
  //MASSB[2]=1.e+12; MASSB[1]=5.e+11; MASSB[0]=1.e+11; 
  
  int subcube=0;
  //load_SScube OB2(subcube, path_gen, SP, nsamp, rng_seed);
  //OB2.read_species(0);
  
  load_SScube OB1(subcube, path_gen, SP, rng_seed, mags, iMAG);
  OB1.read_species(0);
  //OB1.format_data();
  vector<particle> VV1=OB1.return_formatted_data(5);
  for(int i=0; i<50; i++)
    cerr<<VV1[i].pos[0]<<' '<<VV1[i].pos[1]<<' '<<VV1[i].pos[2]<<endl;
  cerr<<"size of formatted vector: "<<VV1.size()<<endl;
   
  return 0;

}
