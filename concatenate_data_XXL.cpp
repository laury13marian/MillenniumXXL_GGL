#include "general.h"
#include "tree_code.h"

int main(int argc, char* argv[])

{

  if(argc != 2) {
    cerr<<"please introduce task number: "<<endl;
    exit(1);
  }

  const int id = atoi(argv[1]);
  const int tasks = 32;
  const int nsampM=150000;
  const int snap=63;
  const int NS=20;
  int slice[NS] = {0, 200, 400, 600, 800, 1000, 1200, 1400, 1500, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3070, 3071, 3072};
  ostringstream oss[NS], osNSAMP, osSNAP, osNTASK;
  osNSAMP<<nsampM; osSNAP<<0; osSNAP<<snap; osNTASK<<tasks;
  for(int i=0; i<NS; i++)
    oss[i]<<slice[i];

  const string path_XXL = "/u/lmarian/LauraTreeCode/LauraMPI/XXL_DATA/";
  const string file_gen = path_XXL + "particles_XXL_snap_" + osSNAP.str() + "_nsamp_" 
    + osNSAMP.str() + "_ntasks_" + osNTASK.str() + "_task_" + argv[1] + "_slices_";
  const string file_out = path_XXL + "COMBINED/particles_XXL_snap_" + osSNAP.str() + "_nsamp_" 
    + osNSAMP.str() + "_ntasks_" + osNTASK.str() + "_task_" + argv[1] + ".dat";

  vector<particle> V;

  for(int i=0; i<NS-1; i++) {
    string file_in = file_gen + oss[i].str() + '_' + oss[i+1].str() + ".dat";
    ifstream in(file_in.c_str(), ios::in);
    check_stream(in, file_in);
    while(!in.eof()) {
      vector<float> T(3);
      in>>T[0]>>T[1]>>T[2];
      if( T[0]==0. && T[1]==0. && T[2]==0.)
	continue;
      particle OT(T);
      V.push_back(OT);
    }
    in.close();
    cerr<<"done reading and copying file " + file_in<<endl;
  }
  cerr<<"new size of particle vector from all slices: "<<V.size()<<endl;

  ofstream out(file_out.c_str(), ios::out);
  for(long unsigned i=0; i<V.size(); i++)
    out<<V[i].pos[0]<<' '<<V[i].pos[1]<<' '<<V[i].pos[2]<<endl;
  out.close();
  cerr<<"done writing particle output file for task: "<<id<<endl;
  cerr<<file_out<<endl;

  return 0;

} 
