#include "general.h"

using namespace std;

int main()

{

const int Nsub=216;

  /*const string path_gen = "/Users/lmarian/W_I_P/MillenniumXXL/P_CORR/MM/";
  const string file_in1 = "CFw_snap_054_sub_";
  const string file_in2 = "_binsR_30_binsZ_10_NRAN_1000000_nsamp_4000_zmax_100_NT_32.dat"; 

  for(int sub=0; sub<Nsub; sub++) {
    vector<float> B; 
    vector<double> V1, V2, V3;
    ostringstream oss;
    oss<<sub;
    string file_in = path_gen + "V0/" + file_in1 + oss.str() + file_in2;
    ifstream in(file_in, ios::in);
    check_stream(in, file_in);
    while(!in.eof()) {
      double d1, d2, d3; 
      float f1;
      in>>f1>>d1>>d2>>d3;
      if(f1>0.03) { 
	B.push_back(f1);
	V1.push_back(d1);
	V2.push_back(d2);
	V3.push_back(d3);
      }
    }
    in.close();

    string file_out = path_gen + file_in1 + oss.str() + file_in2;

    ofstream out(file_out, ios::out);
    for(int i=0; i<B.size()-1; i++) 
      out<<B[i]<<' '<<V1[i]<<' '<<V2[i]<<' '<<V3[i]<<endl;
    out.close();

    cerr<<"done writing new output file " + file_out <<endl;
    }*/

  /*const string path_gen = "/Users/lmarian/W_I_P/MillenniumXXL/P_CORR/GG/";
  const string file_in1 = "CFw_snap_054_sub_";
  const string file_in2 = "_binsM_9_binsR_30_binsZ_10_NRAN_1000000_zmax_100_NT_32.dat";*/

 const string path_gen = "/Users/lmarian/W_I_P/MillenniumXXL/P_CORR/GM/";
 const string file_in1 = "CFw_snap_054_sub_";
 const string file_in2 = "_binsM_9_binsR_30_binsZ_10_NRAN_1000000_nsamp_4000_zmax_100_NT_32.dat";

  const int NMAG=9; const int NMAG2=4;
  for(int sub=0; sub<Nsub; sub++) {
    vector<float> B; 
    vector<double> V1[NMAG], V2[NMAG], V3[NMAG], M1[NMAG], M2[NMAG];
    ostringstream oss;
    oss<<sub;
    string file_in = path_gen + "V0/" + file_in1 + oss.str() + file_in2;
    ifstream in(file_in, ios::in);
    check_stream(in, file_in);
    while(!in.eof()) {
      float f1;
      in>>f1;
      if(f1>0.03) B.push_back(f1);
      for(int im=0; im<NMAG; im++) {
	double d1, d2, d3, m1, m2;
	in>>m1>>m2>>d1>>d2>>d3;
	if(f1>0.03) { 
	  V1[im].push_back(d1);
	  V2[im].push_back(d2);
	  V3[im].push_back(d3);
	  M1[im].push_back(m1);//=m1;
	  M2[im].push_back(m2); //=m2;
	}
      }
    }
    in.close();
  
    string file_out = path_gen + file_in1 + oss.str() + file_in2;

    ofstream out(file_out, ios::out);
    for(int i=0; i<B.size()-1; i++) { 
      out<<B[i];
      for(int im=0; im<NMAG2; im++)
	out<<' '<<M1[im][i]<<' '<<M2[im][i]<<' '<<V1[im][i]<<' '<<V2[im][i]<<' '<<V3[im][i];
      out<<endl;
    }
    out.close();

    cerr<<"done writing new output file " + file_out <<endl;
  }


  return 0;

}
    


