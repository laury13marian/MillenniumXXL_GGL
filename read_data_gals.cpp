#include "read_data.h"

using namespace std;

int main()

{

  const string path_in = "/galformod/scratch/bmh20/SAM/Henriques2014/snaps/Guo10_bug_fix/MR/";
  const double red=0.24;
  read_galaxy_Millennium OB(red, path_in);
  int sub=100;
  OB.open_subfile(sub);
  OB.deallocate(sub);
  return 0;

}
