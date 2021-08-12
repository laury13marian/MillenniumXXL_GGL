#include "read_data.h"

int main()

{

  string path_in = "/galformod/data/mxxl/";
  int snapshot = 63;
  int subsnap = 3071;
  //const signed long rng_seed_gen = -static_cast<signed long>(std::time(NULL)); 
  const int nsamp=150000;
  const signed long rng_seed_gen = -1381718569;
  read_data OB(path_in, snapshot);
  OB.read_particles(subsnap);
  vector<particle> V=OB.stream_particles(rng_seed_gen, nsamp);
  cerr<<V.size()<<endl;

  return 0;

}

