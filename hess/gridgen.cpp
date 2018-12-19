#include <cstdio>
#include <cmath>
#include <string>

#define OPEN_AND_CHECK(fn) \
{ \
  f = fopen((fn).c_str(), "w"); \
  if(!f) { \
   printf("Failed to open %s\n", (fn).c_str()); \
   return 1; \
  }; \
}

int main(int argc, char* argv[])
{
  if(argc < 4)
  {
    printf("Usage: %s <output dir> <n> <m>\n\tGenerating grid n x m and asking for n districts, equal population everywhere\n", argv[0]);
    return 0;
  }

  std::string outdir = argv[1];
  int n = std::atoi(argv[2]);
  int m = std::atoi(argv[3]);

  if(n < 1 || m < 1 || n > 100000 || m > 100000)
  {
    printf("Bad values for n/m\n");
    return 1;
  }

  std::string prefix = outdir + "/grid_" + argv[2] + "_" + argv[3];
  FILE* f = 0;

  // dimacs
  OPEN_AND_CHECK(prefix+".dimacs");
  fprintf(f, "p edge %d %d\n", m*n, 2*m*n - m - n);
  for(int i = 0; i < m*n; ++i)
  {
    if( (i+1) % m != 0) // all but the last column
      fprintf(f, "e %d %d\n", i, i+1);
    if( i < n*m - m ) // all but the last row)
      fprintf(f, "e %d %d\n", i, i+m);
  }
  fclose(f);
  // _distances.csv
  OPEN_AND_CHECK(prefix+"_distances.csv");
  // header
  fprintf(f, "Node_ID");
  for(int i = 0; i < m*n; ++i)
    fprintf(f, ",%d", i);
  for(int i = 0; i < m*n; ++i)
    for(int j = 0; j < m*n; ++j)
    {
      if(j == 0)
        fprintf(f, "\n%d", i);
      int delta_x = abs((i % m) - (j % m));
      int delta_y = abs((i % n) - (j % n));
      int dist = static_cast<int>(sqrt(delta_x*delta_x + delta_y*delta_y)*100.);
      fprintf(f, ",%d", dist);
    }
  fclose(f);
  // .hash
  // .population
  OPEN_AND_CHECK(prefix+".population");
  fprintf(f, "total pop = %d\n", m*n);
  for(int i = 0; i < m*n; ++i)
    fprintf(f, "%d 1\n", i);
  fclose(f);

  return 0;
}
