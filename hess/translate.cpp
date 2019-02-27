#include <cstdio>
#include <vector>

using namespace std;

void scan_file(vector<int>& out, const char* fname)
{
  //TODO newbie checks?
  FILE *f = fopen(fname, "r");
  if(!f)
  {
    printf("Cannot open %s!\n", fname);
    return;
  }

  int dummy;
  while(fscanf(f, "%d ", &dummy) != EOF)
  {
    int data;
    fscanf(f, "%d ", &data);
    out.push_back(data);
  }
  fclose(f);
}

int main(int argc, char* argv[])
{
  if(argc < 3)
  {
    printf("Translate the solution for QGIS drawing script\n");
    printf("Usage: %s <.out> <.hash> [output]\n", argv[0]);
    return 0;
  }

  vector<int> geoids;
  vector<int> districts;

  scan_file(districts, argv[1]);
  scan_file(geoids, argv[2]);

  if(geoids.size() != districts.size())
  {
    printf("Sizes do not match: hash has %lu while result has %lu!\n", geoids.size(), districts.size());
    return 1;
  }

  printf("Size: %lu\n", geoids.size());

  FILE* out;
  if(argc > 3)
    out = fopen(argv[3], "w");
  else
    out = fopen("translation.out", "w");

  for(unsigned int i = 0; i < geoids.size(); ++i)
    fprintf(out, "%d %d\n", geoids[i], districts[i]);

  fclose(out);

  printf("Done\n");

  return 0;
}
