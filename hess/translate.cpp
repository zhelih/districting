#include <cstdio>
#include <vector>
#include <string>

using namespace std;

void scan_file(vector<string>& out, const char* fname)
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
    char data[1023];
    fscanf(f, "%s ", data);
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

  vector<string> geoids;
  vector<string> districts;

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
    fprintf(out, "%s %s\n", geoids[i].c_str(), districts[i].c_str());

  fclose(out);

  printf("Done\n");

  return 0;
}
