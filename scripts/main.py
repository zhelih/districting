#!/usr/bin/python3
import sys

def read_hash(filename):
  return dict([line.replace("\n", "").split(" ")[::-1] for line in open(filename).readlines()])

def read_table(filename_in, filename_out, id_map):
  fin = open(filename_in)
  fout = open(filename_out, "w")
  header = fin.readline().replace("\n", "").split(",")[1:]
  header = [id_map[i] for i in header]
  sorted_header = list(map(str, sorted(map(int, header))))
  #print(list(sorted_header))
  print("Node_ID,"+",".join(sorted_header), file = fout)
  data = {}
  for i in header:
    data[i] = {}

  for line in fin.readlines():
    line = line.replace("\n", "").split(",")
    hashed = line[0].replace('"', "")
    values = line[1:]
    lid = id_map[hashed]
    for i in range(len(values)):
      data[header[i]][lid] = '{:d}'.format(round(float(values[i]) / 1000.0))

  for iid in sorted_header:
    values = iid+","+",".join([data[iid][jid] for jid in sorted_header])
    print(values, file = fout)


if __name__ == "__main__":
  read_table(sys.argv[2], sys.argv[3], read_hash(sys.argv[1]))
