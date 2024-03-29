#!/usr/bin/python
import sys

def read_hash(filename):
  return dict([line.replace("\n", "").split(" ")[::-1] for line in open(filename).readlines()])

def read_sol(filename):
  return dict([line.replace("\n", "").split(" ") for line in open(filename).readlines()])

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: %s <hash> <sol> -> outputs blk to stdout" % sys.argv[0])
        exit(0)
    h = read_hash(sys.argv[1])
    out = read_sol(sys.argv[2])

    for k in h:
        print("%s %s" % (k, out[h[k]]))
