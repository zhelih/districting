#!/usr/bin/python
import sys
f = open(sys.argv[1], "r")
l = f.read()
f.close()
(a,b,c) = l.split("\n", 2)
f = open(sys.argv[1], "w")
f.write("\n".join([b,a,c]))
f.close()
