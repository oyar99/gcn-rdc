#!/usr/bin/env python
import sys

file1 = sys.argv[1]
file2 = sys.argv[2]

file1contents = set(open(file1).readlines())
file2contents = set(open(file2).readlines())
if file1contents == file2contents:
    print("File contents are the same!")
else:
    print("In file2, not file1:\n")
    for diffLine in file2contents - file1contents:
        print("\t", diffLine)
    print("\nIn file1, not file2:\n")
    for diffLine in file1contents - file2contents:
        print ("\t", diffLine)