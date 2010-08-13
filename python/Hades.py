#!/usr/bin/env python

# Code made by Bruno Turcksin
# Interface of the code Acheron

import os
import sys

from data import *

try :
    file_path = sys.argv[1]
    a = data(file_path)
    a.create_mesh()

except IndexError :
    print ("You need to give an input file")
except IOError :
    print ("Cannot open the input file")
