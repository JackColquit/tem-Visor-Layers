#!/bin/bash
 #This source code is licensed under the MIT license found in the
 #MIT-LICENSE file in the root directory of this source tree.

#This file is use to check the temVisorLayers.py on linux
file1=3pab.cif
file2=scatteringTables10.csv
python3 temVisorLayers.py $file1 $file2 80000
