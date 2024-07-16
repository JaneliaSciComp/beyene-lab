# Code to do colocalization analysis of DRD1 and TDT channels

Two main scripts can be used:
1. `drd1Colocalization.m` which will do colocalization analysis on a set of DAPI, DRD1 and TDT images.
2. `drd1ColocalizationRandomShift.m` which will do the same as above, except it will do so over a set number of random trials within which the TDT image will be randomly translated.
In both cases, `doMask` can be toggled to perform the analyis with or without masking of DAPI with DRD1.

Output results are saved in a `/tifs/` directory.