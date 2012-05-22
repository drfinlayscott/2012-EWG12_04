###############################################################################
# EJ(20120507)
# Analysis of Sardinella aurita data for STECF WG 12 04 
#
###############################################################################

#==============================================================================
# Get data from fishbase
#==============================================================================

library(rfishbase)

saurita <- fish.data[findSpecies("Sardinella aurita", fish.data)]

load(url("http://www.fishbase.org/Reproduction/MaturityList.php?ID=1043&GenusName=Sardinella&SpeciesName=aurita&fc=43"))

read.table("http://www.fishbase.org/Reproduction/MaturityList.php?ID=1043&GenusName=Sardinella&SpeciesName=aurita&fc=43", skip=100)

