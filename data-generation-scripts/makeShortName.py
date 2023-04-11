# Makes panacota compatible short name from species name
# Assumes Genus_species
import sys
import re

input_name = sys.argv[1]

input_name_split = input_name.split('_')
output_name = input_name_split[0][0:2]+input_name_split[1][0].upper()+input_name_split[1][1:2]
print(output_name)
