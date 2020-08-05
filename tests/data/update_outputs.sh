#!/bin/bash

rm -fr MAPK8_output
thoraxe -i MAPK8 -o MAPK8_output -l species_list.txt -y

rm -fr POLR3B_output
thoraxe -i POLR3B -o POLR3B_output -y
