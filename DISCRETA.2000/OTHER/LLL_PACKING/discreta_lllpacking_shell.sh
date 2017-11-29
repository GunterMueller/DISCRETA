#!/bin/sh
discreta_deletelines $1 $2 $3 $4 $5 >infile 
discreta_lllpacking $1 $2 $3 $4 infile 
rm infile 
echo "finished."
