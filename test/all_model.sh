#!/bin/sh
#
# Copyright (C) 2011 Hannes Loeffler, STFC Daresbury, UK
#
# simple shell script to add hydrogens to all models of a PDB file
#
# $Id: all_model.sh 86 2012-05-10 12:18:45Z hhl $
#



inpdb=2KJJ.pdb
outpdb=all_model.pdb
prog=../src/molprep
start=1
end=20


rm -f $outpdb

for n in `seq $start $end`; do
  echo "Processing number $n"

  $prog <<_EOF
    inPDB = $inpdb
    outPDB = ${n}_$inpdb
    top_file = ../data/top.dat
    model_no = $n
    no_cryst_record = y
    no_ter_record = y
    no_end_record = y
    remove_H	= y
_EOF

  cat ${n}_$inpdb >> $outpdb
  rm ${n}_$inpdb
done
