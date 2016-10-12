#!/bin/sh
#
#

pdbroot=/c/cbg/mdw/pdb/pdb
prog=../src/molprep
top=$HOME/projects/dl_field/hbuild/branches/exp2/data/top.dat

rm -f exit_code.log

date

for pdb in `find $pdbroot -name '*.ent.gz'`; do
  echo
  echo "Processing $pdb"

  $prog <<_EOF
    inPDB = $pdb
    outPDB = test.pdb
    top_file = $top
    output_format = std
    remove_H = y
    write_ssbond = y
    keep_ss_name = n
_EOF

  exit_code=$?

  if [ $exit_code != 0 ]; then
    echo "$pdb: $prog exited with code $exit_code" >> exit_code.log
  fi
done

date
