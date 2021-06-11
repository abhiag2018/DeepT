#!/usr/bin/env bash

gzip -df $input_feature_part_gz

i=`python ${projectDir}/templates/convert_tar_xz.py $cellType $input_feature_part`

tar cJvfh features.${cellType}.\${i}.tar.xz data_\${i}

echo \$i
rm $input_feature_part



