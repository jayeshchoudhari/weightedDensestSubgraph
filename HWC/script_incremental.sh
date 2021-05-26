#!/usr/bin/env bash

tag="tags-stack-overflow"
r="5"
time="1"

input_f="../data/${tag}/${tag}.hg.I.${time}.txt"

out_dir="../data/output/${tag}/HWC/I"
mkdir -p $out_dir

for eps in "0.1" "0.3" "0.5";do
    output_f="${out_dir}/out.HWC.r=${r}.eps=${eps}.${tag}.hg.I.${time}.txt"
    echo $eps $input_f $output_f
    ./cmake-build-debug/add.out $eps $input_f $r > $output_f
done
