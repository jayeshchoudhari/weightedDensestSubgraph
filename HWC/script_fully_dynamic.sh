#!/usr/bin/env bash

tag="tags-ask-ubuntu"
r=2
time="1"

input_f="../data/${tag}/${tag}.hg.I.${time}.txt"

out_dir="../data/output/${tag}/HWC/D"
mkdir -p $out_dir

for eps in "0.1" "0.3" "0.5";do
    output_f="${out_dir}/out.HWC.eps=${eps}.${tag}.hg.D.${time}.txt"
    echo $eps $input_f $output_f
    ./cmake-build-debug/addRemove.out $eps $input_f $r > $output_f
done