#!/bin/bash
maped=/mn/stornext/u3/trygvels/compsep/cmbco/src/f90/map_editor/map_editor

map2png=/mn/stornext/u3/hke/local/bin/map2png

simdir=rmsvalues2

freq="040 050 060 068_1 068_2 078_1 078_2 089_1 089_2 100_1 100_2 119_1 119_2 140_1 140_2 166 195_1 195_2 235 280 337 402"


for f in $freq; do

    $maped print_stats LB_v28_noise_256_${f}.fits  > stats/stats_${f}.txt
    rms1=$(sed -n 17p stats/stats_${f}.txt | gawk -F ' ' '{print $3}')
    rms2=$(sed -n 25p stats/stats_${f}.txt | gawk -F ' ' '{print $3}')
    rms=$(echo $rms1 + $rms2 | bc) #dosent work with E-002!
#    echo $rms
#    rms=$(python -c "print ($rms1+$rms2)/2.")
    echo $f,$rms1,$rms2, $rms

    $maped scale LB_v28_noise_256_${f}.fits LB_v28_rms_256_${f}.fits 0
    $maped add_offset LB_v28_rms_256_${f}.fits LB_v28_rms_256_${f}.fits ${rms1}
    $map2png LB_v28_rms_256_${f}.fits -bar
    $map2png LB_v28_noise_256_${f}.fits -bar
done

