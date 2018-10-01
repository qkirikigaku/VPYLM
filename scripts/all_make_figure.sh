#!bin/bash

array=("biliary_tract" "bone" "breast" "cervix" "endoemtrium" \
        "haematopoietic_and_lymphoid_tissue" "large_intestine" "liver" \
        "lung" "oesophagus" "pancreas" "prostate" "skin" "soft_tissue" \
        "stomach" "thyroid" "upper_aerodigestive_tract" "urinary_tract")

for i in "${array[@]}"; do
    bash scripts/make_figure.sh ${i}
done
