#!/usr/bin/env bash
for item in Batch*
do
rm -rfv -- $item/final/*/*/final.*.vcf ###Batch1/final/PRAD-CA/DO51159/final.*.vcf gibi
rm -rfv -- $item/final/randombeds
done


