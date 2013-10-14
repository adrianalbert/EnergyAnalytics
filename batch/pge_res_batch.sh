#!/bin/sh
#cat zips.txt | time parallel -u -j-2 --eta echo {}
cat zips.txt | time parallel -u -j-2 --eta ./pge_res_batch_process.R {}
