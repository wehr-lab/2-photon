#!/bin/bash -x

source /etc/profile.d/modules.sh

module load python3
echo "Done!"

conda activate 2-photon 

which python 