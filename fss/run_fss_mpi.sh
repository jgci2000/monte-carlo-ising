#!/bin/bash

#
# Script to run FSS method
# João Inácio, Apr. 17th, 2021
#

# File name
input_file=$1

# Number of processes
shift
n_cores=$1

# Runtime arguments
shift
args="$*"

# Run the program
mpiexec -np $n_cores ./$input_file $args
exit $?

exit $rc
