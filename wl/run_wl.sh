#!/bin/bash

#
# Script to compile and run WL method
# João Inácio, Mar. 25th, 2021
#

# File names
input_file=$1

# Runtime arguments
shift
args="$*"

# Run the program
./$input_file $args
exit $?

exit $rc
