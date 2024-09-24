#!/bin/bash

# Runscript for the BlueCloud pipeline
# Author: Urs Hofmann Elizondo, 01/07/2024
# Usage: run_cephalopod.sh -n 5 -t 12:00:00 -m 8G euler_test_phytobase.R cephalopod_small_test.out

# Function to display usage information
usage() {
  echo "Usage: $0 [-n CORES] [-t TIME] [-m MEM_PER_CPU] script_file output_file"
  echo "  -n CORES: Number of cores (default: 120)"
  echo "  -t TIME: Time limit (default: 24:00:00)"
  echo "  -m MEM_PER_CPU: Memory per CPU (default: 10G)"
  echo "  If everything fails, contact Alexandre Schickele (alexandre.schickele@usys.ethz.ch) to get help"
  exit 1
}

# Default values
cores=120
time="24:00:00"
mem_per_cpu="10G"

# Parse optional arguments
while getopts "n:t:m:" opt; do
  case $opt in
    n) cores=$OPTARG ;;
    t) time=$OPTARG ;;
    m) mem_per_cpu=$OPTARG ;;
    *) usage ;;
  esac
done

# Shift away parsed options, leaving the script file and output file as arguments
shift $((OPTIND -1))

# Ensure script and output file are provided
if [ $# -ne 2 ]; then
  usage
fi

script_file=$1
output_file=$2

# Remove modules
module purge

# Load all needed modules
module load stack/.2024-06-silent  gcc/12.2.0  openmpi/4.1.6
module load gdal/3.4.3
module load proj/9.2.1
module load geos/3.9.1
module load udunits/2.2.28

module load cmake/3.27.7
module load sqlite/3.43.2
module load postgresql/15.2
module load netcdf-c/4.9.2
module load libpng/1.6.39-fz4tvmr
module load libjpeg-turbo/3.0.0
module load poppler/0.79.0-7gdm7cu
module load poppler-data/0.4.12-u2qsqwg
module load r/4.3.2

# --vanilla means R starts with a clean slate, without saving or restoring anything from previous sessions
# --no-echo supresses the printing commands
sbatch -n "$cores" --time="$time" --mem-per-cpu="$mem_per_cpu" --wrap "R --vanilla --args cores $cores < $script_file > $output_file"
