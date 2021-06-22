#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --output=#{OUT_DIR}/output.log
#SBATCH --error=#{OUT_DIR}/output.err
#SBATCH --mem-per-cpu=2000
#SBATCH --partition=haswell
#SBATCH --mail-type=end
#SBATCH --mail-user=simon.praetorius@tu-dresden.de
#SBATCH -J #{POSTFIX}
#SBATCH -A wir

#{COMMAND}
