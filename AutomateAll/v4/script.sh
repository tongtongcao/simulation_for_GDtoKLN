#!/bin/csh -f
setenv ENVIRONMENT $1
setenv RUNNO $2
setenv RUNFILE $3

setenv PPARAM 0x7e
setenv APARAM 1.0
setenv BPARAM 1.0
setenv CPARAM 1.0
setenv FPARAM 1.0

source $ENVIRONMENT

echo pwd = $PWD
printenv

#generator
./generator -R generated_${RUNFILE}.root -M 50000 -E 0.9:2.6 
echo "LS:"
ls -al

#PutInBOS
./BOSWrite generated_${RUNFILE}.root
echo "LS:"
ls -al

# gsim 
gsim_bat -ffread gsim.ffread -mcin generated_${RUNFILE}.evt -bosout gsim_file.bos -tg /group/clas/parms/bgrid_T67to33.fpk
echo "LS:"
ls -al

# gpp
gpp -P${PPARAM} -a${APARAM} -b${BPARAM} -c${CPARAM} -f${FPARAM} -R${RUNNO} -s -Y -oinfile gsim_file.bos
h2root gpp.hbook
echo "LS:"
ls -al

# user_ana
user_ana -t cook.tcl
h2root anamonhist
mv outfile cooked_${RUNFILE}.bos
echo "LS:"
ls -al

# skim
./skim_gsim cooked_${RUNFILE}.bos
echo "LS:"
ls -al
