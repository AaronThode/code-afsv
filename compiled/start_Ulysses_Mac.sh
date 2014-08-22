#!/bin/bash
exe_name=$0
exe_dir=`dirname "$0"`

# Matlab2012b installation path
  MATLABROOT='/Applications/MATLAB_R2013b.app'
# Default MCR install path on mac
  MCRROOT='/Applications/MATLAB/MATLAB_Compiler_Runtime/v82'
	
# Check for existence
if [[ -d "$MCRROOT" ]] ; then
  echo "Using MCR installation in "$MCRROOT""
  MYROOT="$MCRROOT"
elif [[ -d "$MATLABROOT" ]] ; then
  echo "Using Matlab installation in "$MATLABROOT""
  MYROOT="$MATLABROOT"
else
  echo "No valid Matlab or MCR installation found"
  exit 1 # terminate the script with a nonzero exit status (failure)
fi


# Run bundled script
  "${exe_dir}"/run_Ulysses_Mac.sh "$MYROOT"

exit 0
