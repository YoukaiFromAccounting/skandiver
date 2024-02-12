#!/bin/bash

#Switch working directory to active script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$DIR"

EXAMPLE_DIR="test_files"
EXAMPLE_FILE="GCA_018123985.1_PDT001010899.1_genomic.fna"
EXAMPLE_URL="https://github.com/YoukaiFromAccounting/skandiver/tree/main/test_files"

#Check for Python Installation
if ! [ -x  "$(command -v python3)" ]; then
	echo "Python3 not installed; please install Python 3 before continuing."
	exit 1
fi 

#Import Python packages
python3 -c "import pandas, Bio"
if [ "$?" -eq 1 ]; then
	echo "Python packages pandas and/or Bio not found. Please install."
	echo "e.g., if using pip installer, use \"pip3 install biopython pandas\""
	echo "Note, if you are using Ubuntu and do not have pip installed, you can install it using \"sudo apt-get install python3-pip\"."
	echo "Otherwise, see instructions here: https://pip.pypa.io/en/stable/installing/"
    exit 1
fi

#Check for skani installation
if ! [ -x "$(command -v skani)" ]; then
	echo "Error: skani not found. Please install."
	echo "e.g., if using conda installer, use \"conda install -c bioconda skani\""
	exit 1
fi

#Download example datasets
if [ -e "test_files/abaumannii/abaumanniiWGS.fna" ]; then
	echo "Examples already installed"
else
	mkdir -p "${EXAMPLE_DIR}"
	wget -c "${EXAMPLE_URL}" || curl -0 "${EXAMPLE_URL}"
	echo "Examples installed successfully"
fi

#Check if sufficient memory
totalk=$(awk '/^MemTotal:/{print $2}' /proc/meminfo)

if [ "$?" -ne 0 ]; then
	echo "WARNING: Cannot access /proc/meminfo to check for system memory."
	echo "skani's default parameters may not work as intended with less than 32 GB of RAM."
elif [ "$totalk" -lt 3355000 ]; then
	echo "WARNING: Operating system has less than 32 GB of RAM."
	echo "skani may run out of memory during queries using standard parameters."
else
	echo "You have sufficient memory (at least 32 GB) for this installation."
fi

#Complete setup phase
echo ========================================
echo Installation complete!
echo Test out skandiver on a sample whole genome assembly with:
echo ./skandiver.sh test_files 10000 [PATH_TO_REFERENCE_GENOME_DATABASE]
