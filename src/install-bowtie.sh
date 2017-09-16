#!/bin/sh
if [ $# -eq 0 ]
then
echo "Usage: $0 output_directory"
echo "Builds the (64-bit) Bowtie package and copies executables to output_directory and overwrites the old files"
exit 1
fi

cd bowtie-0.12.7

echo ">>> Building the (64-bit) Bowtie package..."

make BITS=64 2> /dev/null

cp bowtie bowtie-build bowtie-inspect ..


echo ">>> Cleaning package..."

make clean 2> /dev/null


cd .. 

echo ">>> Copying the Bowtie package to $1"
mv -v bowtie bowtie-build bowtie-inspect $1


