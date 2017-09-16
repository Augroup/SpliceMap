#!/bin/sh
if [ $# -eq 0 ]
then
echo "Usage: $0 output_directory"
echo "Builds the SpliceMap package and copies executables to output_directory and overwrite the old files"
exit 1
fi

echo ">>> Cleaning package..."

make clean 2> /dev/null

echo ">>> Building the (64-bit) SpliceMap package..."

make install

echo ">>> Copying the SpliceMap package to $1"

cp -v SpliceMap runSpliceMap colorJunction findNovelJunctions neighborFilter sortsam nnrFilter randomJunctionFilter wig2barwig statSpliceMap subseq uniqueJunctionFilter countsam amalgamateSAM precipitateSAM $1





