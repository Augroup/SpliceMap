#!/bin/sh
if [ $# -eq 0 ]
then
echo "Usage: $0 output_directory"
echo "Builds the (64-bit) SpliceMap package [including Barloader] and copies executables to output_directory and overwrites the old files"
exit 1
fi

cd SpliceMap-src

echo ">>> Building the (64-bit) SpliceMap package..."

make install

cp SpliceMap runSpliceMap sortsam nnrFilter neighborFilter uniqueJunctionFilter randomJunctionFilter wig2barwig colorJunction subseq findNovelJunctions statSpliceMap countsam amalgamateSAM precipitateSAM ..

echo ">>> Cleaning package..."

make clean 2> /dev/null

cd ..

echo ">>> Copying the SpliceMap package to $1"

mv -v SpliceMap runSpliceMap sortsam nnrFilter neighborFilter uniqueJunctionFilter randomJunctionFilter wig2barwig colorJunction subseq findNovelJunctions statSpliceMap countsam amalgamateSAM precipitateSAM $1



cd barloader-src

echo ">>> Building the (64-bit) Barloader..."

g++ -m64 -O3 -o barloader *.cpp

cp barloader ..

echo ">>> Cleaning package..."

rm barloader

cd ..

echo ">>> Copying Barloader to $1"

mv -v barloader $1






