CC = g++
CFLAGS-32 = -m32 -O3 -Wall
CFLAGS-64 = -m64 -O3 -Wall

install: SpliceMap

install-32: SpliceMap-32

clean: 
	rm *.o SpliceMap runSpliceMap sortsam nnrFilter neighborFilter uniqueJunctionFilter randomJunctionFilter wig2barwig colorJunction subseq findNovelJunctions statSpliceMap countsam amalgamateSAM precipitateSAM

SpliceMap: runSpliceMap sortsam neighborFilter nnrFilter uniqueJunctionFilter randomJunctionFilter wig2barwig colorJunction subseq findNovelJunctions statSpliceMap countsam amalgamateSAM precipitateSAM params.o cfgfile.o ostools.o SpliceMap_utils.o main.cpp main.h
	$(CC) $(CFLAGS-64) -o SpliceMap main.cpp SpliceMap_utils.o ostools.o params.o cfgfile.o

runSpliceMap: SpliceMap_utils.o params.o cfgfile.o ostools.o runSpliceMap.cpp runSpliceMap.h 
	$(CC) $(CFLAGS-64) -o runSpliceMap SpliceMap_utils.o ostools.o cfgfile.o runSpliceMap.cpp params.o

findNovelJunctions: SpliceMap_utils.o ostools.o params.o findNovelJunctions.cpp findNovelJunctions.h
	$(CC) $(CFLAGS-64) -o findNovelJunctions findNovelJunctions.cpp SpliceMap_utils.o params.o ostools.o

statSpliceMap: SpliceMap_utils.o ostools.o params.o statSpliceMap.cpp statSpliceMap.h
	$(CC) $(CFLAGS-64) -o statSpliceMap statSpliceMap.cpp SpliceMap_utils.o params.o ostools.o

subseq: subseq.cpp subseq.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o subseq subseq.cpp SpliceMap_utils.o params.o ostools.o

wig2barwig: wig2barwig.cpp wig2barwig.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o wig2barwig wig2barwig.cpp SpliceMap_utils.o params.o ostools.o

colorJunction: colorJunction.cpp colorJunction.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o colorJunction colorJunction.cpp SpliceMap_utils.o ostools.o params.o

uniqueJunctionFilter: uniqueJunctionFilter.cpp uniqueJunctionFilter.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o uniqueJunctionFilter uniqueJunctionFilter.cpp SpliceMap_utils.o ostools.o params.o

randomJunctionFilter: randomJunctionFilter.cpp randomJunctionFilter.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o randomJunctionFilter randomJunctionFilter.cpp SpliceMap_utils.o ostools.o params.o

neighborFilter: neighborFilter.cpp neighborFilter.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o neighborFilter neighborFilter.cpp SpliceMap_utils.o ostools.o params.o

nnrFilter: nnrFilter.cpp nnrFilter.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o nnrFilter nnrFilter.cpp SpliceMap_utils.o ostools.o params.o

sortsam: sortsam.cpp sortsam.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o sortsam sortsam.cpp SpliceMap_utils.o ostools.o params.o

countsam: countsam.cpp countsam.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o countsam countsam.cpp SpliceMap_utils.o ostools.o params.o

amalgamateSAM: amalgamateSAM.cpp amalgamateSAM.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o amalgamateSAM amalgamateSAM.cpp SpliceMap_utils.o ostools.o params.o

precipitateSAM: precipitateSAM.cpp precipitateSAM.h SpliceMap_utils.o ostools.o params.o
	$(CC) $(CFLAGS-64) -o precipitateSAM precipitateSAM.cpp SpliceMap_utils.o ostools.o params.o

SpliceMap_utils.o: SpliceMap_utils.cpp SpliceMap_utils.h
	$(CC) $(CFLAGS-64) -c SpliceMap_utils.cpp
   
ostools.o: ostools.cpp ostools.h
	$(CC) $(CFLAGS-64) -c ostools.cpp

params.o: params.cpp params.h
	$(CC) $(CFLAGS-64) -c params.cpp

cfgfile.o: cfgfile.cpp cfgfile.h
	$(CC) $(CFLAGS-64) -c cfgfile.cpp



SpliceMap-32: runSpliceMap-32 sortsam-32 nnrFilter-32 neighborFilter-32 uniqueJunctionFilter-32 randomJunctionFilter-32 wig2barwig-32 colorJunction-32 subseq-32 findNovelJunctions-32 statSpliceMap-32 countsam-32 amalgamateSAM-32 precipitateSAM-32 params.o-32 cfgfile.o-32 SpliceMap_utils.o-32 ostools.o-32 main.cpp main.h
	$(CC) $(CFLAGS-32) -o SpliceMap main.cpp SpliceMap_utils.o ostools.o params.o cfgfile.o

runSpliceMap-32: SpliceMap_utils.o-32 params.o-32 cfgfile.o-32 ostools.o-32 runSpliceMap.cpp runSpliceMap.h 
	$(CC) $(CFLAGS-32) -o runSpliceMap SpliceMap_utils.o ostools.o cfgfile.o runSpliceMap.cpp params.o

findNovelJunctions-32: SpliceMap_utils.o-32 ostools.o-32 params.o-32 findNovelJunctions.cpp findNovelJunctions.h
	$(CC) $(CFLAGS-32) -o findNovelJunctions findNovelJunctions.cpp SpliceMap_utils.o ostools.o params.o

statSpliceMap-32: SpliceMap_utils.o-32 ostools.o-32  params.o-32 statSpliceMap.cpp statSpliceMap.h
	$(CC) $(CFLAGS-32) -o statSpliceMap statSpliceMap.cpp SpliceMap_utils.o ostools.o  params.o

subseq-32: subseq.cpp subseq.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o subseq subseq.cpp SpliceMap_utils.o ostools.o  params.o

wig2barwig-32: wig2barwig.cpp wig2barwig.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o wig2barwig wig2barwig.cpp SpliceMap_utils.o ostools.o  params.o

colorJunction-32: colorJunction.cpp colorJunction.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o colorJunction colorJunction.cpp SpliceMap_utils.o ostools.o  params.o

uniqueJunctionFilter-32: uniqueJunctionFilter.cpp uniqueJunctionFilter.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o uniqueJunctionFilter uniqueJunctionFilter.cpp SpliceMap_utils.o ostools.o  params.o

randomJunctionFilter-32: randomJunctionFilter.cpp randomJunctionFilter.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o randomJunctionFilter randomJunctionFilter.cpp SpliceMap_utils.o ostools.o  params.o

neighborFilter-32: neighborFilter.cpp neighborFilter.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o neighborFilter neighborFilter.cpp SpliceMap_utils.o ostools.o  params.o

nnrFilter-32: nnrFilter.cpp nnrFilter.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o nnrFilter nnrFilter.cpp SpliceMap_utils.o ostools.o  params.o

sortsam-32: sortsam.cpp sortsam.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o sortsam sortsam.cpp SpliceMap_utils.o ostools.o  params.o

countsam-32: countsam.cpp countsam.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o countsam countsam.cpp SpliceMap_utils.o ostools.o  params.o

amalgamateSAM-32: amalgamateSAM.cpp amalgamateSAM.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o amalgamateSAM amalgamateSAM.cpp SpliceMap_utils.o ostools.o  params.o

precipitateSAM-32: precipitateSAM.cpp precipitateSAM.h SpliceMap_utils.o-32 ostools.o-32  params.o-32
	$(CC) $(CFLAGS-32) -o precipitateSAM precipitateSAM.cpp SpliceMap_utils.o ostools.o  params.o

SpliceMap_utils.o-32: SpliceMap_utils.cpp SpliceMap_utils.h
	$(CC) $(CFLAGS-32) -c SpliceMap_utils.cpp
    
ostools.o-32: ostools.cpp ostools.h
	$(CC) $(CFLAGS-32) -c ostools.cpp

params.o-32: params.cpp params.h
	$(CC) $(CFLAGS-32) -c params.cpp

cfgfile.o-32: cfgfile.cpp cfgfile.h
	$(CC) $(CFLAGS-32) -c cfgfile.cpp


