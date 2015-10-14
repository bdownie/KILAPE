CC=g++
#CFLAGS= -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -pg
CFLAGS= -O3 -Wall 
#INCLUDE=-I/home/lakatos/bdownie/inc
#BOOST_INCLUDE=-I/home/lakatos/bdownie/inc/boost
OMP_FLAGS=-fopenmp -lgomp
#OMP_LINK_FLAGS=-lgomp

all: adjust_graph PrepareLocalAssemblies BuildGraphFromSAM BuildScaffolds Kmasker filterSamOnCount countSamSites RestoreFullContigs scaffolder_int DistributeOrphanedContigs
	mkdir -p ../bin
	mv scaffolder_int adjust_graph BuildScaffolds BuildGraphFromSAM PrepareLocalAssemblies Kmasker filterSamOnCount countSamSites RestoreFullContigs DistributeOrphanedContigs ../bin

adjust_graph: adjust_graph.cpp
	${CC} ${CFLAGS} ${INCLUDE} -o adjust_graph adjust_graph.cpp 

DistributeOrphanedContigs: DistributeOrphanedContigs.cpp
	${CC} ${CFLAGS} ${INCLUDE} -o DistributeOrphanedContigs DistributeOrphanedContigs.cpp 

BuildScaffolds: BuildScaffolds.cpp
	${CC} ${CFLAGS} ${INCLUDE} -o BuildScaffolds BuildScaffolds.cpp 

BuildGraphFromSAM: BuildGraphFromSAM.cpp
	${CC} ${CFLAGS} ${INCLUDE} -o BuildGraphFromSAM BuildGraphFromSAM.cpp 

PrepareLocalAssemblies: PrepareLocalAssemblies.cpp
	${CC} ${CFLAGS} ${OMP_FLAGS} ${INCLUDE} -o PrepareLocalAssemblies PrepareLocalAssemblies.cpp 

Kmasker: Kmasker.cpp
	${CC} ${CFLAGS} ${OMP_FLAGS} ${INCLUDE} -o Kmasker Kmasker.cpp 

RestoreFullContigs: RestoreFullContigs.cpp
	${CC} ${CFLAGS} ${INCLUDE} -o RestoreFullContigs RestoreFullContigs.cpp 

filterSamOnCount: filterSamOnCount.cpp
	${CC} ${CFLAGS} ${INCLUDE} -o filterSamOnCount filterSamOnCount.cpp 

countSamSites: countSamSites.cpp
	${CC} ${CFLAGS} ${INCLUDE} -o countSamSites countSamSites.cpp 

scaffolder_int: scaffolder.cpp
	${CC}  -DID_TYPE=2 ${CFLAGS} ${INCLUDE} -o scaffolder_int scaffolder.cpp 

clean: 
	rm -rf *.d *.o 
	rm -rf ../Kmasker ../PrepareLocalAssemblies ../BuildGraphFromSAM ../BuildScaffolds ../adjust_graph ../countSamSites
