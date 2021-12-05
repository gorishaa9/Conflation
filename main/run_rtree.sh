P=../index
#P2=/home/osboxes/libhilbert/src
P2=/scratch/gorishaa/kvpairs/libhilbert/src
g++ -w -g -I /usr/local/include -o read OsmRead.cpp $P/Insert.cc $P/RTree.cc $P/Rectangle.cc $P/Node.cc $P/NodeEntry.cc $P/RTreeHelper.cc $P/NonLeafEntry.cc $P/LeafEntry.cc $P/HilbertValue.cc $P/Common.cc $P2/BigBitVec.o $P2/FixBitVec.o $P2/Hilbert.o -lpthread -lz -lexpat -lbz2 -lgeos -I/cs/local/lib/pkg/libhilbert-0.2.1/include -I/cs/local/lib/pkg/libosmium-2.15.5/include -I/scratch/gorishaa/include/protozero/include -I/usr/include/boost
