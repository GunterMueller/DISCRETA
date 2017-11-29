#!/bin/sh 
TMPDIR="makeKM_tmp"
echo $TMPDIR
if [ -d $TMPDIR ] ; then
	echo "rm -f $TMPDIR/*"
	rm -f $TMPDIR/*
else
	echo "mkdir $TMPDIR"
	mkdir $TMPDIR
fi

echo "cp $1 $TMPDIR"
cp $1 $TMPDIR

echo "cp $DISCRETA_HOME/OTHER/MAKE_KM/SRC/* $TMPDIR"
cp $DISCRETA_HOME/OTHER/MAKE_KM/SRC/* $TMPDIR

echo "cd $TMPDIR"
cd $TMPDIR

echo "make"
make kmgen pres
#chmod u+x makeKM

echo "makeKM $*"
makeKM $*

cd ..
mv $TMPDIR/$5 .

echo "rm -rf $TMPDIR"
rm -rf $TMPDIR

echo "finished"
