#Shell script to update golden test files
for SUBDIRECTORY in *
	do
		echo $SUBDIRECTORY
		rm $SUBDIRECTORY/\*.csv.golden
		rm $SUBDIRECTORY/\*.tre.golden
		cp $SUBDIRECTORY/$SUBDIRECTORY.tre $SUBDIRECTORY/$SUBDIRECTORY.tre.golden
		cp $SUBDIRECTORY/$SUBDIRECTORY.csv $SUBDIRECTORY/$SUBDIRECTORY.csv.golden
		cp $SUBDIRECTORY/$SUBDIRECTORY.ia.txt $SUBDIRECTORY/$SUBDIRECTORY.ia.txt.golden
        cp $SUBDIRECTORY/$SUBDIRECTORY.dot $SUBDIRECTORY/$SUBDIRECTORY.dot.golden
        rm $SUBDIRECTORY/\*.txt*
	done
