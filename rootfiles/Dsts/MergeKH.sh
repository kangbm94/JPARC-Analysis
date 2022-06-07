DKH="DstKuramaHodoscope0"
addstring=""
space=" "
for i in {5421..5426}
	do
		addstring=$addstring$space$DKH$i.root
#		echo $addstring
	done
hadd BigKuramaHodoscope.root $addstring
