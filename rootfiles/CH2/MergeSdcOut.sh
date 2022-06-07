SOT="SdcOutTracking0"
addstring=""
space=" "
for i in {5641..5646}
	do
		addstring=$addstring$space$SOT$i.root
#		echo $addstring
	done
hadd BigSdcOut.root $addstring
