SOT="SdcOutTracking0"
SIT="SdcInTracking0"
BOT="BcOutTracking0"
sotstring=""
sitstring=""
botstring=""
space=" "
for i in {5641..5646}
	do
		sotstring=$sotstring$space$SOT$i.root
		sitstring=$sitstring$space$SIT$i.root
		botstring=$botstring$space$BOT$i.root
#		echo $addstring
	done
hadd BigSdcOut.root $sotstring
hadd BigSdcIn.root $sitstring
hadd BigBcOut.root $botstring
