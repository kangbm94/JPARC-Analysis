HODO="Hodoscope0"
EA0C="Easiroc0"
BcOut="BcOutTracking0"
SdcIn="SdcInTracking0"
SdcOut="SdcOUtTracking0"
hodostring=""
ea0cstring=""
BcOutstring=""
SdcInstring=""
SdcOutstring=""
space=" "
for i in {5641..5646}
	do
		hodostring=$hodostring$space$HODO$i.root
		ea0cstring=$ea0cstring$space$EA0C$i.root
		BcOutstring=$BcOutstring$space$BcOut$i.root
		SdcInstring=$SdcInstring$space$SdcIn$i.root
		SdcOutstring=$SdcOutstring$space$SdcOut$i.root
#		echo $addstring
	done
hadd AllHodoscope.root $hodostring
hadd AllEasiroc.root $ea0cstring
hadd AllBcOutTracking.root $BcOutstring
hadd AllSdcInTracking.root $SdcInstring
hadd AllSdcOutTracking.root $SdcOutstring
