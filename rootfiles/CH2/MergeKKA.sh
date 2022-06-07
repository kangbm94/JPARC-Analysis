DKH="DstKKAna0"
addstring=""
space=" "
for i in {5641..5666}
	do
		addstring=$addstring$space$DKH$i.root
#		echo $addstring
	done
hadd AllKKAna.root $addstring
