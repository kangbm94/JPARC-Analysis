sleep 2000

rsync	-av -P kek:/home/had/kangbm/k18-analyzer_2/rootfiles/Callibration/Bc*.root .
rsync	-av -P kek:/home/had/kangbm/k18-analyzer_2/rootfiles/Callibration/SdcIn*.root .
rsync	-av -P kek:/home/had/kangbm/k18-analyzer_2/rootfiles/Callibration/SdcOut*.root .
rsync	-av -P kek:/home/had/kangbm/k18-analyzer_2/rootfiles/Callibration/Easiroc*.root .
rsync	-av -P kek:/home/had/kangbm/k18-analyzer_2/rootfiles/Callibration/Hodo*.root .
