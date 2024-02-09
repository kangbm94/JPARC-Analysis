#sleep 1000
#rsync	-avh -P kek:/home/had/kangbm/ws/k18-analyzer/rootfiles/CH2/run056*_DstTPC*.root .
rsync	-avh -P kek:/home/had/kangbm/ws/k18-analyzer/rootfiles/CH2/run056*_DstKScat*.root ./CH2
rsync	-avh -P kek:/home/had/kangbm/ws/k18-analyzer/rootfiles/Prod/run056*_DstKScat*.root ./Prod
#rsync	-avh -P kek:/home/had/kangbm/ws/e42_k18-analyzer_ws/rootfiles/CH2/run056*Genfit*.root .
