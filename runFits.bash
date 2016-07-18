rawFile=$1
correctedFile=$2
fitFile=$3
resolution=$4
q2OrMult=$5
q2MultBin=$6
zdcBin=$7

root -b -q 'doCorrectHistogram.C("'$rawFile'","'$correctedFile'",'$resolution')'
# for((i=0; i <= 0; i++))
# for((i=0; i <= 7; i++))
# do
#     # root -b -q 'fitManager.C("'$correctedFile'","'$fitFile'", '$zdcBin', 0, '$q2OrMult', '$q2MultBin', '$i')'
#     root -b -q 'fitManager.C("'$rawFile'","'$fitFile'", '$zdcBin', 0, '$q2OrMult', '$q2MultBin', '$i')'
# done
