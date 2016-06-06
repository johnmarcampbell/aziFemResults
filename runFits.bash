rawFile=$1
correctedFile=$2
fitFile=$3
resolution=$4
q2OrMult=$5
q2MultBin=$6

# root -b -q 'doCorrectHistogram.C("'$rawFile'","'$correctedFile'",'$resolution')'
for((i=3; i <= 3; i++))
# for((i=0; i <= 7; i++))
do
    root -b -q 'fitManager.C("'$correctedFile'","'$fitFile'", 0, 0, '$q2OrMult', '$q2MultBin', '$i')'
done
