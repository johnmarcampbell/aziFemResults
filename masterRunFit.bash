auauOrUu=0
q2OrMult=1
zdcBin=1

if [[ $auauOrUu -eq 1 ]]; then
    species=UUFemto193
    dataDirectory=193GeVData
else
    species=AuAuFemto200
    dataDirectory=200GeVData
fi

if [[ $q2OrMult -eq 1 ]]; then
    tag=mult
else
    tag=q2
fi

if [[ $zdcBin -eq 1 ]]; then
    binDirectory=0to0x5PercZdc
else
    binDirectory=0to0x25PercZdc
fi
echo $species 
echo $tag
echo $binDirectory

# declare -a res=(0.260516 0.26218 0.256748 0.270716 0.287465) # Au - 0.0% to 0.25% - q2
# declare -a res=(0.336965 0.331986 0.338615 0.354593 0.367271) # U - 0.0% to 0.25% - q2
# declare -a res=(0.341097 0.328521 0.326276 0.326144 0.332257) # Au - 0.0% to 0.25% - mult
# declare -a res=(0.425452 0.424873 0.420491 0.41498 0.406762) # U - 0.0% to 0.25% - mult
# declare -a res=(0.261268 0.265055 0.264835 0.272727 0.282356) # Au - 0.0% to 0.5% - q2
# declare -a res=(0.336164 0.336912 0.343711 0.349486 0.372916) # U - 0.0% to 0.5% - q2
declare -a res=(0.338668 0.330646 0.330272 0.327157 0.329236) # Au - 0.0% to 0.5% - mult
# declare -a res=(0.426443 0.424254 0.422858 0.420137 0.41072) # U - 0.0% to 0.5% - mult
# declare -a res=(1 1 1 1 1) # Dummy res 

for((i=2; i <= 2; i++))
# for((i=0; i <= 4; i++))
do
    inFile=$dataDirectory/$binDirectory/${species}_${tag}_$i.root
    correctedFile=$dataDirectory/$binDirectory/${species}_${tag}_${i}_CorrectedTest.root
    # fitFile=$dataDirectory/$binDirectory/${species}_${tag}_${i}_Fit.root
    fitFile=$dataDirectory/$binDirectory/${species}_${tag}_${i}_FitNoResTest.root
    logFile=$dataDirectory/$binDirectory/${species}_${tag}_${i}.log

    ( ./runFits.bash $inFile $correctedFile $fitFile ${res[$i]} $q2OrMult $i $zdcBin >& $logFile & )

done
