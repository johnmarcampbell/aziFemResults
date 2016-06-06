#declare -a resMult=(0.347914 0.33366 0.331359 0.32536 0.330831) # AuAu 0 to 1% Mult
#declare -a resq2=(0.262934 0.267016 0.259908 0.267436 0.281443) # AuAu 0 to 1% q2
# AuAu 0% to 0.5% - Mult
# declare -a resMult=(0.341188 0.332383 0.329073 0.326508 0.330878)
# AuAu 0.5% to 1% - Mult
# declare -a resMult=(0.348646 0.3345 0.333354 0.324128 0.330758)
# AuAu 0% to 0.5% - q2
#declare -a resq2=(0.25436 0.26111 0.25547 0.263246 0.274286)
# AuAu 0.5% to 1% - q2
# declare -a resq2=(0.269405 0.271206 0.262421 0.269567 0.286708)

# UU 0% to 0.5% - Mult
# declare -a resMult=(0.425093 0.42712 0.42153 0.415694 0.406869)
# UU 0.5% to 1% - Mult
#declare -a resMult=(0.432553 0.425807 0.424784 0.422587 0.416744)
# UU 0% to 0.5% - q2
#declare -a resq2=(0.33675 0.332326 0.337774 0.353558 0.367908)
# UU 0.5% to 1% - q2
#declare -a resq2=(0.335093 0.341606 0.346129 0.35192 0.375451)

# declare -a resMult=(0.341097 0.328521 0.326276 0.326144 0.332257) # Au - 0.0% to 0.25% - mult
# declare -a resq2=(0.260516 0.26218 0.256748 0.270716 0.287465) # Au - 0.0% to 0.25% - q2
# declare -a resMult=(0.338668 0.330646 0.330272 0.327157 0.329236) # Au - 0.0% to 0.5% - mult
# declare -a resq2=(0.261268 0.265055 0.264835 0.272727 0.282356) # Au - 0.0% to 0.5% - q2
# declare -a resMult=(0.425452 0.424873 0.420491 0.41498 0.406762) # U - 0.0% to 0.25% - mult
# declare -a resq2=(0.336965 0.331986 0.338615 0.354593 0.367271) # U - 0.0% to 0.25% - q2
declare -a resMult=(0.426443 0.424254 0.422858 0.420137 0.41072) # U - 0.0% to 0.5% - mult
# declare -a resq2=(0.336164 0.336912 0.343711 0.349486 0.372916) # U - 0.0% to 0.5% - q2

# directory=200GeVData/0to0x25PercZdc
# directory=200GeVData/0to0x5PercZdc
# directory=193GeVData/0to0x25PercZdc
directory=193GeVData/0to0x5PercZdc

# tag=q2
# q2OrMult=0
tag=mult
q2OrMult=1

# for((i=3; i <= 3; i++))
for((i=0; i <= 4; i++))
do
    inFile=$directory/UUFemto193_${tag}_$i.root
    correctedFile=$directory/UUFemto193_${tag}_${i}_Corrected.root
    fitFile=$directory/UUFemto193_${tag}_${i}_Fit.root
    logFile=$directory/UUFemto193_${tag}_${i}.log

    # ( ./runFits.bash $inFile $correctedFile $fitFile ${resq2[$i]} $q2OrMult $i >& $logFile & )
    ( ./runFits.bash $inFile $correctedFile $fitFile ${resMult[$i]} $q2OrMult $i >& $logFile & )

done
