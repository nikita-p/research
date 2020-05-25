JOBS_COUNTER=0
MAX_CHILDREN=10
MY_PID=$$
KEY="19"

source KinScript.sh

>log.txt
while read LINE
do
    echo Washing: $LINE
    JOBS_COUNTER=$((`ps ax -Ao ppid | grep $MY_PID | wc -l`))
    while [ $JOBS_COUNTER -ge $MAX_CHILDREN ]
    do
        JOBS_COUNTER=$((`ps ax -Ao ppid | grep $MY_PID | wc -l`))
        echo Jobs counter: $JOBS_COUNTER
        sleep 10
    done
    root -l -b "macro2.cpp(\"$LINE\", \"$KEY\")" &
done < /spoolA/petrov/research/inputs/$KEY/trees
wait
