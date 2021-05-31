#!
for s in */
do
       cd $s

        for m in */
        do
                cd $m
                echo $m
                screen -d -m ~/migrate-n
                cd ..
        done
       cd ..
done