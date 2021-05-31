#!
for s in */
do
       cd $s

        for m in */
        do
                cd $m
                echo $m
                mv parmfile.txt parmfile
                cd ..
        done
       cd ..
done