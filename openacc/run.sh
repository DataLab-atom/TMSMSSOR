for numworks in 1 2 4 8 16 32 64 128 256 512 ;  
do  
    for m in 256 512 1024 2048 4096;  
    do
    ./main_acc $numworks $m 1 4 & ./main_acc $numworks $m 2 5 ;
    done
done  