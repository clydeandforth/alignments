# change contig header from spades
awk '/NODE/{print ">Contig"(++n);next} { print }' 
