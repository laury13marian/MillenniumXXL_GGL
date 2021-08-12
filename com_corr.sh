#!/bin/csh

echo COMPUTE THE CORRELATION FUNCTION FROM PAIR COUNTS

set NSUB = 215

set NTASK = 31

#set iTASK = 0

#set iSC = 0

#while ( $iSC <= $NSUB )

set iTASK = 0

while ( $iTASK <= $NTASK )

#echo COMPUTING CORRELATIONS FOR TASK $iTASK AND SUBCUBE $iSC
#echo COMPUTE DELTA_SIGMA FOR $iTASK AND SUBCUBE $iSC
#./gds $iSC 100 
./stat $iTASK
#./weight $iSC 100 $iTASK
#./cot $iSC $iTASK

@ iTASK++

end

#@ iSC++

#end


