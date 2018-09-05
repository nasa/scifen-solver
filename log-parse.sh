#!/bin/bash
#
if [ "$1" != "" ]; then
   BASENAME=$1
else
   BASENAME=time
fi
#
HOSTNAME=$(hostname)
#
rm -f $HOSTNAME/*.tmp.csv
for RTOL in "1e-3" "1e-5" "1e-4" ; do
   for KSP in \
   symmlq cr    cgs    lcd    cg     fcg \
   bicg   bcgs  ibcgs  bcgsl             \
   minres gmres lgmres dgmres tfqmr ;
   do
   sed -n -e 's/^Time (sec):          //p' \
       $HOSTNAME/$BASENAME-ksp-$KSP-r$RTOL.log 2>/dev/null \
     | cut -c 1-10 >> $HOSTNAME/time.tmp.csv
   grep -s 'jacobi,' $HOSTNAME/$BASENAME-ksp-$KSP-r$RTOL.log \
      >> $HOSTNAME/data.tmp.csv
   done
done
#
paste -d ',' $HOSTNAME/data.tmp.csv $HOSTNAME/time.tmp.csv \
   > $HOSTNAME-$BASENAME.csv
#
