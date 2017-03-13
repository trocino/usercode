#!/bin/csh

setenv locdir `pwd` 
setenv afsdir AFSTEMPLDIR

set pytfile=PYTTEMPLFILE
set logfile=LOGTEMPLFILE

cd ${afsdir} 
eval `scramv1 runtime -csh` 
cp ${pytfile} ${locdir}
cd -
cmsRun ${pytfile} >&! ${logfile}
\ls -lrt
foreach rootfile ( `\ls *.root` ) 
    cp ${rootfile} ${afsdir} 
end 
foreach outlogfile ( `\ls *.log` ) 
    cp ${outlogfile} ${afsdir} 
end 

exit 
