#! /bin/csh

setenv inputdir ${1} 
if ( "${inputdir}" == "" ) then 
    echo " *** ERROR: input EOS directory required! " 
    exit 1 
endif 

setenv itsnow `echo $inputdir | gawk -F"/" '{print $(NF)}'`
if ( "$itsnow" == "" ) then 
    setenv itsnow `echo $inputdir | gawk -F"/" '{print $(NF-1)}'`
endif 
setenv itsnow "optim_"${itsnow} 
mkdir ${itsnow} 
cd ${itsnow} 
setenv afsdir `pwd` 

set queue='2nd' 

set cfgfile='TEMPL_me0muonAnalyzer_optim_cfg.py'
set jobfile='TEMPL_singleJob.csh' 

set filesperjob=5 
set ijob=0 
set div=0 
set rem=0 

set filelist="" 

foreach ifile ( `eos ls ${inputdir} | \grep root` ) 
    @ ijob++ 
    if ( "${filelist}" != "" ) then 
	set filelist=${filelist}","
    endif 
    set filelist=${filelist}"'root://eoscms/"${inputdir}"/"${ifile}"'"  
    #echo $filelist
    @ div = $ijob / $filesperjob 
    @ rem = $ijob % $filesperjob 
    if ( "${rem}" == "0" ) then 
	set newcfgfile='me0muonAnalyzer_'${div}'_cfg.py' 
	set newjobfile='singleJob_'${div}'.csh' 
	set newoutfile='anMe0Muons_muGun_noYcuts_'${div}'.root' 
	set newlogfile='log_'${div}'.log' 
	set newtmpfile='tmp_'${div}'.tmp' 
	cat ../${cfgfile} | sed "s?TEMPLFILELIST?${filelist}?g;s?TEMPLOUTFILE?${newoutfile}?g" > ${newcfgfile} 
	cat ../${jobfile} | sed "s?AFSTEMPLDIR?${afsdir}?g;s/PYTTEMPLFILE/${newcfgfile}/g;s/LOGTEMPLFILE/${newlogfile}/g" > ${newjobfile} 
	chmod +x ${newjobfile} 
	bsub -q ${queue} -o ${newtmpfile} ${newjobfile} 
	set filelist=""
    endif 
end 

# One last time, if the last gorup of files contained less than ${filesperjob} files 
if ( "${rem}" != "0" ) then 
    @ div++ 
    set newcfgfile='me0muonAnalyzer_'${div}'_cfg.py' 
    set newjobfile='singleJob_'${div}'.csh' 
    set newoutfile='anMe0Muons_muGun_noYcuts_'${div}'.root' 
    set newlogfile='log_'${div}'.log' 
    set newtmpfile='tmp_'${div}'.tmp' 
    cat ../${cfgfile} | sed "s?TEMPLFILELIST?${filelist}?g;s?TEMPLOUTFILE?${newoutfile}?g" > ${newcfgfile} 
    cat ../${jobfile} | sed "s?AFSTEMPLDIR?${afsdir}?g;s/PYTTEMPLFILE/${newcfgfile}/g;s/LOGTEMPLFILE/${newlogfile}/g" > ${newjobfile} 
    chmod +x ${newjobfile} 
    bsub -q ${queue} -o ${newtmpfile} ${newjobfile} 
    set filelist=""
endif 

cd .. 
 
exit 
