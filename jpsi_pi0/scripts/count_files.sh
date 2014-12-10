for t in jpsi pip_pim;
do
    for p in 5.5 8.0 12.0;
    do
	listf=lists/pid_${t}_${p}.list
	echo $listf
	rm $listf
	touch $listf
	for f in `ls output/pid/pbar_p_${t}_pi0_plab${p}_*.root | sort`;
	do
	    fid=`echo $f| sed s/.root//| awk -Fplab '{print $2}' | awk -F_ '{print $2}'`
	    echo $fid $f >> $listf;
	done;
    done
done
