
bkpdir=new_test/.bkp$(expr 0 + $(ls -d new_test/.bkp* 2>/dev/null| wc -l))
mkdir $bkpdir

mv new_test/anav2_* $bkpdir

mv anav2_pi0pipm_dpm_ip0_brem.root new_test/anav2_pip_pim_brem_plab5.5.root
mv anav2_pi0pipm_dpm_ip1_brem.root new_test/anav2_pip_pim_brem_plab8.0.root
mv anav2_pi0pipm_dpm_ip2_brem.root new_test/anav2_pip_pim_brem_plab12.0.root

mv anav2_pi0jpsi_tda_ip0_brem.root new_test/anav2_jpsi_brem_plab5.5.root
mv anav2_pi0jpsi_tda_ip1_brem.root new_test/anav2_jpsi_brem_plab8.0.root
mv anav2_pi0jpsi_tda_ip2_brem.root new_test/anav2_jpsi_brem_plab12.0.root
