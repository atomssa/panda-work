source /vol0/panda/svn/pandaroot-trunk-2014.02.21/build/config.sh

echo -n @@@ Step 1 Start:"  "  ; date 
root -b -q sim_complete.C"(10000)"
echo -n @@@ Step 1   End:"  "  ; date

echo -n @@@ Step 2 Start:"  "  ; date
root -b -q digi_complete.C
echo -n @@@ Step 2   End:"  "  ; date

echo -n @@@ Step 3 Start:"  "  ; date
root -b -q reco_complete.C
echo -n @@@ Step 3   End:"  "  ; date

echo -n @@@ Step 4 Start:"  "  ; date
root -b -q pid_complete.C"(0)"
echo -n @@@ Step 4   End:"  "  ; date
