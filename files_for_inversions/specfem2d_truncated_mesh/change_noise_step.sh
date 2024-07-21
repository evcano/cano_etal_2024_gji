# noise simulation step
step=$1


case $step in
1) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 1:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" DATA/Par_file
   ;;
2) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 2:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .true.:" DATA/Par_file
   ;;
3) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 3:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 3:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" DATA/Par_file
   ;;
*) echo "step not recognized: $step"; echo "please use as step number 1, 2 or 3"; exit 1 ;;
esac
echo
