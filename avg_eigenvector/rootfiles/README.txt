Respin 2 eigenvector cleanup analysis.
Obtained these outputs with the following scripts:


$ root -l -b -q SortEigenvectors_CREX_all.C'("mini_eigen_reg_5bpms")' > & rootfiles/5bpm_sorting.out
$ root -l -b -q SortEigenvectors_CREX_all.C'("mini_eigen_reg_allbpms")' > & rootfiles/allbpm_sorting.out
$ root -l -b -q merge_sorted_trees.C'("","")'

-rw-r--r-- 1 cameronc a-parity 206M May 30 20:47 rcdb_eigenvectors.root                                           Unsorted eigenvector regression outputs
-rw-r--r-- 1 cameronc a-parity 4.0M May 30 20:32 mini_eigen_reg_5bpms_sorted.root                                 Regular sort (no fancy cuts or anything) - 5bpm
-rw-r--r-- 1 cameronc a-parity 1.9K May 30 20:32 mini_eigen_reg_5bpms_skeleton_matrix.root                        Regular sort (no fancy cuts or anything) - 5bpm skeleton
-rw-r--r-- 1 cameronc a-parity 4.0M May 30 20:32 5bpm_sorting.out                                                 Regular sort (no fancy cuts or anything) - 5bpm output file
-rw-r--r-- 1 cameronc a-parity  15M May 30 20:33 mini_eigen_reg_allbpms_sorted.root                               Regular sort (no fancy cuts or anything) - allbpm
-rw-r--r-- 1 cameronc a-parity 3.7K May 30 20:33 mini_eigen_reg_allbpms_skeleton_matrix.root                      Regular sort (no fancy cuts or anything) - allbpm skeleton
-rw-r--r-- 1 cameronc a-parity  13M May 30 20:33 allbpm_sorting.out                                               Regular sort (no fancy cuts or anything) - allbpm output file

-rw-r--r-- 1 cameronc a-parity 205M May 30 20:48 rcdb_eigenvectors_sorted.root                                    Regular sort (no fancy cuts or anything) - collection of all trees (master file here)
                                                                                                                  Note... reg_asym_us_dd not included in the regular sort outputs...

$ root -l -b -q SortEigenvectors_CREX_stabilize.C'("mini_eigen_reg_5bpms","_sorted",5)' > & rootfiles/5bpm_sorting_round3.out
$ root -l -b -q SortEigenvectors_CREX_stabilize.C'("mini_eigen_reg_allbpms","_sorted",12)' > & rootfiles/allbpm_sorting_round3.out
$ root -l -b -q merge_sorted_trees.C'("","_sorted")'

-rw-r--r-- 1 cameronc a-parity 4.3M May 31 18:00 mini_eigen_reg_5bpms_sorted_sorted_parts.root                    Initial try to segment the parts out in double sorting
-rw-r--r-- 1 cameronc a-parity 4.0M May 31 18:00 5bpm_sorting_round2.out                                          output file of that try
-rw-r--r-- 1 cameronc a-parity 4.3M May 31 18:11 mini_eigen_reg_5bpms_sorted_sorted_reversed_parts.root           Initial try to reverse order do the averaging/ring mapping (not friendable)
-rw-r--r-- 1 cameronc a-parity 4.3M May 31 18:47 mini_eigen_reg_5bpms_sorted_sorted_stabilized_parts.root         Double sort (fancy cuts and such) - 5bpm - initial version (likely the same)
-rw-r--r-- 1 cameronc a-parity 4.3M May 31 20:12 mini_eigen_reg_5bpms_sorted_sorted.root                          Double sort (fancy cuts and such) - 5bpm 
-rw-r--r-- 1 cameronc a-parity 1.9K May 31 20:12 mini_eigen_reg_5bpms_sorted_sorted_skeleton_matrix.root          Double sort (fancy cuts and such) - 5bpm skeleton
-rw-r--r-- 1 cameronc a-parity 4.2M May 31 20:12 5bpm_sorting_round3.out                                          Double sort (fancy cuts and such) - 5bpm output file
-rw-r--r-- 1 cameronc a-parity  13M May 31 20:13 mini_eigen_reg_allbpms_sorted_sorted.root                        Double sort (fancy cuts and such) - allbpm 
-rw-r--r-- 1 cameronc a-parity 3.3K May 31 20:13 mini_eigen_reg_allbpms_sorted_sorted_skeleton_matrix.root        Double sort (fancy cuts and such) - allbpm skeleton
-rw-r--r-- 1 cameronc a-parity  12M May 31 20:13 allbpm_sorting_round3.out                                        Double sort (fancy cuts and such) - allbpm output file

-rw-r--r-- 1 cameronc a-parity 204M May 31 20:26 rcdb_eigenvectors_sorted_sorted.root                             Double sort (fancy cuts and such) - collection of all trees (master file here)
