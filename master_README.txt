##########################################
## How to run the model ##
##########################################


This folder contains the model files and all the subroutines. Run the m-file "trait_evo_sim_v6_hpcc_start_clone_reps.m" to set up the inputs to the actual model contained in "trait_evo_driver_v6_hpcc_start_clone_reps.m". "trait_evo_driver_v6_hpcc_start_clone_reps.m" is the large model that calls all the other subroutine functions that are also contained within that folder. 

We have generated random ancestral and evolved trait values that you can load in to run the model. They are in the text files ?ancestral_traits.txt? and ?evolved_traits.txt?, respectively. These are supposed to be in place of empirical trait values that you can feed into the model. For this run, there are 4 traits across 29 ?empirical? individuals that were ?measured? in both the ancestral and evolved environments. In the model in the paper, we chose the traits, 1) growth 2) respiration 3) cell size and 4) daughter cell production. So, you will see those as placeholders here. 

1) The best way to get familiar is to first open "trait_evo_sim_v6_hpcc_start_clone_reps.m" and manually run lines of code (or chunks of code) so you can see how things are starting to set up.

2) Eventually, when you near the bottom of "trait_evo_sim_v6_hpcc_start_clone_reps.m", you will see that ?trait_evo_driver_v6_hpcc_start_clone_reps.m? is called once all the inputs are set up.

3) Then open up "trait_evo_driver_v6_hpcc_start_clone_reps.m" and run lines manually to see how the model runs with the subroutines.

4) Results are eventually output back to "trait_evo_sim_v6_hpcc_start_clone_reps.m" and saved within the same folder.