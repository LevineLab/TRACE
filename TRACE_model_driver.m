%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIALIZE VALUES FOR MODEL VARIABLES  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ancestral_traits = dlmread('ancestral_traits.txt');
evolved_traits = dlmread('evolved_traits.txt');
all_traits_num = {'1', 'growth';'4', 'dark.resp.cell';'5', 'cell.size';'6', 'daughter'};
manually_selected_traits = [1,4,5,6];

%input_data = {ancestral_traits, evolved_traits, all_traits_num, manually_selected_traits};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIALIZE VALUES FOR MODEL VARIABLES  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_mod = 1000; %population size
time = 2; %number of generations

%The simulated population will travel towards a single evolutionary point
%in PCA trait space. Here, we specify the index (row) of the evolved individual
%we will be using as the evolutionary endpoint. In this case, we chose 
%individual 20 from evolved_z dataset (see L 995 in TRACE_model.m) 
%Each row in this file represents one evolved individual
%and their evolved trait values. This population is then multiplied by the
%evolved factor loadings (Ar_e_dat; See L 37 in TRACE_model.m)
%to calculate each evolved individual's
%coordinates in the evolved PCA trait space (saved in Sr_e). Then, the
%first individual specified by the variable sim_clone_num (L 995 in TRACE_model.m)
%is selected to be the evolutionary
%endpoint for which the starting simulated ancestral population is
%traveling towards.
sim_e_index = 20;  %specifies which individual to select as the evolutionary endpoint in the evolved population

%trait changes are randomly drawn from a normal distribution using function
%randn. tgrad specifies the standard deviation of the normal distribution.
%It is a vector to test a range of standard deviations
tgrad = [0.05]; %see L 475 in TRACE_model.m

%Similarly, correlation changes are randomly drawn from a normal distribution.
%cgrad specifies the standard deviation of the normal distribution.
%It is a vector to test a range of standard deviations
cgrad = [0.05];%see L 701 in TRACE_model.m

%nu_trait is the number of individuals to experience a random trait change
nu_trait = N_mod*(.9);
%nu_trait2 is the number of individuals to experience both a random trait
%change and correlation change
nu_trait2 = N_mod*(.1);

%clonal is a flag that forces the initial population to retain all of the
%same trait values. The trait values in the first row in ancestral_starting_pop.txt are used
%to generate the clonal population
clonal = "TRUE";
%clonal = myclone;

%mode is a flag that denotes whether an individual in a population is randomly sampled to
%the next generation weighted by fitness or not. If mode = "random", then
%individuals are randomly sampled without weighting. If mode is any other
%string, then individuals are randomly sampled to the next generation
%weighted by fitness.

mode = "negative";
%mode = mymode; 

%corr_all is a flag and if corr_all = "TRUE", then any of the correlations
%can be randomly selected to be changed. Always leave as "TRUE" for now.
corr_all = "TRUE";
%corr_all = mycorr; %if true, select any correlation to change; if false, only select significant correlations to change

%patch is a flag where if patch = 'TRUE', then the model can fluctuate
%between two different selection regimes. If patch = 'FALSE', then the
%model only runs in 1 selection mode. Always leave 'FALSE' for now.
patch = 'FALSE';
%patch = mypatch; %If FALSE, only sample with weighting


%If the variable patch = 'FALSE' and mode ~= 'random', then tfact sets the
%selection differential value when the population enters the 2nd selection
%regime in the model. It is not used for now.
%tfact = mytfact; %no more than tfact (%) selection differential
tfact = 0.1;

%select is a flag. If select = 0, then individuals in the population are
%randomly resampled not weighted by fitness. Always leave select = 1.
select = 1;

%If trait_flag = 1, then allow trait values to change. If trait_flag = 0,
%then trait values cannot change.
%trait_flag = mytflag;
trait_flag = 1;

%if corr_flag = 1, allow correlation values to change If corr_flag = 0,
%then correlation values cannot change.
%corr_flag = mycflag;
corr_flag = 1;

%runnumber keeps track of the run number of this model. We usually run 100
%replicate runs of this model. 
%runnumber = myrun;
runnumber = 1;

%Driver is a vector whose rows equal the number of generations saved in
%variable "time". 

Driver=ones(time,1);
%for q=t_split:t_split*2:time
%    Driver(q+1:q+t_split-1)=0; %generate time intervals to switch back and forth from
%end

[anc_sig_corrs_cell, pop_track_nu_all,track_ave_cor_changes_all,anc_sig_corrs_cell_all,Sr_ae_sim_cell_all,Sr_aa_sim_cell_all,...
    trait_a_cell, trait_b_cell, trait_c_cell, trait_d_cell, pop_cell_all] = TRACE_model(sim_e_index,N_mod,time,Driver, ancestral_traits, evolved_traits,...
all_traits_num, manually_selected_traits,...
nu_trait, nu_trait2, clonal, corr_all, patch, tfact, mode, tgrad, cgrad, select, trait_flag, corr_flag, runnumber);


%save data structures
save('pop_track_nu_all.mat','pop_track_nu_all');
save('track_ave_cor_changes_all.mat','track_ave_cor_changes_all');
save('anc_sig_corrs_cell_all.mat','anc_sig_corrs_cell_all');
save('anc_sig_corrs_cell.mat','anc_sig_corrs_cell');
save('Sr_ae_sim_cell_all.mat','Sr_ae_sim_cell_all');
save('Sr_aa_sim_cell_all.mat','Sr_aa_sim_cell_all');
save('trait_a_cell.mat','trait_a_cell');
save('trait_b_cell.mat','trait_b_cell');
save('trait_c_cell.mat','trait_c_cell');
save('trait_d_cell.mat','trait_d_cell');
save('pop_cell_all.mat','pop_cell_all');

