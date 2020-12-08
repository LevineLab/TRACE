%read in chlamy trait data from Sinead Collins 2017 experiment located in
%file mean_evopheno_slowgrowth_chlamy_selExp_2017_RL_forNate. Read the master_README.txt
%and selection_expt_description.docx files for description of data and experiments 
fformat='%f%s%f%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s';
headerlines=1;
delimiter=',';
filename='mean_evopheno_slowgrowth_chlamy_selExp_2017_RL_forNate.csv';
fid=fopen(filename);
TMP=textscan(fid, fformat, 'headerlines', headerlines, 'Delimiter', delimiter);
fclose(fid);

%create cell arrays representing traits of interest
all_traits = {'growth'; 'chla'; 'net-photo-cell'; 'dark.resp.cell';'cell.size';'daughter';'ROS.tol';'stress.recov'};
all_traits_num = {'1', 'growth';'2', 'chla';'3', 'net.photo.cell';'4', 'dark.resp.cell';'5', 'cell.size';'6', 'daughter';'7', 'ROS.tol';'8', 'stress.recov'};
all_traits_vec = [1 2 3 4 5 6 7 8];

%gather all of the numeric trait values across conditions and save to
%data_all table
for i=6:19
    tmp_vec = abs(cellfun(@str2double,TMP{i}));
    data_all(:,i-5)=tmp_vec;
end

%only select traits of interest
data_select = data_all(:,[1:4 11:14]); 

all_ids = TMP{2}; %gather all replicate ids
all_strains = TMP{3}; %isolate strain names
all_select = TMP{4}; %isolate selection co2
all_co2 = TMP{5}; %isolate assay co2

%create strings to match against
selection = 'ambient';
co2 = 'ambient';

%clear any data structures in the current work environment
clear tmp_strain_mat_ancestral
clear tmp_strain_mat_evolved
clear current_ids_ancestral
clear current_ids_evolved
clear current_strains_ancestral
clear current_strains_evolved

%Loop through all strains to create cell arrays containing strains and ids
for i=1:length(all_strains)

    select_tf = strcmp(all_select{i}, selection); %match selection regime
    co2_tf = strcmp(all_co2{i}, co2); %match assay regime

    if (select_tf == 1) && (co2_tf == 1) %ambient/ancestral strains

        tmp_strain_mat_ancestral(i,:) = data_select(i,:); %ancestral trait data

        current_ids_ancestral{i,1} = cellstr(all_ids{i}); %save ids 
        current_strains_ancestral{i,1} = all_strains(i); %save strains

    elseif (select_tf == 0) && (co2_tf == 0) %evolved strains

        tmp_strain_mat_evolved(i,:) = data_select(i,:); %evolved trait data 

        current_ids_evolved{i,1} = cellstr(all_ids{i}); %save ids
        current_strains_evolved{i,1} = all_strains(i); %save strains

    end



end

%remove any rows with zeros and create a matrix with ancestral trait values
indices = find(tmp_strain_mat_ancestral(:,1)==0);
tmp_strain_mat_ancestral(indices,:) = [];
ancestral_traits = tmp_strain_mat_ancestral; %final ancestral trait matrix

%remove any rows with zeros and create a matrix with evolved trait values
indices = find(tmp_strain_mat_evolved(:,1)==0);
tmp_strain_mat_evolved(indices,:) = [];
evolved_traits = tmp_strain_mat_evolved; %final evolved trait matrix

%remove any empty cells in cell arrays
current_ids_ancestral = current_ids_ancestral(~cellfun('isempty',current_ids_ancestral));
current_strains_ancestral = current_strains_ancestral(~cellfun('isempty',current_strains_ancestral));

current_ids_evolved = current_ids_evolved(~cellfun('isempty',current_ids_evolved));
current_strains_evolved = current_strains_evolved(~cellfun('isempty',current_strains_evolved));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%START SELECT 4 TRAITS WHERE PC1 + PC2 > 70% AND PC2 > 30%
%For the model, we are going to manually select 4 traits and their correlations
%according to the following criteria. In this code block, we will perform a PCA analysis 
%with all possible trait combinations of four. 
%We want the first 2 principal components to explain
%>70% of the variance with PC2 explaining >30%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%standardize trait data
ancestral_z=colstd(ancestral_traits);
evolved_z=colstd(evolved_traits);
dlmwrite('evolved_z.txt',evolved_z);

%save original ancestral and evolved standardized matrices
ancestral_z_og = ancestral_z; %make a copy of the standardized ancestral trait matrix
ancestral_z_og(19,:) = []; %the ancestral trait matrix has one more row than the evolved trait matrix. So, remove one row to match the evolved trait matrix.
evolved_z_og = evolved_z; %make a copy of the evolved trait matrix


all_trait_coms = nchoosek(all_traits_vec,4); %generate all possible combinations of 4 traits from 8 traits total
all_coms = length(all_trait_coms); %number of rows in all_trait_coms
selected_combs = [];
selected_combs_PoV_a = [];
selected_combs_PoV_e = [];
counter = 0;
for i=1:all_coms
    
    %take the ith row of traits
    current_traits = all_trait_coms(i,:);
    
    %extract the selected traits from the ancestral trait matrix
    current_trait_mat_a = ancestral_z(:,[current_traits(1) current_traits(2) current_traits(3) current_traits(4)]); 
    
    %extract the selected traits from the evolved matrix
    current_trait_mat_e = evolved_z(:,[current_traits(1) current_traits(2) current_traits(3) current_traits(4)]); 
    
    %%%%%%%Ancestral PCA of traits%%%%%%%%%%%
    [N,M]=size(current_trait_mat_a);
    R=current_trait_mat_a'*current_trait_mat_a/(length(current_trait_mat_a)-1);   % same as cov(X)
    [ V,Lambda ] = eigsort( R ); %V contains eigenvectors and Lambda eigenvalues
    %the total variance is the trace (sum of the diagonal elements)
    PoV_a=100*diag(Lambda)/trace(Lambda); % percent of variance
    Ar_a=V*Lambda.^0.5;   % Factor loading // each column corresponds to the factor loadings of a princicpal component on each of the original variables
    Sr_a=current_trait_mat_a*Ar_a;        % Ancestral Factor Score
    %%%%%%%Ancestral PCA of traits%%%%%%%%%%%
    
    
    %%%%%%%Evolved PCA of traits%%%%%%%%%%%
    [N,M]=size(current_trait_mat_e);
    R=current_trait_mat_e'*current_trait_mat_e/(length(current_trait_mat_e)-1);   % same as cov(X)
    [ V,Lambda ] = eigsort( R ); %V contains eigenvectors and Lambda eigenvalues
    %the total variance is the trace (sum of the diagonal elements)
    PoV_e=100*diag(Lambda)/trace(Lambda); % percent of variance
    Ar_e=V*Lambda.^0.5;   % Factor loading // each column corresponds to the factor loadings of a princicpal component on each of the original variables
    Sr_e=current_trait_mat_e*Ar_e;        % Ancestral Factor Score

    %%%%%%%Evolved PCA of traits%%%%%%%%%%%
    
    if ((PoV_a(1)+PoV_a(2)) >= 70 && PoV_a(2) >= 30) %do ancestral principal components meet the PC criteria?
        if ((PoV_e(1)+PoV_e(2)) >= 70 && PoV_e(2) >= 30) %do evolved principal components meet the PC criteria?
            counter = counter + 1;
            
            tmp_selected_combs_PoV_a = [PoV_a(1) PoV_a(2) PoV_a(3) PoV_a(4)]; %PoV for ancestral
            tmp_selected_combs_PoV_e = [PoV_e(1) PoV_e(2) PoV_e(3) PoV_e(4)]; %PoV for evolved
            
            %save trait combinations that pass criteria
            selected_combs = [selected_combs; current_traits]; %selected traits
            selected_combs_PoV_a = [selected_combs_PoV_a; tmp_selected_combs_PoV_a]; %selected traits PoV ancestral
            selected_combs_PoV_e = [selected_combs_PoV_e; tmp_selected_combs_PoV_e]; %selected traits PoV evolved
            
        end
    end
    

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END SELECT 4 TRAITS WHERE PC1 + PC2 > 70% AND PC2 > 20%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%START MANUALLY CHECKED PCA'S AND SELECTED TRAIT COMBINATION 
%Now we have the trait combinations that passed the above criteria. The
%traits and the numbers that correspond to them can be found in the cell array
%"all_traits_num". Now, we put our biology hats on and choose which traits
%we would like to insert into the model to represent our trait 
%space. We decided to choose the trait combination containing growth
%(trait 1), respiration (trait 4), cell size (trait 5), and daughter cell
%production (trait 6). This trait combination is in row 6 (i.e.
%selected_combs(6,:) -- 1, 4, 5, 6). Below, we select these traits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
manually_selected_traits = selected_combs(6,:);
ancestral_traits = ancestral_traits(:,selected_combs(6,:)); %pull selected ancestral traits
evolved_traits = evolved_traits(:,selected_combs(6,:)); %pull selected evolved traits
all_traits_num = all_traits_num(selected_combs(6,:),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END MANUALLY CHECKED PCA'S AND SELECTED TRAIT COMBINATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIALIZE VALUES FOR MODEL VARIABLES  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_mod = 1000; %population size
time = 2000; %number of generations

%The simulated population will travel towards a single evolutionary point
%in PCA trait space. Here, we specify the index (row) of the evolved individual
%we will be using as the evolutionary endpoint. In this case, we chose 
%individual 20 from evolved_z dataset (see L 995 in trait_evo_driver_v6_hpcc_start_clone_reps.m) 
%Each row in this file represents one evolved individual
%and their evolved trait values. This population is then multiplied by the
%evolved factor loadings (Ar_e_dat; See L 37 in trait_evo_driver_v6_hpcc_start_clone_reps.m)
%to calculate each evolved individual's
%coordinates in the evolved PCA trait space (saved in Sr_e). Then, the
%first individual specified by the variable sim_clone_num (L 995 in trait_evo_driver_v6_hpcc_start_clone_reps.m)
%is selected to be the evolutionary
%endpoint for which the starting simulated ancestral population is
%traveling towards.
sim_e_index = 20;  %specifies which individual to select as the evolutionary endpoint in the evolved population

%trait changes are randomly drawn from a normal distribution using function
%randn. tgrad specifies the standard deviation of the normal distribution.
%It is a vector to test a range of standard deviations
tgrad = [0.05]; %see L 475 in trait_evo_driver_v6_hpcc_start_clone_reps.m

%Similarly, correlation changes are randomly drawn from a normal distribution.
%cgrad specifies the standard deviation of the normal distribution.
%It is a vector to test a range of standard deviations
cgrad = [0.05];%see L 701 in trait_evo_driver_v6_hpcc_start_clone_reps.m

%nu_trait is the number of individuals to experience a random trait change
nu_trait = N_mod*(.9);
%nu_trait2 is the number of individuals to experience both a random trait
%change and correlation change
nu_trait2 = N_mod*(.1);

%clonal is a flag that forces the initial population to retain all of the
%same trait values. The trait values in the first row in ancestral_starting_pop.txt are used
%to generate the clonal population
%clonal = "TRUE";
clonal = myclone;

%mode is a flag that denotes whether an individual in a population is randomly sampled to
%the next generation weighted by fitness or not. If mode = "random", then
%individuals are randomly sampled without weighting. If mode is any other
%string, then individuals are randomly sampled to the next generation
%weighted by fitness.

%mode = "negative";
mode = mymode; 

%corr_all is a flag and if corr_all = "TRUE", then any of the correlations
%can be randomly selected to be changed. Always leave as "TRUE" for now.
%corr_all = "TRUE";
corr_all = mycorr; %if true, select any correlation to change; if false, only select significant correlations to change

%patch is a flag where if patch = 'TRUE', then the model can fluctuate
%between two different selection regimes. If patch = 'FALSE', then the
%model only runs in 1 selection mode. Always leave 'FALSE' for now.
%patch = 'FALSE';
patch = mypatch; %If FALSE, only sample with weighting


%If the variable patch = 'FALSE' and mode ~= 'random', then tfact sets the
%selection differential value when the population enters the 2nd selection
%regime in the model. It is not used for now.
tfact = mytfact; %no more than tfact (%) selection differential
%tfact = 0.1;

%select is a flag. If select = 0, then individuals in the population are
%randomly resampled not weighted by fitness. Always leave select = 1.
select = 1;

%If trait_flag = 1, then allow trait values to change. If trait_flag = 0,
%then trait values cannot change.
trait_flag = mytflag;
%trait_flag = 1;

%if corr_flag = 1, allow correlation values to change If corr_flag = 0,
%then correlation values cannot change.
corr_flag = mycflag;
%corr_flag = 1;

%runnumber keeps track of the run number of this model. We usually run 100
%replicate runs of this model. 
runnumber = myrun;
%runnumber = 1;


%The evolved population stored in evolved_traits only contains 19
%individuals. To match this number, remove 1 row from the ancestral
%population saved in ancestral_traits.
ancestral_traits(19,:) = []; 
current_strains_ancestral{19} = [];
current_strains_ancestral = current_strains_ancestral(~cellfun('isempty',current_strains_ancestral));

%Driver is a vector whose rows equal the number of generations saved in
%variable "time". 

dlmwrite('ancestral_traits.txt',ancestral_traits, 'delimiter', '\t');
dlmwrite('evolved_traits.txt',evolved_traits, 'delimiter', '\t');

Driver=ones(time,1);
%for q=t_split:t_split*2:time
%    Driver(q+1:q+t_split-1)=0; %generate time intervals to switch back and forth from
%end

%trait_evo_driver_v5_hpcc_start_clone_reps is the driver mfile that runs
%the model. Input are the variables generated above and output is a range
%of matrices.
[anc_sig_corrs_cell, pop_track_nu_all,track_ave_cor_changes_all,anc_sig_corrs_cell_all,Sr_ae_sim_cell_all,Sr_aa_sim_cell_all,...
    trait_a_cell, trait_b_cell, trait_c_cell, trait_d_cell, pop_cell_all] = trait_evo_driver_v6_hpcc_start_clone_reps(sim_e_index,N_mod,time,Driver, ancestral_traits, evolved_traits,...
all_traits_num, ancestral_z_og, evolved_z_og, all_traits, current_strains_evolved, manually_selected_traits,...
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

