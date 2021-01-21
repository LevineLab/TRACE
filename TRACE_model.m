function [anc_sig_corrs_cell, pop_track_nu_all,track_ave_cor_changes_all,anc_sig_corrs_cell_all,Sr_ae_sim_cell_all,Sr_aa_sim_cell_all,...
    trait_a_cell, trait_b_cell, trait_c_cell, trait_d_cell, pop_cell_all] = TRACE_model(sim_e_index,N_mod,time,Driver, ancestral_traits, evolved_traits,...
all_traits_num, manually_selected_traits,...
nu_trait, nu_trait2, clonal, corr_all, patch, tfact, mode, tgrad, cgrad, select, trait_flag, corr_flag, runnumber)

%%%%%%%Ancestral PCA of traits%%%%%%%%%%%
X=colstd(ancestral_traits); %standardize traits
ancestral_z = X; %make a copy of the standardized traits
[N,M]=size(X);
R=X'*X/(length(X)-1);   % same as cov(X)
[ V,Lambda ] = eigsort( R ); %V contains eigenvectors and Lambda eigenvalues
%the total variance is the trace (sum of the diagonal elements)
PoV_a_dat=100*diag(Lambda)/trace(Lambda); % percent of variance
Ar=V*Lambda.^0.5;   % Factor loading // each column corresponds to the factor loadings of a princicpal component on each of the original variables
%plot factor loadings vs the original variables
Ar_a_dat = Ar;

Sr_aa=X*Ar_a_dat;        % Ancestral Factor Score

dlmwrite('ancestral_z.txt',ancestral_z);
%%%%%%%End Ancestral PCA of traits%%%%%%%%%%%


%%%%%%%Evolved PCA of traits%%%%%%%%%%%
X=colstd(evolved_traits);
evolved_z = X;
[N,M]=size(X);
R=X'*X/(length(X)-1);   % same as cov(X)
[ V,Lambda ] = eigsort( R ); %V contains eigenvectors and Lambda eigenvalues
%the total variance is the trace (sum of the diagonal elements)
PoV_e_dat=100*diag(Lambda)/trace(Lambda); % percent of variance
Ar=V*Lambda.^0.5;   % Factor loading // each column corresponds to the factor loadings of a princicpal component on each of the original variables
Ar_e_dat = Ar;

dlmwrite('Ar_e_dat.txt',Ar_e_dat,'\t');

Sr_e=X*Ar_e_dat;        % Evolved Factor Score

dlmwrite('evolved_z.txt',evolved_z);
%%%%%%%End Evolved PCA of traits%%%%%%%%%%%


%%%%%%%Ancestral projected onto Evolved PCA of traits%%%%%%%%%%%
Y=ancestral_z;
Sr_ae=Y*Ar_e_dat; % Factor Score
%%%%%%%End Ancestral projected onto Evolved PCA of traits%%%%%%%%%%%


%calculate distance between ancestral and evolved points in the ancestral and evolved PC
%space

Sr_dist_a = [];
Sr_dist = [];
for i=1:length(Sr_ae)
    
    %tmp_a = [Sr_aa(i,1),Sr_aa(i,2);Sr_e(i,1),Sr_e(i,2)];
    %Sr_dist_a(i) = pdist(tmp_a,'euclidean'); %distances between ancestral and evolved points in ancestral PCA
    
    tmp = [Sr_ae(i,1),Sr_ae(i,2);Sr_e(i,1),Sr_e(i,2)];
    Sr_dist(i) = pdist(tmp,'euclidean'); %distances between ancestral and evolved points in evolved PCA
end


%%%%%%%%%%%%%%ANCESTRAL CORRELATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make cell array of all traits
%all_traits_num = {'t1', 'chla';'t2', 'cell.size';'t3', 'daughter';'t4', 'stress.recov'};
%all_traits_num = {'t1', 'growth';'t2', 'chla';'t3', 'net.photo.cell';'t4', 'dark.resp.cell';'t5', 'cell.size';'t6', 'daughter';'t7', 'ROS.tol';'t8', 'stress.recov'};
%all_traits_num = {'1', 'growth';'2', 'chla';'3', 'net.photo.cell';'4', 'dark.resp.cell';'5', 'cell.size';'6', 'daughter';'7', 'ROS.tol';'8', 'stress.recov'};
    
anc_corr_cell_sorted = create_anc_corr_mat(ancestral_z, all_traits_num);

%randomize
c=clock;
pause(c(6)*runnumber/100)

%%%%%%%%%%%%%%EVOLVED CORRELATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evo_corr_cell_sorted = create_evo_corr_mat(evolved_z,  all_traits_num);

%set ancestral correlations%%%%%%%%%%%%%%
%randomize
c=clock;
pause(c(6)*runnumber/100)

dlmwrite('anc_corr_cell_sorted.txt',anc_corr_cell_sorted,'\t');

dlmwrite('evo_corr_cell_sorted.txt',evo_corr_cell_sorted,'\t');


%randomize
c=clock;
pause(c(6)*runnumber/100)

grad_all = [1];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%MODEL BEGINS%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lincol = {'m';'c';'r';'g';'b';};
%pop_track_nu = zeros(time, 8, N_mod); %3D matrix to keep trait and correlation changes over time

pop_track_nu = cell(time,1); %create cell array to hold time-step matrices of population traits
tmp_pop_track_nu = zeros(N_mod, 10);
%C1 - 1 trait change for 90% of individuals // 
%C2 - delta of 1 trait change for 90% of individuals //
%C3 - trait 2 change for 10% of individuals //
%C4 - delta trait 2 change for 10% of individuals // 
%C5 - correlation that changed for 10% // 
%C6 - %correlation change for 10% of individuals // 
%C7 - %delta correlation change for 10% // 
%C8 - %distance of each individual from evolved destination //
%C9 - %fitness of individual
%C10 - individual index
track_ave_cor_changes = cell(time+1,1); %cell array to track population-wide correlation changes
ind_cor_changes = cell(N_mod,time+1); %cell array to track individual correlation changes 
sig_corr_index_cell = cell(N_mod,time+1); %cell array to track individual significant correlations
anc_sig_corrs_cell = cell(N_mod, time+1); %cell array to track individual significant correlation matrices
pop_cell = cell(time+1,1); %cell array for each population at every time step
Sr_ae_sim_cell = cell(time+1,1); %cell array for each projection at every time step in evolved
Sr_aa_sim_cell = cell(time+1,1); %cell array for each projection at every time step in ancestral
%sig_corr_cut_mat = zeros(N_mod, time+1); %matrix to keep track of cutoff for significant correlations
trait_a = zeros(N_mod, time+1); %track trait a changes in each individual over time
trait_b = zeros(N_mod, time+1); %track trait b changes in each individual over time
trait_c = zeros(N_mod, time+1); %track trait c changes in each individual over time
trait_d = zeros(N_mod, time+1); %track trait d changes in each individual over time

%%GLOBAL CELL MATRICES FOR ALL tgrad and cgrad values%%%%
pop_track_nu_all = cell(length(tgrad),length(cgrad));
track_ave_cor_changes_all = cell(length(tgrad),length(cgrad));
anc_sig_corrs_cell_all = cell(length(tgrad),length(cgrad));
pop_cell_all = cell(length(tgrad),length(cgrad));
Sr_ae_sim_cell_all = cell(length(tgrad),length(cgrad));
Sr_aa_sim_cell_all = cell(length(tgrad),length(cgrad));
trait_a_cell = cell(length(tgrad),length(cgrad));
trait_b_cell = cell(length(tgrad),length(cgrad));
trait_c_cell = cell(length(tgrad),length(cgrad));
trait_d_cell = cell(length(tgrad),length(cgrad));
%%%%%%%%%


if (strcmp(clonal,"TRUE") == 1) %create clonal population 
    
    %create clonal population
    %pop = dlmread('ancestral_starting_pop.txt');
    
    %clone_num = randsample(index_clone,1,false);
    clone_num = 20; %see above for notes on selecting this
    clone = ancestral_z(clone_num,:);
    %create 2 columns so that traits 1, 4, 5, and 6 align with column
    %numbers
    clone = [-0.7038 0 0 1.3192    1.3514   -0.3538];
    pop = repelem(clone, N_mod, 1); %clonal population
    %create new individual indices for clonal populations
    index_clone=(1:length(pop))';
    pop(:,end+1) = index_clone;
    
    %make a copy
    pop_og = pop;
    
    
    for i=1:1 
        trait_a(:,1) = pop(:,manually_selected_traits(i)); %store initial trait a values
        trait_b(:,1) = pop(:,manually_selected_traits(i+1)); %store initial trait b values
        trait_c(:,1) = pop(:,manually_selected_traits(i+2)); %store initial trait c values
        trait_d(:,1) = pop(:,manually_selected_traits(i+3)); %store initial trait d values
    end
    
elseif (strcmp(clonal,"TRUE") == 0) %create genetically diverse population
    
    for i=1:1 
        trait_a(:,1) = pop(:,manually_selected_traits(i)); %store initial trait a values
        trait_b(:,1) = pop(:,manually_selected_traits(i+1)); %store initial trait b values
        trait_c(:,1) = pop(:,manually_selected_traits(i+2)); %store initial trait c values
        trait_d(:,1) = pop(:,manually_selected_traits(i+3)); %store initial trait c values
    end
    
end


%store which traits in which individuals changed for this time step
%C1 - 1 trait change for 90% of individuals // 
%C2 - delta of 1 trait  change for 90% of individuals //
%C3 - trait 2 change for 10% of individuals //
%C4 - delta trait 2 change for 10% of individuals // 
%C5 - correlation that changed for 10% // 
%C6 - %correlation change for 10% of individuals // 
%C7 - %delta correlation change for 10% // 
%C8 - %distance of each individual from evolved destination //
%C9 - %fitness of individual
%C10 - Individual index

%create variable to save trait evolution figures
rep_count = 0;
for x=1:length(tgrad)  
    
    for y=1:length(cgrad)
        
        rep_count = rep_count + 1;
        
        for t=1:length(Driver)

            c=clock;
            pause(c(6)*runnumber/100)
            rng('shuffle')

            patch_tf = strcmp(patch,'FALSE'); %If true, employ 2 different selection regimes

            if patch_tf == 1

                select = 1; %only sample with weighting by fitness

            elseif (patch_tf == 0) && (Driver(i) == 1)


                select = 1; %resample population with weight by fitness


            elseif (patch_tf == 0) && (Driver(i) == 0)


                select = 0; %resample population with other weighting parameters


            end


            for k=1:length(grad_all)

                %%%Assign trait changes to individuals%%%
                index_mut=1:N_mod;

                %Sampling mutated individuals without replacement (multiple hits are not allowed)
                nu_total = nu_trait;
                index_mutated=randsample(index_mut,nu_total,false); 

                %generate vector of random trait changes
                tchange = tgrad(x)*randn(nu_total,1) + 0; %randn with mean of 0 and sd of 0.05 
                %tchange(tchange<0.01) = 0.01; %make sure none are below 0.01
                
                %if trait_flag == 0, set tchange equal to zero to not
                %change trait values
                if (trait_flag == 0)
                    
                    tchange = tchange.*0;
                    
                end

                %retrieve mutants from population 
                pop_nu_trait=pop(index_mutated,:); %change traits


                if (t == 1)
                    %find significantly correlated traits
                    sig_corr_index_cut = find(anc_corr_cell_sorted(:,4) < 0.05);

                    if (strcmp(corr_all,"TRUE") == 1)
                        %change all correlations
                        sig_corr_index = 1:length(anc_corr_cell_sorted);
                    else
                        %change only significant correlations
                        sig_corr_index = find(anc_corr_cell_sorted(:,4) < 0.05);

                    end
                    
                    %saves the correlations that we are allowing to change for every
                    %individual at every time step
                    for ind=1:N_mod
                        sig_corr_index_cell{ind,t} = sig_corr_index;
                        %repelem(sig_corr_index,1,N_mod)';
                    end
                    
                    
                    for ind=1:N_mod
                        
                       
                        
                        %%%%%%%%%%%%%generate random correlations between -1 and 1
                        a = -1;
                        b = 1;
                        anc_rand = (b-a).*rand(6,1) + a;
                        anc_corr_cell_sorted(:,3) = anc_rand;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %save individual correlation values
                        ind_cor_changes{ind,t} = anc_corr_cell_sorted;
                        
                        %save population-wide correlation values
                        anc_sig_corrs = anc_corr_cell_sorted(sig_corr_index,:);
                        anc_sig_corrs_cell{ind,t} = anc_sig_corrs;
                        
                        
                    end
                    
                    %obtain largest absolute correlation value just below
                    %significance cutoff to set a threshold for significant correlation
                    sig_corr_cut = abs(round(anc_corr_cell_sorted(length(sig_corr_index_cut)+1,3),1));


                elseif ( t > 1 )



                    %find significant correlation values above sig_corr_cut
                    for ind=1:N_mod
                        tmp_ind_cor_changes = ind_cor_changes{ind,t};

                        if (strcmp(corr_all,"TRUE") == 1)
                            %change all correlations
                            sig_cut_index = sig_corr_index; 
                        else
                            %change only significant correlations
                            sig_cut_index = find(tmp_ind_cor_changes(abs(tmp_ind_cor_changes(:,3)) > sig_corr_cut));
                        end



                        sig_corr_index_cell{ind,t} = sig_cut_index;

                        anc_sig_corrs_cell{ind,t} = tmp_ind_cor_changes(sig_cut_index,:);


                    end

                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %change trait a cell and recalculate all other traits
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %pull individual correlation values
                ind_cor_changes_nu_trait = ind_cor_changes(index_mutated,t);

                %pull individual significant correlation values
                ind_sig_corr_index = sig_corr_index_cell(index_mutated,t);

                %randomly choose which traits to change with replacement
                traits_to_change = (randsample(manually_selected_traits,nu_total,true))';


                %pop_delta_nu_trait contains delta change in traits
                %trait change contains new trait values 
                %pop_nu_trait contains the population with new trait values
                [pop_nu_trait, pop_delta_nu_trait, trait_change] = change_traits(nu_total,traits_to_change,pop_nu_trait, tchange, all_traits_num, ind_sig_corr_index, ind_cor_changes_nu_trait);

                %store new correlations for each individual
                ind_cor_changes(index_mutated,t+1) = ind_cor_changes_nu_trait;

                %store changed traits in their respective trait matrices
                for tr=1:length(traits_to_change)
                    if (traits_to_change(tr) == manually_selected_traits(1))
                        trait_a(index_mutated(tr),t+1) = trait_change(tr);
                    elseif (traits_to_change(tr) == manually_selected_traits(2))
                        trait_b(index_mutated(tr),t+1) = trait_change(tr);
                    elseif (traits_to_change(tr) == manually_selected_traits(3))
                        trait_c(index_mutated(tr),t+1) = trait_change(tr);
                    elseif (traits_to_change(tr) == manually_selected_traits(4))
                        trait_d(index_mutated(tr),t+1) = trait_change(tr);
                    end
                end

                %store which traits in which individuals changed for this time step
                %tmp_pop_track_nu = squeeze(pop_track_nu(t,:,:))';
                %C1 - 1 trait change for 90% of individuals // 
                %C2 - delta of 1 trait  change for 90% of individuals //
                %C3 - trait 2 change for 10% of individuals //
                %C4 - delta trait 2 change for 10% of individuals // 
                %C5 - correlation that changed for 10% // 
                %C6 - %correlation change for 10% of individuals // 
                %C7 - %delta correlation change for 10% // 
                %C8 - %distance of each individual from evolved destination //
                %C9 - %fitness of individual
                %C10 - %index
                tmp_pop_track_nu = zeros(N_mod, 10);
                tmp_pop_track_nu(index_mutated,1) = traits_to_change; %store traits that changed for each individual
                tmp_pop_track_nu(index_mutated,2) = pop_delta_nu_trait; %store delta change of traits for each individual
                tmp_pop_track_nu(index_mutated,10) = pop_nu_trait(:,end); %store individual indices

                %insert mutated individuals back into population
                pop(index_mutated,:) = pop_nu_trait;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %END change trait a cell and recalculate all other traits
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%CHANGE TRAIT AND CORRELATION%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %identify rest of individuals to change both correlation and traits
                index_mutated_hold = setdiff(index_mut,index_mutated);
                index_mutated_2 =randsample(index_mutated_hold,nu_trait2,false);
                
                %pull individual significant correlation values
                ind_sig_corr_index_2 = sig_corr_index_cell(index_mutated_2,t);

                %randomly choose the indices in sig_corr_index to determine which correlations to change with replacement
                cors_to_change = zeros(length(index_mutated_2),1);
                cors_that_changed = zeros(length(index_mutated_2),1);
                for ind=1:length(index_mutated_2)


                    if (strcmp(corr_all,"TRUE") == 1) %randomly select correlation to change

                        cors_to_change(ind) = randsample(length(anc_sig_corrs),1,true);
                        cors_that_changed(ind) = cors_to_change(ind);
                        ind_sig_corr_index_2{ind} = cors_to_change(ind);

                    elseif (isempty(ind_sig_corr_index_2{ind}) && strcmp(corr_all,"FALSE") == 1) %if no significant correlation, randomly select a correlation %to change

                        cors_to_change(ind) = randsample(length(anc_sig_corrs),1,true);
                        cors_that_changed(ind) = cors_to_change(ind);
                        ind_sig_corr_index_2{ind} = cors_to_change(ind);

                    elseif (~isempty(ind_sig_corr_index_2{ind}) && strcmp(corr_all,"FALSE") == 1) %select significant correlation

                        cors_to_change(ind) = randsample(length(ind_sig_corr_index_2{ind}),1,true);
                        tmp_ind_sig_corr_index_2 = ind_sig_corr_index_2{ind};
                        cors_that_changed(ind) = tmp_ind_sig_corr_index_2(cors_to_change(ind));
                    end

                end

                %generate vector of random correlation changes
                cchange = cgrad(y)*randn(nu_trait2,1) + 0;
                
                if (corr_flag == 0)
                    
                    cchange = cchange.*0;
                    
                end
                
                
                %generate vector of random trait changes
                tchange = tgrad(x)*randn(nu_trait2,1) + 0; %randn with mean of 1 and sd of 0.3 with a min of 0.01
                
                if (trait_flag == 0)
                    
                    tchange = tchange.*0;
                    
                end
                
                %retrieve from population
                pop_nu_trait2=pop(index_mutated_2,:); %change trait and correlation coefficient

                %pull individual correlation values
                ind_cor_changes_nu_trait_2 = ind_cor_changes(index_mutated_2,t);

                tmp_traits = zeros(length(cors_to_change),1);

                tmp_anc_sig_corrs_2 = anc_sig_corrs_cell(index_mutated_2,t);
                for r=1:length(cors_to_change)

                    %if tmp_anc_sig_corrs_2 empty
                    if isempty(tmp_anc_sig_corrs_2{r,1})
                        tmp_c = ind_cor_changes_nu_trait_2{r};
                        tmp_tmp_anc = tmp_c(cors_to_change(r),:);
                        tmp_traits(r) = tmp_tmp_anc(1);
                    else
                        %grab the trait in column1 of anc_sig_corrs of the randomly
                        %chosen indices in cors_to_change
                        tmp_tmp_anc = tmp_anc_sig_corrs_2{r,1};
                        tmp_traits(r) = tmp_tmp_anc(cors_to_change(r),1);
                    end
                end

                %store which traits in which individuals changed for this time step
                tmp_pop_track_nu(index_mutated_2,1) = tmp_traits; 

                %store correlations that changed for each individual
                tmp_pop_track_nu(index_mutated_2,5) = cors_that_changed; 

                %store individual indices
                tmp_pop_track_nu(index_mutated_2,10) = pop_nu_trait2(:,end); 

                %%%%change traits and correlations systematically%%%%%
                %pop_nu_trait2 contains the population with new trait values
                %pop_delta_nu_trait2_trait_1 contains delta change in 1st randomly
                %anc_delta_sig_corrs vector to save delta changes for correlations
                %selected trait
                %pop_delta_nu_trait2_trait_2 contains delta change in 2nd trait
                %that is correlated to the 1st randomly selected trait
                %trait_change_1 - trait change contains new trait 1 values 
                %trait_change_2 - trait change contains new trait 2 values 
                %tmp_traits_2 - vector to save 2nd changed trait
                %saved cor changes - save all correlation changes
                [pop_nu_trait2, pop_delta_nu_trait2_trait_1, pop_delta_nu_trait2_trait_2, anc_delta_sig_corrs, trait_change_1, trait_change_2, tmp_traits_2, ind_cor_changes_nu_trait_2] = change_traits_cors(nu_trait2,tmp_traits,pop_nu_trait2, tchange, cchange, cors_to_change, ind_sig_corr_index_2, ind_cor_changes_nu_trait_2);

                %store changed traits in their respective trait matrices
                for tr=1:length(tmp_traits)
                    if (tmp_traits(tr) == manually_selected_traits(1))
                        trait_a(index_mutated_2(tr),t+1) = pop_nu_trait2(tr,manually_selected_traits(1)); %store 1st changed trait
                    elseif (tmp_traits(tr) == manually_selected_traits(2))
                        trait_b(index_mutated_2(tr),t+1) = pop_nu_trait2(tr,manually_selected_traits(2)); %store 1st changed trait
                    elseif (tmp_traits(tr) == manually_selected_traits(3))
                        trait_c(index_mutated_2(tr),t+1) = pop_nu_trait2(tr,manually_selected_traits(3)); %store 1st changed trait
                    elseif (tmp_traits(tr) == manually_selected_traits(4))
                        trait_d(index_mutated_2(tr),t+1) = pop_nu_trait2(tr,manually_selected_traits(4)); %store 1st changed trait
                    end

                    if (tmp_traits_2(tr) == manually_selected_traits(1))
                        trait_a(index_mutated_2(tr),t+1) = pop_nu_trait2(tr,manually_selected_traits(1)); %store 2nd changed trait
                    elseif (tmp_traits_2(tr) == manually_selected_traits(2))
                        trait_b(index_mutated_2(tr),t+1) = pop_nu_trait2(tr,manually_selected_traits(2)); %store 2nd changed trait
                    elseif (tmp_traits_2(tr) == manually_selected_traits(3))
                        trait_c(index_mutated_2(tr),t+1) = pop_nu_trait2(tr,manually_selected_traits(3)); %store 2nd changed trait
                    elseif (tmp_traits_2(tr) == manually_selected_traits(4))
                        trait_d(index_mutated_2(tr),t+1) = pop_nu_trait2(tr,manually_selected_traits(4)); %store 2nd changed trait
                    end

                end

                %store delta change of traits for each individual
                tmp_pop_track_nu(index_mutated_2,2) = pop_delta_nu_trait2_trait_1;

                %store 2nd set of traits that changed
                tmp_pop_track_nu(index_mutated_2,3) = tmp_traits_2;

                %store delta change of 2nd set of traits for each individual
                tmp_pop_track_nu(index_mutated_2,4) = pop_delta_nu_trait2_trait_2;

                %store new correlations for each individual
                ind_cor_changes(index_mutated_2,t+1) = ind_cor_changes_nu_trait_2;
                
                %average population-wide correlation coefficients
                if (t==1)

                    %average of initial matrix
                    track_ave_cor_changes{t} = anc_corr_cell_sorted;

                    %average of first changes
                    %tmp_mat_ave = [];
                    c1 = zeros(N_mod,1);
                    c2 = zeros(N_mod,1);
                    c3 = zeros(N_mod,1);
                    c4 = zeros(N_mod,1);
                    c5 = zeros(N_mod,1);
                    c6 = zeros(N_mod,1);
                    for a=1:N_mod
                        tmp_mat = ind_cor_changes{a,t+1};
                        c1(a) = tmp_mat(1,3);
                        c2(a) = tmp_mat(2,3);
                        c3(a) = tmp_mat(3,3);
                        c4(a) = tmp_mat(4,3);
                        c5(a) = tmp_mat(5,3);
                        c6(a) = tmp_mat(6,3);
                    end

                    tmp_mat(1,3) = mean(c1);
                    tmp_mat(2,3) = mean(c2);
                    tmp_mat(3,3) = mean(c3);
                    tmp_mat(4,3) = mean(c4);
                    tmp_mat(5,3) = mean(c5);
                    tmp_mat(6,3) = mean(c6);

                    %store changed matrix
                    track_ave_cor_changes{t+1} = tmp_mat;

                elseif (t>1)

                    %average of first changes
                    %tmp_mat_ave = [];
                    c1 = zeros(N_mod,1);
                    c2 = zeros(N_mod,1);
                    c3 = zeros(N_mod,1);
                    c4 = zeros(N_mod,1);
                    c5 = zeros(N_mod,1);
                    c6 = zeros(N_mod,1);
                    for a=1:N_mod
                        tmp_mat = ind_cor_changes{a,t+1};
                        c1(a) = tmp_mat(1,3);
                        c2(a) = tmp_mat(2,3);
                        c3(a) = tmp_mat(3,3);
                        c4(a) = tmp_mat(4,3);
                        c5(a) = tmp_mat(5,3);
                        c6(a) = tmp_mat(6,3);
                    end

                    tmp_mat(1,3) = mean(c1);
                    tmp_mat(2,3) = mean(c2);
                    tmp_mat(3,3) = mean(c3);
                    tmp_mat(4,3) = mean(c4);
                    tmp_mat(5,3) = mean(c5);
                    tmp_mat(6,3) = mean(c6);

                    %store changed matrix
                    track_ave_cor_changes{t+1} = tmp_mat;


                end


                %store delta correlation change for each individual
                for c=1:length(anc_delta_sig_corrs)

                    if (anc_delta_sig_corrs(c,1) ~= 0)
                        tmp_pop_track_nu(index_mutated_2(c),7) = anc_delta_sig_corrs(c,1);
                    elseif (anc_delta_sig_corrs(c,2) ~= 0)
                        tmp_pop_track_nu(index_mutated_2(c),7) = anc_delta_sig_corrs(c,2);
                    elseif (anc_delta_sig_corrs(c,3) ~= 0)
                        tmp_pop_track_nu(index_mutated_2(c),7) = anc_delta_sig_corrs(c,3);
                    elseif (anc_delta_sig_corrs(c,4) ~= 0)
                        tmp_pop_track_nu(index_mutated_2(c),7) = anc_delta_sig_corrs(c,4);    
                    elseif (anc_delta_sig_corrs(c,5) ~= 0)
                        tmp_pop_track_nu(index_mutated_2(c),7) = anc_delta_sig_corrs(c,5);     
                    elseif (anc_delta_sig_corrs(c,6) ~= 0)
                        tmp_pop_track_nu(index_mutated_2(c),7) = anc_delta_sig_corrs(c,6);     

                    end
                end

                %insert mutated individuals back into population
                pop(index_mutated_2,:) = pop_nu_trait2;


            end %tgrad end

            %EVOLVED SIMULATED POPULATION
            if (t==1)
             
                
                pop_og_proj = pop_og;
                pop_og_proj( :, ~any(pop_og_proj,1) ) = [];  %delete columns with all zeros
                pop_og_proj = pop_og_proj(:,1:length(all_traits_num));
                
                Sr_aa_sim_og=pop_og_proj*Ar_a_dat; %ancestral PCA space
                Sr_ae_sim_og=pop_og_proj*Ar_e_dat; %evolved PCA space
                

            end

            %PROJECT POPULATION ONTO EVOLVED AXES%%%%%
            clear sim_pca_strains
            for s=1:N_mod
                sim_pca_strains{s,1} = pop(s,end);
            end
            pop_proj = pop;
            pop_proj( :, ~any(pop_proj,1) ) = [];  %delete columns with all zeros
            pop_proj = pop_proj(:,1:length(all_traits_num));

            Sr_aa_sim=pop_proj*Ar_a_dat;      % Ancestral Factor Score
            Sr_ae_sim=pop_proj*Ar_e_dat;        % Ancestral plotted onto evolved Factor Score

            %save 1st projection before reproduction
            if (t==1)
                pop_cell{t} = pop_og_proj; %save trait values of starting population
                Sr_aa_sim_cell{t} = Sr_aa_sim_og; %ancestral space
                Sr_ae_sim_cell{t} = Sr_ae_sim_og; %evolved space

            end
            %choose a single individual from the simulated evolved populations
            %calculate distances from simulated population to single individual in
            %evolved PC space
            if (strcmp(clonal,"TRUE")==1)


                if (t==1)

                    %Based on PCA above, we selected individual from row 20
                    % of evolved_z population
                    %sim_e_index=1:length(Sr_e_sim);
                    %sim_e_index = 1;
                    sim_clone_num = sim_e_index;
                    sim_e_clone = Sr_e(sim_clone_num,:);  
                    dlmwrite('sim_e_clone.txt',sim_e_clone);

                end

                %calculate distances of current population from selected evolved
                %individual
                Sr_sim_dist_a = []; %ancestral space
                Sr_sim_dist = []; %evolved space
                for i=1:length(Sr_ae_sim)
                    
                    tmp_a = [Sr_aa_sim(i,1),Sr_aa_sim(i,2);sim_e_clone(1,1),sim_e_clone(1,2)]; %ancestral
                    Sr_sim_dist_a(i) = pdist(tmp_a,'euclidean'); %distances between ancestral and evolved points in ancestral PCA
                    
                    tmp = [Sr_ae_sim(i,1),Sr_ae_sim(i,2);sim_e_clone(1,1),sim_e_clone(1,2)];
                    Sr_sim_dist(i) = pdist(tmp,'euclidean'); %distances between ancestral and evolved points in evolved PCA
                end

                tmp_pop_track_nu(:,8) = Sr_sim_dist';

            end



            %CALCULATE FITNESS HERE FIRST
            tmp_pop_track_nu(:,9) = exp((-Sr_sim_dist'.^2)/2);
            
            %insert final index into tmp_pop_track_nu
            %tmp_pop_track_nu(index_mutated_3,10) = pop(index_mutated_3,end);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Population regulation and reproduction
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Soft selection is acting -> population size stays the same
            %Sampling individuals to the next time point with replacement

            ofsp_index = 1:N_mod;

            if select == 0

                mode_type = strcmp(mode,'random'); 

                if mode_type == 1

                    %Sampling the next generation, no weighting
                    c=clock;
                    pause(c(6)*runnumber/100)
                    rng('shuffle')
                    sampled_index = randsample(ofsp_index, N_mod, true);

                else

                    %Sampling the next generation, negative weighting on distance
                    %to optimum
                    tmp = 1./(tmp_pop_track_nu(:,8)+0.1);
                    tmp2=(tmp-mean(tmp))./std(tmp);
                    reciprocal_weights=(tmp2-min(tmp2))./max(tmp2-min(tmp2))*tfact;

                    %w=tmp/sum(tmp);
                    %max(w);

                    c=clock;
                    pause(c(6)*runnumber/100)
                    rng('shuffle')

                    sampled_index = randsample(ofsp_index, N_mod, true, reciprocal_weights);

                end

            elseif select == 1
                
                mode_type = strcmp(mode,'random'); %If true, sample without weighting

                if mode_type == 1

                    %Sampling the next generation, no weighting
                    c=clock;
                    pause(c(6)*runnumber/100)
                    rng('shuffle')
                    sampled_index = randsample(ofsp_index, N_mod, true);
                else
                    %Sampling the next generation, weighting with fitness
                    c=clock;
                    pause(c(6)*runnumber/100)
                    rng('shuffle')
                    sampled_index = randsample(ofsp_index, N_mod, true, tmp_pop_track_nu(:,9));
                    
                end
            end

            %Population after reproduction and regulation
            tmp_pop_track_nu = tmp_pop_track_nu(sampled_index,:);
            pop = pop(sampled_index,:);
            trait_a = trait_a(sampled_index,:);
            trait_b = trait_b(sampled_index,:);
            trait_c = trait_c(sampled_index,:);
            trait_d = trait_d(sampled_index,:);
            ind_cor_changes = ind_cor_changes(sampled_index,:);
            %track_ave_cor_changes = track_ave_cor_changes(sampled_index,:);
            sig_corr_index_cell = sig_corr_index_cell(sampled_index,:);
            anc_sig_corrs_cell = anc_sig_corrs_cell(sampled_index,:);

            %insert time step matrix after reproduction into global cell array
            pop_track_nu{t} = tmp_pop_track_nu;

            %calculate new ancestral projects following reproduction
            pop_proj = pop;
            pop_proj( :, ~any(pop_proj,1) ) = [];  %delete columns with all zeros
            pop_proj = pop_proj(:,1:length(all_traits_num));
            Sr_aa_sim=pop_proj*Ar_a_dat;        % Ancestral plotted onto ancestral Factor Score
            Sr_ae_sim=pop_proj*Ar_e_dat;        % Ancestral plotted onto evolved Factor Score

            %save each projection into a cell array
            pop_cell{t+1} = pop_proj; %save population trait values
            Sr_aa_sim_cell{t+1} = Sr_aa_sim; %ancestral
            Sr_ae_sim_cell{t+1} = Sr_ae_sim; %evolved


            %%%END DO SUBSEQUENT PCAS WITH THE SIMULATED POPULATION DATA%%%%%

        end
        
    end %end cgrad
    
    %save all pop track nu
    tmp_pop_track_nu_comb = vertcat(pop_track_nu{:});
    pop_track_nu_all{x,y} = tmp_pop_track_nu_comb;

    %save all track_ave_cor_changes
    tmp_track_ave_cor_changes = vertcat(track_ave_cor_changes{:});
    track_ave_cor_changes_all{x,y} = tmp_track_ave_cor_changes;

    %save all anc_sig_corrs_cell
    tmp_anc_sig_corrs_cell = vertcat(anc_sig_corrs_cell{:});
    anc_sig_corrs_cell_all{x,y} = tmp_anc_sig_corrs_cell;

    %save all pop_cell
    tmp_pop_cell = vertcat(pop_cell{:});
    pop_cell_all{x,y} = tmp_pop_cell;

    %save all Sr_ae_sim_cell
    tmp_Sr_ae_sim_cell = vertcat(Sr_ae_sim_cell{:});
    Sr_ae_sim_cell_all{x,y} = tmp_Sr_ae_sim_cell;

    %save all Sr_aa_sim_cell
    tmp_Sr_aa_sim_cell = vertcat(Sr_aa_sim_cell{:});
    Sr_aa_sim_cell_all{x,y} = tmp_Sr_aa_sim_cell;

    %save all trait a
    trait_a_cell{x,y} = trait_a;

    %save all_trait b
    trait_b_cell{x,y} = trait_b;

    %save all trait c
    trait_c_cell{x,y} = trait_c;

    %save all trait d
    trait_d_cell{x,y} = trait_d;

    
    
    
end %end tgrad



