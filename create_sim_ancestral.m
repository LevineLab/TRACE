function pop = create_sim_ancestral(N_mod, all_traits_num, ancestral_traits_change, anc_corr_cell_sorted, ancestral_trait_stats)


pop = zeros(N_mod,length(all_traits_num));
traits_already_done = [];
%traits_already_done(1) = trait_to_change;
for i=1:length(anc_corr_cell_sorted)

    %current_pair = strsplit(anc_corr_cell_sorted{i,1}, '-');

    %tmp_trait_a = str2double( regexp(current_pair{1},'t(\d+)','tokens','once') );
    %tmp_trait_b = str2double( regexp(current_pair{2},'t(\d+)','tokens','once') );
    tmp_trait_a = anc_corr_cell_sorted(i,1);
    tmp_trait_b = anc_corr_cell_sorted(i,2);
    
    trait_a_check = ismember(tmp_trait_a, traits_already_done); %check if trait a is already analyzed
    trait_b_check = ismember(tmp_trait_b, traits_already_done); %check if trait a is already analyzed

    if (trait_a_check == 1 && trait_b_check == 1)
        %i
        %disp('hi')

        continue

    end

    if (anc_corr_cell_sorted(i,4) < 0.05)

        if (trait_a_check == 1 && trait_b_check == 0)

            ta_rand = pop(:,tmp_trait_a); %pull already calculated trait values

            %Linear regression from empirical data to generate correlated trait
            tmp_a = ancestral_traits_change(:,tmp_trait_a); 
            tmp_b = ancestral_traits_change(:,tmp_trait_b); 
            p = polyfit(tmp_a,tmp_b,1); %polyfit to compute a linear regression that predicts y from x
            %p(1) is slope and p(2) is y-intercept 
            tb_rand = polyval(p,ta_rand); %calculate y from the model
            pop(:,tmp_trait_b) = tb_rand; %populate appropriate column with randomized trait values

            if (i == 1) 
                trait_index_holder = i+1;
            elseif (i > 1)
                trait_index_holder = length(traits_already_done);
            end

            if (i > 1)
                traits_already_done(i+trait_index_holder) = tmp_trait_b; %document trait measured
            else
                traits_already_done(i+1) = tmp_trait_b; %document trait measured
            end




        elseif (trait_a_check == 0 && trait_b_check == 1) %if trait a is uknown but trait b is known

            tb_rand = pop(:,tmp_trait_b); %pull already calculated trait values

            %Linear regression from empirical data to generate correlated trait
            tmp_a = ancestral_traits_change(:,tmp_trait_a); 
            tmp_b = ancestral_traits_change(:,tmp_trait_b); 
            p = polyfit(tmp_a,tmp_b,1); %polyfit to compute a linear regression that predicts y (second number) from x
            ta_rand = (tb_rand - p(2))/p(1);%get ta_rand from tb_rand
            pop(:,tmp_trait_a) = ta_rand; %populate appropriate column with randomized trait values

            if (i == 1) 
                trait_index_holder = i+1;
            elseif (i > 1)
                trait_index_holder = length(traits_already_done);
            end

            if (i > 1)
                traits_already_done(i+trait_index_holder) = tmp_trait_a; %document trait measured
            else
                traits_already_done(i+1) = tmp_trait_a; %document trait measured
            end



        elseif (trait_a_check == 0 && trait_b_check == 0) %if trait a is uknown but trait b is known
            
            %generate normal distribution with sd and mean from
            %standardized ancestral data
            mean_a = ancestral_trait_stats(tmp_trait_a,1);
            sd_a = ancestral_trait_stats(tmp_trait_a,2);
            ta_rand = sd_a.*randn(N_mod,1) + mean_a;
            pop(:,tmp_trait_a) = ta_rand;

            if (i == 1) 
                trait_index_holder = i+1;
            elseif (i > 1)
                trait_index_holder = length(traits_already_done);
            end

            if (i > 1)
                traits_already_done(i+trait_index_holder) = tmp_trait_a; %document trait measured
            else
                traits_already_done(i+1) = tmp_trait_a; %document trait measured
            end



            %Linear regression from empirical data to generate correlated trait
            tmp_a = ancestral_traits_change(:,tmp_trait_a); 
            tmp_b = ancestral_traits_change(:,tmp_trait_b); 
            p = polyfit(tmp_a,tmp_b,1); %polyfit to compute a linear regression that predicts y from x
            %p(1) is slope and p(2) is y-intercept 
            tb_rand = polyval(p,ta_rand); %calculate y from the model
            pop(:,tmp_trait_b) = tb_rand;

            if (i == 1) 
                trait_index_holder = i+2;
            elseif (i > 1)
                trait_index_holder = length(traits_already_done);
            end

            traits_already_done(i+trait_index_holder) = tmp_trait_b; %document trait measured



        end

    elseif (anc_corr_cell_sorted(i,4) >= 0.05) %generate random numbers if significant correlation doesnt exist

        
        if (trait_a_check == 1 && trait_b_check == 0)

            mean_b = ancestral_trait_stats(tmp_trait_b,1);
            sd_b = ancestral_trait_stats(tmp_trait_b,2);
            tb_rand = sd_b.*randn(N_mod,1) + mean_b;
            pop(:,tmp_trait_b) = tb_rand;

            if (i == 1) 
                trait_index_holder = i+1;
            elseif (i > 1)
                trait_index_holder = length(traits_already_done);
            end

            if (i > 1)
                traits_already_done(i+trait_index_holder) = tmp_trait_b; %document trait measured
            else
                traits_already_done(i+1) = tmp_trait_b; %document trait measured
            end




        elseif (trait_a_check == 0 && trait_b_check == 1) %if trait a is uknown but trait b is known

            mean_a = ancestral_trait_stats(tmp_trait_a,1);
            sd_a = ancestral_trait_stats(tmp_trait_a,2);
            ta_rand = sd_a.*randn(N_mod,1) + mean_a;
            pop(:,tmp_trait_a) = ta_rand;

            if (i == 1) 
                trait_index_holder = i+1;
            elseif (i > 1)
                trait_index_holder = length(traits_already_done);
            end

            if (i > 1)
                traits_already_done(i+trait_index_holder) = tmp_trait_a; %document trait measured
            else
                traits_already_done(i+1) = tmp_trait_a; %document trait measured
            end


        elseif (trait_a_check == 0 && trait_b_check == 0) %if trait a is uknown but trait b is known

            mean_a = ancestral_trait_stats(tmp_trait_a,1);
            sd_a = ancestral_trait_stats(tmp_trait_a,2);
            ta_rand = sd_a.*randn(N_mod,1) + mean_a;
            pop(:,tmp_trait_a) = ta_rand;

            if (i == 1) 
                trait_index_holder = i+1;
            elseif (i > 1)
                trait_index_holder = length(traits_already_done);
            end

            if (i > 1)
                traits_already_done(i+trait_index_holder) = tmp_trait_a; %document trait measured
            else
                traits_already_done(i+1) = tmp_trait_a; %document trait measured
            end



            mean_b = ancestral_trait_stats(tmp_trait_b,1);
            sd_b = ancestral_trait_stats(tmp_trait_b,2);
            tb_rand = sd_b.*randn(N_mod,1) + mean_b;
            pop(:,tmp_trait_b) = tb_rand;

            if (i == 1) 
                trait_index_holder = i+1;
            elseif (i > 1)
                trait_index_holder = length(traits_already_done);
            end

            if (i > 1)
                traits_already_done(i+trait_index_holder) = tmp_trait_b; %document trait measured
            else
                traits_already_done(i+1) = tmp_trait_b; %document trait measured
            end




        end


    end



end