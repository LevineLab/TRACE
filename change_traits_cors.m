function [pop_nu_trait2, pop_delta_nu_trait2_trait_1, pop_delta_nu_trait2_trait_2, anc_delta_sig_corrs, trait_change_1, trait_change_2, tmp_traits_2, ind_cor_changes_nu_trait_2] = change_traits_cors(nu_trait2,tmp_traits,pop_nu_trait2, tchange, cchange, cors_to_change, ind_sig_corr_index_2, ind_cor_changes_nu_trait_2)

pop_delta_nu_trait2_trait_1 = zeros(nu_trait2,1); %vector to save delta changes for trait 1
pop_delta_nu_trait2_trait_2 = zeros(nu_trait2,1); %vector to save delta changes for trait 2
%saved_cor_changes = zeros(length(cors_to_change),max(cors_to_change)); %save all correlation changes
anc_delta_sig_corrs = zeros(nu_trait2,max(cors_to_change)); %vector to save delta changes for correlations
trait_change_1 = zeros(nu_trait2,1); %vector to save trait changes
trait_change_2 = zeros(nu_trait2,1); %vector to save trait 2 changes
tmp_traits_2 = zeros(nu_trait2,1); %vector to save 2nd changed trait 

for t1=1:nu_trait2
            
    %keep track inside below loop of systematically changing traits
    done_traits = []; 

    %randomly choose trait
    done_traits(1) = tmp_traits(t1);
    
    %pull significantly correlations to change
    sig_corr_index = ind_sig_corr_index_2{t1,1};
    
    %make copy of original correlation coefficients to systematically
    %change correlation coefficients for each individual
    %if (t1 == 1)
    anc_sig_corrs_hold = ind_cor_changes_nu_trait_2{t1,1};
    anc_sig_corrs_hold = anc_sig_corrs_hold(sig_corr_index,:);
    
    
        %anc_sig_corrs_hold = anc_sig_corrs;
        %ind_cor_changes_nu_trait_2_hold = ind_cor_changes_nu_trait_2;
    %end
    
    %change correlation coefficient
    cor_change = anc_sig_corrs_hold(1,3) + cchange(t1);
    %anc_sig_corrs(cors_to_change(t1),3) = anc_sig_corrs(cors_to_change(t1),3)*cchange(t1);
    anc_sig_corrs_hold(1,3) = cor_change;
    
    %save changed correlation coefficient
    %saved_cor_changes(t1,cors_to_change(t1)) = cor_change;
    
    %save the delta cor change
    anc_delta_sig_corrs(t1, cors_to_change(t1)) = cchange(t1);
    %anc_delta_sig_corrs(t1, cors_to_change(t1)) = abs(anc_sig_corrs_hold(cors_to_change(t1),3) - cor_change);
        
    %change trait
    trait_change_1(t1) = pop_nu_trait2(t1,tmp_traits(t1)) + tchange(t1);
    
    %pop_nu_trait2(t1,tmp_traits(t1)) = pop_nu_trait2(t1,tmp_traits(t1))*tchange(t1);

    %save the delta trait change
    pop_delta_nu_trait2_trait_1(t1) = tchange(t1);
    %pop_delta_nu_trait2_trait_1(t1) = abs(pop_nu_trait2(t1,tmp_traits(t1)) - trait_change_1(t1));
    
    %insert changed trait
    pop_nu_trait2(t1,tmp_traits(t1)) = trait_change_1(t1); 
    %change other trait values based on significantly correlated
    %traits

    if (find(ismember(tmp_traits(t1),anc_sig_corrs_hold(:,1))) >= 1)

        match_index = find(ismember(anc_sig_corrs_hold(:,1),tmp_traits(t1)));

        for m=1:length(match_index)
            
            tmp_traits_2(t1) = anc_sig_corrs_hold(match_index(m),2); %vector to save second changed trait
            trait_change_2(t1) = pop_nu_trait2(t1,tmp_traits(t1))*cor_change; %vector to save changed values for trait2
            pop_delta_nu_trait2_trait_2(t1) = trait_change_2(t1);
            %pop_delta_nu_trait2_trait_2(t1) = abs(pop_nu_trait2(t1,tmp_traits_2(t1)) - trait_change_2(t1)); %save delta trait change
            
            pop_nu_trait2(t1,anc_sig_corrs_hold(match_index(m),2)) = trait_change_2(t1); 
            dt_index_tmp = length(done_traits);
            done_traits(dt_index_tmp+m) = anc_sig_corrs_hold(match_index(m),2);

        end
    end

    if (find(ismember(tmp_traits(t1), anc_sig_corrs_hold(:,2))) >= 1)

        match_index = find(ismember(anc_sig_corrs_hold(:,2),tmp_traits(t1)));


        for m=1:length(match_index)

            if (ismember(anc_sig_corrs_hold(match_index(m),1),done_traits))
                continue
            else
                tmp_traits_2(t1) = anc_sig_corrs_hold(match_index(m),2); %vector to save second changed trait
                trait_change_2(t1) = pop_nu_trait2(t1,tmp_traits(t1))/cor_change; %vector to save changed values for trait2
                pop_delta_nu_trait2_trait_2(t1) = trait_change_2(t1);
                %pop_delta_nu_trait2_trait_2(t1) = abs(pop_nu_trait2(t1,tmp_traits_2(t1)) - trait_change_2(t1)); %save delta trait change
                
                pop_nu_trait2(t1,anc_sig_corrs_hold(match_index(m),1)) = trait_change_2(t1);
                dt_index_tmp = length(done_traits);
                done_traits(dt_index_tmp+m) = anc_sig_corrs_hold(match_index(m),1);
            end
        end

    end

    done_traits = done_traits(done_traits~=0);
    
    %for c=1:length(all_traits_num)
    for c=1:length(done_traits) %change all other traits

        if (ismember(done_traits(c),anc_sig_corrs_hold(:,1)))

            done_trait_index = find(ismember(anc_sig_corrs_hold(:,1),done_traits(c)));

            for d=1:length(done_trait_index)

                if (ismember(anc_sig_corrs_hold(done_trait_index(d),2),done_traits))

                    continue

                else
                    tmp_traits_2(t1) = anc_sig_corrs_hold(done_trait_index(d),2); %vector to save second changed trait
                    trait_change_2(t1) = pop_nu_trait2(t1,anc_sig_corrs_hold(done_trait_index(d),1))*cor_change; %vector to save changed values for trait2
                    pop_delta_nu_trait2_trait_2(t1) = trait_change_2(t1);
                    pop_delta_nu_trait2_trait_2(t1) = abs(pop_nu_trait2(t1,tmp_traits_2(t1)) - trait_change_2(t1)); %save delta trait change
                    
                    pop_nu_trait2(t1,anc_sig_corrs_hold(done_trait_index(d),2)) = trait_change_2(t1);
                    dt_index_tmp = length(done_traits);
                    done_traits(dt_index_tmp+d) = anc_sig_corrs_hold(done_trait_index(d),2);
                    done_traits = done_traits(done_traits~=0);

                end

            end


        end 

        if (ismember(done_traits(c),(anc_sig_corrs_hold(:,2))))

            done_trait_index = find(ismember(anc_sig_corrs_hold(:,2),done_traits(c)));

            for d=1:length(done_trait_index)

                if (ismember(anc_sig_corrs_hold(done_trait_index(d),1),done_traits))

                    continue

                else
                    
                    tmp_traits_2(t1) = anc_sig_corrs_hold(done_trait_index(d),2); %vector to save second changed trait
                    trait_change_2(t1) = pop_nu_trait2(t1,anc_sig_corrs_hold(done_trait_index(d),1))/cor_change; %vector to save changed values for trait2
                    pop_delta_nu_trait2_trait_2(t1) = trait_change_2(t1);
                    pop_delta_nu_trait2_trait_2(t1) = abs(pop_nu_trait2(t1,tmp_traits_2(t1)) - trait_change_2(t1)); %save delta trait change
                    
                    pop_nu_trait2(t1,anc_sig_corrs_hold(done_trait_index(d),1)) = trait_change_2(t1);
                    dt_index_tmp = length(done_traits);
                    done_traits(dt_index_tmp+d) = anc_sig_corrs_hold(done_trait_index(d),1);
                    done_traits = done_traits(done_traits~=0);

                end

            end

        end

        %break out of loop if all traits have been changed
        %if (length(done_traits) == length(all_traits_num))
        if (length(done_traits) == length(sig_corr_index))

            break

        end

    end
    
    %replaced changed correlation values back into mutated individual
    tmp_anc_sig_corrs = ind_cor_changes_nu_trait_2{t1,1};
    tmp_anc_sig_corrs(sig_corr_index,:) = anc_sig_corrs_hold;
    ind_cor_changes_nu_trait_2{t1,1} = tmp_anc_sig_corrs;
    
end

%ave_cor_changes = zeros(max(cors_to_change),1);
%for i=1:max(cors_to_change)
    
%    tmp_index = saved_cor_changes(:,i)~=0;
%    ave_cor_changes(i) = mean(saved_cor_changes(tmp_index,i));
   
    
%end


