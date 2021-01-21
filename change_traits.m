function [pop_nu_trait, pop_delta_nu_trait, trait_change] = change_traits(nu_total,traits_to_change,pop_nu_trait, tchange, all_traits_num, ind_sig_corr_index, ind_cor_changes_nu_trait)

pop_delta_nu_trait = zeros(nu_total,1); %vector to save delta changes for traits
trait_change = zeros(nu_total,1); %vector to save changed traits
for t1=1:nu_total
    
    
    %keep track inside below loop of systematically changing traits
    done_traits = []; 

    %randomly choose trait
    done_traits(1) = traits_to_change(t1);

    %change trait
    trait_change(t1) = pop_nu_trait(t1,traits_to_change(t1)) + tchange(t1); 
    
    %save the delta change
    pop_delta_nu_trait(t1) = tchange(t1);
    %pop_delta_nu_trait(t1) = abs(abs(pop_nu_trait(t1,traits_to_change(t1))) - abs(trait_change(t1)));
    
    %insert new trait value
    pop_nu_trait(t1,traits_to_change(t1)) = trait_change(t1);
    
    %change other trait values based on significantly correlated
    %traits
    
    anc_sig_corrs_tmp = ind_cor_changes_nu_trait{t1,1};
    sig_corr_index = ind_sig_corr_index{t1,1};
    
    anc_sig_corrs_tmp = anc_sig_corrs_tmp(sig_corr_index,:);
    
    if (find(ismember(traits_to_change(t1),anc_sig_corrs_tmp(:,1))) >= 1)

        match_index = find(ismember(anc_sig_corrs_tmp(:,1),traits_to_change(t1)));

        for m=1:length(match_index)

            pop_nu_trait(t1,anc_sig_corrs_tmp(match_index(m),2)) = pop_nu_trait(t1,traits_to_change(t1))*anc_sig_corrs_tmp(match_index(m),3);
            dt_index_tmp = length(done_traits);
            done_traits(dt_index_tmp+m) = anc_sig_corrs_tmp(match_index(m),2);

        end
    end

    if (find(ismember(traits_to_change(t1), anc_sig_corrs_tmp(:,2))) >= 1)

        match_index = find(ismember(anc_sig_corrs_tmp(:,2),traits_to_change(t1)));


        for m=1:length(match_index)

            if (ismember(anc_sig_corrs_tmp(match_index(m),1),done_traits))
                continue
            else
                pop_nu_trait(t1,anc_sig_corrs_tmp(match_index(m),1)) = pop_nu_trait(t1,traits_to_change(t1))/anc_sig_corrs_tmp(match_index(m),3);
                dt_index_tmp = length(done_traits);
                done_traits(dt_index_tmp+m) = anc_sig_corrs_tmp(match_index(m),1);
            end
        end

    end

    done_traits = done_traits(done_traits~=0);

    if length(done_traits) > 1
    %for c=1:length(all_traits_num)
        for c=1:length(done_traits) %change all other traits

            if (ismember(done_traits(c),anc_sig_corrs_tmp(:,1)))

                done_trait_index = find(ismember(anc_sig_corrs_tmp(:,1),done_traits(c)));

                for d=1:length(done_trait_index)

                    if (ismember(anc_sig_corrs_tmp(done_trait_index(d),2),done_traits))

                        continue

                    else

                        pop_nu_trait(t1,anc_sig_corrs_tmp(done_trait_index(d),2)) = pop_nu_trait(t1,anc_sig_corrs_tmp(done_trait_index(d),1))*anc_sig_corrs_tmp(done_trait_index(d),3);
                        dt_index_tmp = length(done_traits);
                        done_traits(dt_index_tmp+d) = anc_sig_corrs_tmp(done_trait_index(d),2);
                        done_traits = done_traits(done_traits~=0);

                    end

                end


            end 

            if (ismember(done_traits(c),(anc_sig_corrs_tmp(:,2))))

                done_trait_index = find(ismember(anc_sig_corrs_tmp(:,2),done_traits(c)));

                for d=1:length(done_trait_index)

                    if (ismember(anc_sig_corrs_tmp(done_trait_index(d),1),done_traits))

                        continue

                    else

                        pop_nu_trait(t1,anc_sig_corrs_tmp(done_trait_index(d),1)) = pop_nu_trait(t1,anc_sig_corrs_tmp(done_trait_index(d),2))/anc_sig_corrs_tmp(done_trait_index(d),3);
                        dt_index_tmp = length(done_traits);
                        done_traits(dt_index_tmp+d) = anc_sig_corrs_tmp(done_trait_index(d),1);
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
    end
    
end