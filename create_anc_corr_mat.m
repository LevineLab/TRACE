function anc_corr_cell_sorted = create_anc_corr_mat(ancestral_z,  all_traits_num)

[R,P]=corrcoef(ancestral_z); %R has correlation coefficients; P has Ps-values
R = triu(R,1); %upper triangle only
P = triu(P,1); %upper triangle only
[I, J] = find(R ~= 0); %Find indices 
cor_mat = [I, J];
cor_mat = sortrows(cor_mat,1);
anc_corr_cell = []; %cell array for correlated traits coefficients
abs_anc_corr_cell = {}; %cell array for absolute values of correlated traits coefficients

for i=1:length(cor_mat)
    
    if (cor_mat(i,1) == 1) %trait 1 vs the world 
        
        if (cor_mat(i,2)==i+1)
            
            %anc_corr_cell(i,:) = {[all_traits_num{1,1} '-' all_traits_num{i+1,1}],R(1,i+1), P(1,i+1)};
            anc_corr_cell(i,1) = str2double(all_traits_num{1,1}); 
            anc_corr_cell(i,2) = str2double(all_traits_num{i+1,1});
            anc_corr_cell(i,3) = R(1,i+1);
            anc_corr_cell(i,4) = P(1,i+1);
            abs_anc_corr_cell(i,:) = {[all_traits_num{1,1} '-' all_traits_num{i+1,1}],abs(R(1,i+1)), P(1,i+1)};
            
        end
        
    elseif (cor_mat(i,1) == 2) %trait 2 vs the world
        
            
        %anc_corr_cell(i,:) = {[all_traits_num{2,1} '-' all_traits_num{cor_mat(i,1)+1,1}],R(2,cor_mat(i,1)+1), P(2,cor_mat(i,1)+1)};
        %anc_corr_cell(i,:) = {[all_traits_num{2,1} '-' all_traits_num{cor_mat(i,2),1}],R(2,cor_mat(i,2)), P(2,cor_mat(i,2))};
        anc_corr_cell(i,1) = str2double(all_traits_num{2,1});
        anc_corr_cell(i,2) = str2double(all_traits_num{cor_mat(i,2),1});
        anc_corr_cell(i,3) = R(2,cor_mat(i,2));
        anc_corr_cell(i,4) = P(2,cor_mat(i,2));
        abs_anc_corr_cell(i,:) = {[all_traits_num{2,1} '-' all_traits_num{cor_mat(i,2),1}],abs(R(2,cor_mat(i,2))), P(2,cor_mat(i,2))};
        
    
    elseif (cor_mat(i,1) == 3) %trait 3 vs the world
        
        %anc_corr_cell(i,:) = {[all_traits_num{3,1} '-' all_traits_num{cor_mat(i,2),1}],R(3,cor_mat(i,2)), P(3,cor_mat(i,2))};
        anc_corr_cell(i,1) = str2double(all_traits_num{3,1});
        anc_corr_cell(i,2) = str2double(all_traits_num{cor_mat(i,2),1});
        anc_corr_cell(i,3) = R(3,cor_mat(i,2));
        anc_corr_cell(i,4) = P(3,cor_mat(i,2));
        abs_anc_corr_cell(i,:) = {[all_traits_num{3,1} '-' all_traits_num{cor_mat(i,2),1}],abs(R(3,cor_mat(i,2))), P(3,cor_mat(i,2))};
        %anc_corr_cell(i,:) = {[all_traits_num{3,1} '-' all_traits_num{i-2,1}],R(3,i-2), P(3,i-2)};
        %abs_anc_corr_cell(i,:) = {[all_traits_num{3,1} '-' all_traits_num{i-2,1}],abs(R(3,i-2)), P(3,i-2)};
    
    elseif (cor_mat(i,1) == 4) %trait 4 vs the world
        
        %anc_corr_cell(i,:) = {[all_traits_num{4,1} '-' all_traits_num{cor_mat(i,2),1}],R(4,cor_mat(i,2)), P(4,cor_mat(i,2))};
        anc_corr_cell(i,1) = str2double(all_traits_num{4,1});
        anc_corr_cell(i,2) = str2double(all_traits_num{cor_mat(i,2),1});
        anc_corr_cell(i,3) = R(4,cor_mat(i,2));
        anc_corr_cell(i,4) = P(4,cor_mat(i,2));
        abs_anc_corr_cell(i,:) = {[all_traits_num{4,1} '-' all_traits_num{cor_mat(i,2),1}],abs(R(4,cor_mat(i,2))), P(4,cor_mat(i,2))};
    
    elseif (cor_mat(i,1) == 5) %trait 5 vs the world

        %anc_corr_cell(i,:) = {[all_traits_num{5,1} '-' all_traits_num{cor_mat(i,2),1}],R(5,cor_mat(i,2)), P(5,cor_mat(i,2))};
        anc_corr_cell(i,1) = str2double(all_traits_num{5,1});
        anc_corr_cell(i,2) = str2double(all_traits_num{cor_mat(i,2),1});
        anc_corr_cell(i,3) = R(5,cor_mat(i,2));
        anc_corr_cell(i,4) = P(5,cor_mat(i,2));
        abs_anc_corr_cell(i,:) = {[all_traits_num{5,1} '-' all_traits_num{cor_mat(i,2),1}],abs(R(5,cor_mat(i,2))), P(5,cor_mat(i,2))};
    
    elseif (cor_mat(i,1) == 6) %trait 6 vs the world

        %anc_corr_cell(i,:) = {[all_traits_num{6,1} '-' all_traits_num{cor_mat(i,2),1}],R(6,cor_mat(i,2)), P(6,cor_mat(i,2))};
        anc_corr_cell(i,1) = str2double(all_traits_num{6,1});
        anc_corr_cell(i,2) = str2double(all_traits_num{cor_mat(i,2),1});
        anc_corr_cell(i,3) = R(6,cor_mat(i,2));
        anc_corr_cell(i,4) = P(6,cor_mat(i,2));
        abs_anc_corr_cell(i,:) = {[all_traits_num{6,1} '-' all_traits_num{cor_mat(i,2),1}],abs(R(6,cor_mat(i,2))), P(6,cor_mat(i,2))};
    
    elseif (cor_mat(i,1) == 7) %trait 7 vs the world

        %anc_corr_cell(i,:) = {[all_traits_num{7,1} '-' all_traits_num{cor_mat(i,2),1}],R(7,cor_mat(i,2)), P(7,cor_mat(i,2))};
        anc_corr_cell(i,1) = str2double(all_traits_num{7,1});
        anc_corr_cell(i,2) = str2double(all_traits_num{cor_mat(i,2),1});
        anc_corr_cell(i,3) = R(7,cor_mat(i,2));
        anc_corr_cell(i,4) = P(7,cor_mat(i,2));
        abs_anc_corr_cell(i,:) = {[all_traits_num{7,1} '-' all_traits_num{cor_mat(i,2),1}],abs(R(7,cor_mat(i,2))), P(7,cor_mat(i,2))};
        
    end
    
end

abs_anc_corr_cell_sorted = sortrows(abs_anc_corr_cell, 2, 'descend');
anc_corr_cell_sorted = sortrows(anc_corr_cell, 4, 'ascend');