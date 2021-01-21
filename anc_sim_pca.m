function [Sr_a_sim, PoV_sim_a, pop_anc_pca] = anc_sim_pca(pop, all_traits_num, sim_pca_strains)

figure(2)
pop_pca = pop; %make a copy for pca
pop_pca( :, ~any(pop_pca,1) ) = [];  %delete columns with all zeros
pop_pca = pop_pca(:,1:length(all_traits_num));
[N,M]=size(pop_pca);
R=pop_pca'*pop_pca/(length(pop_pca)-1);   % same as cov(X)
[ V,Lambda ] = eigsort( R ); %V contains eigenvectors and Lambda eigenvalues
%the total variance is the trace (sum of the diagonal elements)
PoV_sim_a=100*diag(Lambda)/trace(Lambda); % percent of variance
Ar=V*Lambda.^0.5;   % Factor loading // each column corresponds to the factor loadings of a princicpal component on each of the original variables
%plot factor loadings vs the original variables
Ar_a_sim = Ar;

Sr_a_sim=pop_pca*Ar_a_sim;        % Ancestral Factor Score

%figure(2);
scatter(Sr_a_sim(1:5,1), Sr_a_sim(1:5,2), 50, 'filled', 'g')
%scatter(Sr_a_sim(:,1), Sr_a_sim(:,2), 50, 'filled', 'g')
dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
%text(Sr_a_sim(:,1)+dx, Sr_a_sim(:,2)+dy, sim_pca_strains)
%text(Sr_a_sim(1:5,1)+dx, Sr_a_sim(1:5,2)+dy, sim_pca_strains{1:5})
hold on

%for z=1:length(Ar)
    %plot([0 Ar_a_sim(z,1)]*5, [0 Ar_a_sim(z,2)]*5, 'k');
%    plot([0 Ar_a_sim(z,1)], [0 Ar_a_sim(z,2)], 'k');
%    text(Ar_a_sim(z,1), Ar_a_sim(z,2), all_traits_num{z})
%end
xlabel(['PC 1 (' num2str(PoV_sim_a(1)) '%)'])
ylabel(['PC 2 (' num2str(PoV_sim_a(2)) '%)'])
%axis([-110 110 -110 110])
tmp_max_Sr_a_sim = max(Sr_a_sim);
max_Sr_a_sim = abs(max(tmp_max_Sr_a_sim));
%axis([-max_Sr_a_sim max_Sr_a_sim -max_Sr_a_sim max_Sr_a_sim])
title(['ancestral_simulated traits_PCA'], 'Interpreter', 'none')

filename = ['PCA_simulated_ancestral_traits'];
hold off
pause(.1)
print('-f2', '-dpdf', filename, '-r0')

pop_anc_pca = pop_pca;