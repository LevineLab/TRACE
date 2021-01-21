function [Ar_e_sim, PoV_sim_e, Sr_e_sim, pop_evo_pca] = evo_sim_pca(pop_evo, all_traits_num, sim_pca_strains)

figure(3)
pop_evo_pca = pop_evo; %make a copy for pca
pop_evo_pca( :, ~any(pop_evo_pca,1) ) = [];  %delete columns with all zeros
pop_evo_pca = pop_evo_pca(:,1:length(all_traits_num));
[N,M]=size(pop_evo_pca);
R=pop_evo_pca'*pop_evo_pca/(length(pop_evo_pca)-1);   % same as cov(X)
[ V,Lambda ] = eigsort( R ); %V contains eigenvectors and Lambda eigenvalues
%the total variance is the trace (sum of the diagonal elements)
PoV_sim_e=100*diag(Lambda)/trace(Lambda); % percent of variance
Ar=V*Lambda.^0.5;   % Factor loading // each column corresponds to the factor loadings of a princicpal component on each of the original variables
%plot factor loadings vs the original variables
Ar_e_sim = Ar;

Sr_e_sim=pop_evo_pca*Ar_e_sim;        % Ancestral Factor Score

%figure(2);
%scatter(Sr_e_sim(:,1), Sr_e_sim(:,2), 50, 'filled', 'g')
scatter(Sr_e_sim(1:5,1), Sr_e_sim(1:5,2), 50, 'filled', 'g')
dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
%text(Sr_e_sim(:,1)+dx, Sr_e_sim(:,2)+dy, sim_pca_strains)
hold on

%for z=1:length(Ar_e_sim)
    %plot([0 Ar_e_sim(z,1)]*5, [0 Ar_e_sim(z,2)]*5, 'k');
%    plot([0 Ar_e_sim(z,1)], [0 Ar_e_sim(z,2)], 'k');
%    text(Ar_e_sim(z,1), Ar_e_sim(z,2), all_traits_num{z})
%end
xlabel(['PC 1 (' num2str(PoV_sim_e(1)) '%)'])
ylabel(['PC 2 (' num2str(PoV_sim_e(2)) '%)'])
%axis([-110 110 -110 110])
tmp_max_Sr_e_sim = max(Sr_e_sim);
max_Sr_e_sim = abs(max(tmp_max_Sr_e_sim));
%axis([-max_Sr_e_sim max_Sr_e_sim -max_Sr_e_sim max_Sr_e_sim])
title(['evolved_simulated traits_PCA'], 'Interpreter', 'none')

filename = ['PCA_simulated_evolved_traits'];
hold off
pause(.1)
print('-f3', '-dpdf', filename, '-r0')