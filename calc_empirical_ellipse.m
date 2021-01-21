function [evo_emp_ellipse, G, evo_strains] = calc_empirical_ellipse(Sr_e, current_strains_evolved)

Sr_e_pc = Sr_e(:,1:2); %first two PC's
%make a matrix that mirrors current_strains_evolved dimensions
G = [1*ones(6,1);2*ones(6,1);3*ones(6,1);4*ones(5,1);5*ones(6,1)]; %factors of ancestral strains
evo_ellipse = {}; %contains ellipse information of all evolved strains
evo_strains = [125 1691 1952 2931 2937]; %array of all strains
figure(1)
for k=1:5
    
    
    %# indices of points in this group
    idx = ( G == k );
    
    %# substract mean
    Mu = mean( Sr_e_pc(idx,:) );
    X0 = bsxfun(@minus, Sr_e_pc(idx,:), Mu);
    
    %Draw ellipse with scale of 2 standard deviations 
    STD = 2;                     %# 2 standard deviations
    conf = 2*normcdf(STD)-1;     %# covers around 95% of population
    scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions
    
    %# eigen decomposition [sorted by eigen values]
    Cov = cov(X0) * scale;
    [V D] = eig(Cov);
    [D order] = sort(diag(D), 'descend');
    D = diag(D);
    V = V(:, order);

    t = linspace(0,2*pi,100);
    e = [cos(t) ; sin(t)];        %# unit circle
    VV = V*sqrt(D);               %# scale eigenvectors
    e = bsxfun(@plus, VV*e, Mu');
    plot(e(1,:), e(2,:), 'Color','k'); hold on
    gscatter(Sr_e_pc(:,1),Sr_e_pc(:,2), cell2mat(current_strains_evolved)); hold on
    quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
    quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')

    %calculate major/minor axis length https://math.stackexchange.com/questions/185954/finding-the-major-and-minor-axes-of-an-n-dimensional-ellipse
    evo_ellipse(k,1) = {evo_strains(k)};
    evo_ellipse(k,2) = {sqrt(D(1))/2}; %divide by 2 to get semi-major 
    evo_ellipse(k,3) = {sqrt(D(4))/2}; %divide by 2 to get semi-minor
    evo_ellipse(k,4) = {Mu}; %center of ellipse
    
    
    %calculate if point is inside ellipse
    %https://math.stackexchange.com/questions/76457/check-if-a-point-is-within-an-ellipse

end

title(['evolved_empirical_PCA_ellipse'], 'Interpreter', 'none')
filename = ['PCA_evolved_ellipse'];
hold off
%pause(.1)
print('-f1', '-dpdf', filename, '-r0')

evo_emp_ellipse = evo_ellipse;