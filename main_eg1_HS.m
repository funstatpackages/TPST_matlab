clear;
format long g;
shape = 'HorseShoe';

load('data/V1HS.mat'); load('data/Th1HS.mat'); 
load('data/V2HS.mat'); load('data/Th2HS.mat'); 
Vs = cell(2,1); Vs{1} = V1; Vs{2} = V2; 
Ths = cell(2,1); Ths{1} = Th1; Ths{2} = Th2; 

psnr = 10; % 5, 10
n = 50000; % sample size; 20000, 50000
type = 1; % 1, 2
Th_ind = 1; % choice of partition;
nx_pop = 105; ny_pop = 69; nz_pop = 19;% population size;
x0 = -0.88:(3.38+0.88)/(nx_pop-1):3.38;
y0 = -0.88:(0.88+0.88)/(ny_pop-1):0.88;
z0 = -0.48:1/(nz_pop-1):0.48;
[xx_grid, yy_grid, zz_grid] = meshgrid(x0, y0, z0);
xx_pop = xx_grid(:); yy_pop = yy_grid(:); zz_pop = zz_grid(:);
Zpop = [xx_pop yy_pop zz_pop];
Npop = size(Zpop, 1); ind_all = 1:Npop;

VT1 = triangulation(Th1, V1);
[ind1_pop, ~] = pointLocation(VT1, Zpop);
ind1_pop = ind_all(~isnan(ind1_pop));
VT2 = triangulation(Th2, V2);
[ind2_pop, ~] = pointLocation(VT2, Zpop);
ind2_pop = ind_all(~isnan(ind2_pop));
ind_pop = intersect(ind1_pop, ind2_pop);    

[~, mupop, Ypop, sigma] = popGeneratorHS(Zpop, psnr, ind_pop, 2020, type);
pop = [Ypop, mupop, Zpop];

V = Vs{Th_ind}; Th = Ths{Th_ind};
indY_pop = 1:size(pop, 1);
indY_pop = indY_pop(~isnan(Ypop));
ind_pop = intersect(ind_pop, indY_pop);
mup = mupop(ind_pop);

d = 3;
r = 1;
K = energyM3D(V, Th, d);
H = smoothness(V, Th, d, r); % smoothness matrix;
Q2 = qrH(H);
[Bpop, ~, ~, ind_tetr_pop] = basis(V, Th, d, Zpop(ind_pop,:));

% candidates for roughness penalty parameter;
index = -6:1:4;
lambda = 10.^index;

% main iteration;
mise_all = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TPST
for iterSeed = 1:100
    rng(iterSeed)
    
    xx = -0.85 + (3.35+0.85) * rand(2 * n, 1); 
    yy = -0.85 + (0.85+0.85) * rand(2 * n, 1);
    zz = -0.45 + (0.45+0.45) * rand(2 * n, 1);
    Z = [xx yy zz];
    ind1 = find(insideVT(V1, Th1, Z));
    ind2 = find(insideVT(V2, Th2, Z));
    ind_in = intersect(ind1, ind2);
    
    [Z, mu, Y] = dataGeneratorRandomHS(Z, sigma, ind_in, iterSeed, type);
    
    ind3 = find(~isnan(Y));
    ind = intersect(ind_in, ind3);
    ind = ind(1:n);
    Zi = Z(ind,:); Yi = Y(ind); mui = mu(ind);
    
    % Estimation;
    [B, ~, ~, ind_tetr] = basis(V, Th, d, Zi);
    [~,gamma,~,~,lamc,gcv,~] = ...
        TPST_est(B, Q2, K, Yi, lambda);
    
    % Prediction;
    Ypred = Bpop * gamma;
    mise = nanmean((mup - Ypred).^2);
    mise_all = [mise_all, mise];
end

mean(mise_all)
