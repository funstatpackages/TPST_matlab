function [Z, mu, v, sigma] = popGeneratorHS(Z, psnr, ind_in, seed, type)
r0 = 0.1; r = 0.5; l = 3; b = 1; 
q = pi * r/2;

rng(seed)
xx = Z(:,1); yy = Z(:,2); zz = Z(:,3);
n = length(xx);
mu = NaN(n, 1);
a = zeros(n, 1); d = zeros(n, 1);
ind_all = 1:n;
ind_out = setdiff(ind_all, ind_in);

% Part 1
ind = (xx >= 0) & (yy > 0);
a(ind) = q + xx(ind);
d(ind) = yy(ind) - r;
    
% Part 2
ind = (xx >= 0) & (yy <= 0);
a(ind) = -q - xx(ind);
d(ind) = -r - yy(ind);

% Part 3
ind = xx < 0;
a(ind) = -atan(yy(ind)./xx(ind)) * r;
d(ind) = sqrt(xx(ind).^2 + yy(ind).^2) - r;

ind = (abs(d) > r - r0) | (xx > l & (xx - l).^2 + d.^2 > (r - r0)^2);
% figure
% scatter3(xx(ind == 1), yy(ind == 1), zz(ind == 1));
% view(40,35)

if type == 1
    mu = (a * b + d.^2) .* (2 - zz.^2);
end

if type == 2
    mu = (a * b + d.^2) .* sin(pi * zz);
end

if type == 3
    mu = (a * b + d.^2) .* cos(pi * zz);
end

if type == 4
    mu = (a * b + d.^2) .* exp(- 8 * zz.^2);
end

mu(ind == 1) = nan;
mu(ind_out) = nan;
% figure
% scatter3(xx(ind_out), yy(ind_out), zz(ind_out));
% view(40,35)

mu_max = nanmax(abs(mu));
sigma = 10^((20 * log10(mu_max) - psnr)/20);
e = normrnd(0, sigma, n, 1);
v = mu + e;