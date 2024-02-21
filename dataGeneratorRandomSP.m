function [Z, mu, v] = dataGeneratorRandomSP(Z, sigma, ind_in, seed, type)
r = 0.5; q = pi * r/2; b = 1;

rng(seed)
xx = Z(:,1); yy = Z(:,2); zz = Z(:,3);
n = length(xx);
mu = NaN(n, 1);
a = zeros(n, 1); d = zeros(n, 1);
ind_all = 1:n;
ind_out = setdiff(ind_all, ind_in);

% Part 1.1
ind = (xx + 3.5 >= 0) & (yy > 0) & (xx + 3.5 < 3.5);
a(ind) = q + xx(ind) + 3.5;
d(ind) = yy(ind) - r;
% Part 1.2
ind = (xx + 3.5 >= 3.5) & (yy > 0) & (xx + 3.5 < 7);
a(ind) = q + abs((xx(ind) + 3.5) - 7);
d(ind) = yy(ind) - r;

% Part 2.1
ind = (xx + 3.5 >= 0) & (yy <= 0) & (xx + 3.5 < 3.5);
a(ind) = -q - (xx(ind) + 3.5);
d(ind) = -r - yy(ind);
% Part 2.2
ind = (xx + 3.5 >= 3.5) & (yy <= 0) & (xx + 3.5 < 7);
a(ind) = -q - abs((xx(ind) + 3.5) - 7);
d(ind) = -r - yy(ind);

% Part 3.1
ind = (xx + 3.5 < 0);
a(ind) = -atan(yy(ind)./(xx(ind) + 3.5)) * r;
d(ind) = sqrt((xx(ind) + 3.5).^2 + yy(ind).^2) - r;
% Part 3.2
ind = (xx + 3.5 >= 7);
a(ind) = atan(yy(ind)./abs(xx(ind) + 3.5 - 7)) * r;
d(ind) = sqrt((xx(ind) + 3.5 - 7).^2 + yy(ind).^2) - r;

switch type
    case 1
        mu = (a * b + d.^2).* (2 - zz.^2);
    case 2
        mu = (a.^3 * 0.25 + d.^2).* cos(pi * (zz + 0.5));
end

% exclude
ind = (yy > -0.1) & (yy < 0.1) & (xx >= -3.5) & (xx <= 3.5);
mu(ind == 1) = nan;
mu(ind_out) = nan;

e = normrnd(0, sigma, n, 1);
v = mu + e;
