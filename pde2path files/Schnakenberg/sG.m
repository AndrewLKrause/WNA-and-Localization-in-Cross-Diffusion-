function r = sG(p,u)
    u1 = u(1:p.np); % identify solution component 1
    u2 = u(p.np + 1:2*p.np); % identify solution component 2
    par = u(p.nu + 1:end); % identify parameters
    lambda = par(1);
    sigma = par(2);
    d = par(3);
    pp = par(4);
    f1 = - u1 + (u1.^2).*u2 + sigma*(u1 - 1./u2).^2; % non-linearity for u1 
    f2 = lambda - (u1.^2).*u2 - sigma*(u1 - 1./u2).^2; % non-linearity for u2
    f = [f1; f2];
    [Kaux, ~, ~] = p.pdeo.fem.assema(p.pdeo.grid, u1./(1 + u1.^2), 1, 1);
    r = kron([[1, 0]; [0, d]], p.mat.K)*u(1:p.nu) + kron([[0, pp]; [0, 0]], Kaux)*u(1:p.nu) - p.mat.M*f; % calculation of the residual
end