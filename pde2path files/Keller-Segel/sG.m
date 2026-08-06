function r = sG(p,u)
    u1 = u(1:p.np); % identify solution component 1
    u2 = u(p.np + 1:2*p.np); % identify solution component 2
    par = u(p.nu + 1:end); % identify parameters
    a = par(1);
    b = par(2);
    chi = par(3);
    r = par(4);
    K = par(5);
    D1 = par(6);
    D2 = par(7);
    f1 = r*u1.*(1 - u1/K); % non-linearity for u1 
    f2 = a*u1 - b*u2; % non-linearity for u2
    f = [f1; f2];
    [Kaux, ~, ~] = p.pdeo.fem.assema(p.pdeo.grid, u1, 1, 1);
    r = kron([[D1, 0]; [0, D2]], p.mat.K)*u(1:p.nu) + kron([[0, - chi]; [0, 0]], Kaux)*u(1:p.nu) - p.mat.M*f; % calculation of the residual
end