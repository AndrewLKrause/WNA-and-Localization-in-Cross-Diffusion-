function r = sG(p,u)
    u1 = u(1:p.np); % identify solution component 1
    u2 = u(p.np + 1:2*p.np); % identify solution component 2
    u3 = u(2*p.np + 1:3*p.np); % identify solution component 3
    u4 = u(3*p.np + 1:4*p.np); % identify solution component 4
    par = u(p.nu + 1:end); % identify parameters
    a = par(1);
    b = par(2);
    c = par(3);
    D = par(4);
    r_par = par(5);
    r1 = par(6);
    r2 = par(7);
    f1 = r_par*u1.*(1 - u1); % non-linearity for u1 
    f2 = u1 - a*u2 - r2*u2.*u3 + b*(u1.^2).*u4; % non-linearity for u2
    f3 = - r2*u2.*u3 + (r1 + b*u1.^2).*u4 + u1.*(1 - u3);
    f4 = r2*u2.*u3 - (r1 + b*u1.^2).*u4;
    f = [f1; f2; f3; f4];
    [Kaux, ~, ~] = p.pdeo.fem.assema(p.pdeo.grid, c*u1./(1 + u1.^2), 0, 0);
    r = kron(diag([1, D, D, D]), p.mat.K)*u(1:p.nu) -...
        kron([0 1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0], Kaux)*u(1:p.nu) - p.mat.M*f; % calculation of the residual
end