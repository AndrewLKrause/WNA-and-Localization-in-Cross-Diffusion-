function r = sG(p,u)
    u1 = u(1:p.np); % identify solution component 1
    u2 = u(p.np + 1:2*p.np); % identify solution component 2
    par = u(p.nu + 1:end); % identify parameters
    r1 = par(1);
    r2 = par(2);
    a1 = par(3);
    a2 = par(4);
    b1 = par(5);
    b2 = par(6);
    d1 = par(7);
    d2 = par(8);
    d11 = par(9);
    d12 = par(10);
    d21 = par(11);
    d22 = par(12);
    f1 = r1*u1.*(1 - a1*u1 - b1*u2); % non-linearity for u1 
    f2 = r2*u2.*(1 - b2*u1 - a2*u2); % non-linearity for u2
    f = [f1; f2];
    [Kaux1, ~, ~] = p.pdeo.fem.assema(p.pdeo.grid, u1, 1, 1);
    [Kaux2, ~, ~] = p.pdeo.fem.assema(p.pdeo.grid, u2, 1, 1);
    r = kron([[d1, 0]; [0, d2]], p.mat.K)*u(1:p.nu) +...
        (kron([[d11, 0]; [0, 0]], Kaux1) + kron([[d12, 0]; [0, 0]], Kaux2))*u(1:p.nu) +...
        + kron([[0, d12]; [0, 0]], Kaux1)*u(1:p.nu) +...
        kron([[0, 0]; [d21, 0]], Kaux2)*u(1:p.nu) +...
        (kron([[0, 0]; [0, d21]], Kaux1) + kron([[0, 0]; [0, d22]], Kaux2))*u(1:p.nu) -...
        p.mat.M*f; % calculation of the residual
end