function Gu=sGjac(p,u)
    u1 = u(1:p.np); % identify solution component 1
    u2 = u(p.np + 1:2*p.np); % identify solution component 2
    par = u(p.nu + 1:end); % identify parameters
    d = par(3);
    pp = par(4);
    n = p.np;
    [f1u, f1v, f2u, f2v] = njac(p, u, par); % loading the non-linear jacobian in nodal form, see below
    Fu = [[spdiags(f1u, 0, n, n), spdiags(f1v, 0, n, n)];
          [spdiags(f2u, 0, n, n), spdiags(f2v, 0, n, n)]]; % nonlinear part of the jacobian
    [Kaux, ~, ~] = p.pdeo.fem.assema(p.pdeo.grid, u1./(1 + u1.^2), 1, 1);
    % x = p.pdeo.grid.p(1, :)'; % Node coordinates (n x 1)
    % t = p.pdeo.grid.t; % 1D elements connectivity (3 x num_elements or 2 x num_elements)
    % 
    % node1 = t(1, :);
    % node2 = t(2, :);
    % 
    % % Element-wise u1 average, d'(u1), and du2/dx
    % u1_elem = 0.5 * (u1(node1) + u1(node2));
    % dprime_elem = (1 - u1_elem.^2) ./ ((1 + u1_elem.^2).^2);
    % du2_elem = (u2(node2) - u2(node1)); % du2 = (u2_right - u2_left)
    % 
    % weight = 0.5 * dprime_elem .* du2_elem; % element weights
    % 
    % % Triplets for sparse matrix assembly (2x2 per element)
    % ii = [node1, node1, node2, node2];
    % jj = [node1, node2, node1, node2];
    % vv = [- weight, -weight, weight, weight];
    % 
    % Kconv = sparse(ii, jj, vv, n, n);
    % Gu = kron([[1, 0]; [0, d]], p.mat.K) + kron([[pp, 0]; [0, 0]], Kconv) + kron([[0, pp]; [0, 0]], Kaux) - p.mat.M*Fu; % assemble the jacobian
    Gu = kron([[1, 0]; [0, d]], p.mat.K) + kron([[0, pp]; [0, 0]], Kaux) - p.mat.M*Fu;
end

function [f1u, f1v, f2u, f2v] = njac(p, u, par) % Jacobian for Schnakenberg, nodal version
    u1 = u(1:p.np); % identify solution component 1
    u2 = u(p.np + 1:2*p.np); % identify solution component 2
    par = u(p.nu + 1:end);
    sigma = par(2);
    % entries of the jacobian
    % f1 = - u1 + (u1.^2).*u2 + sigma*(u1 - 1./u2).^2; % non-linearity for u1 
    % f2 = lambda - (u1.^2).*u2 - sigma*(u1 - 1./u2).^2; % non-linearity for u2
    f1u = - 1 + 2*u1.*u2 + 2*sigma*(u1 - 1./u2);
    f1v = u1.^2 + 2*sigma*(u1 - 1./u2).*(1./u2.^2); 
    f2u = - 2*u1.*u2 - 2*sigma*(u1 - 1./u2); 
    f2v = - u1.^2 - 2*sigma*(u1 - 1./u2).*(1./u2.^2);
end