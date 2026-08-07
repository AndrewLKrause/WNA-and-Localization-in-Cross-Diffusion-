%% 1 - initialising the problem
changedir = 1;

if changedir==1
    sleep()
end

    close all;
    clc
    keep pphome;
    p = [];
    % r1, r2, a1, a2, b1, b2, d1, d2, d11, d12, d21, d22
    par = [1.0, 1.0, 0.9, 0.5, 0.6, 0.2, 1.0, 1.0, 0.0, 0.0, 52.86845318178516, 0.0]; 
    h = 1e-2;
    p = SKTinit(p, h, par); % also use nper=80, 120, 160, 200
    p.sw.jac = 0;

%% 2 - continue trivial branch to find BP 
    tic;
    p = cont(p, 25);
    toc

    % bif point: a2 = 0.554362

    [Gua, Gun] = jaccheck(p);

    max(max(abs(Gun - Gua)));

    % pause

%% 3 - switch to periodic branch and continue. For comparison of \ and 
    % lssbel, switch off stuff not related to lss 
    p = swibra('p', 'bpt1', 'b1', 1e-2);
    p.nc.dsmin = 1e-10;
    p.sw.spcalc = 1;
    p.sw.foldcheck = 1;
    p.sw.bifcheck = 2;
    p.sw.verb = 2;
    p0 = p;
    t1 = tic;
    p = cont(p, 150);
    t1 = toc(t1); % cont with default settings

    p = swibra('b1', 'bpt1', 'snake', 1e-2);
    p.nc.dsmin = 1e-10;
    p.sw.verb = 2;
    p0 = p;
    t1 = tic;
    p = cont(p, 1000);
    t1 = toc(t1); % cont with default settings

    % p = swibra('p', 'bpt46', 'b46', 0.02);
    % p.nc.dsmax = 0.15;
    % p.sw.jac = 0;
    % p.sw.spcalc = 1;
    % p.sw.foldcheck = 1;
    % p.sw.bifcheck = 2; 
    % p.sw.verb = 2;
    % p0 = p;
    % t1 = tic;
    % p = cont(p, 43);
    % t1 = toc(t1); % cont with default settings
    % 
    % p = swibra('p', 'bpt1', 'baux', 0.02);
    % p.nc.dsmax = 0.15; 
    % p.sw.jac = 0;
    % p.sw.spcalc = 1;
    % p.sw.foldcheck = 1;
    % p.sw.bifcheck = 2; 
    % p.sw.verb = 2;
    % p0 = p;
    % t1 = tic;
    % p = cont(p, 81);
    % t1 = toc(t1); % cont with default settings
    % 
    % p = swibra('p', 'bpt3', 'b3', 0.02);
    % p.nc.dsmax = 0.15; 
    % p.sw.spcalc = 1;
    % p.sw.foldcheck = 1;
    % p.sw.bifcheck = 2; 
    % p.sw.verb = 2;
    % p0 = p;
    % t1 = tic;
    % p = cont(p, 252);
    % t1 = toc(t1); % cont with default settings
    % 
    % p = swibra('p', 'bpt4', 'b4', 0.02);
    % p.nc.dsmax = 0.15;
    % p.sw.spcalc = 1;
    % p.sw.foldcheck = 1;
    % p.sw.bifcheck = 2;
    % p.sw.verb = 2;
    % p0 = p;
    % t1 = tic;
    % p = cont(p, 257);
    % t1 = toc(t1); % cont with default settings
    % 
    % p = swibra('p', 'bpt5', 'b5', 0.02);
    % p.nc.dsmax = 0.15; 
    % p.sw.spcalc = 1;
    % p.sw.foldcheck = 1;
    % p.sw.bifcheck = 2;
    % p.sw.verb = 2;
    % p0 = p;
    % t1 = tic;
    % p = cont(p, 266);
    % t1 = toc(t1); % cont with default settings
    % 
    % p = swibra('p', 'bpt6', 'b6', 0.02);
    % p.nc.dsmax = 0.15; 
    % p.sw.spcalc = 1;
    % p.sw.foldcheck = 1;
    % p.sw.bifcheck = 2;
    % p.sw.verb = 2;
    % p0 = p;
    % t1 = tic;
    % p = cont(p, 276);
    % t1 = toc(t1); % cont with default settings
    % 
    % p = swibra('p', 'bpt7', 'b7', 0.01);
    % p.nc.dsmax = 0.15;
    % p.sw.spcalc = 1;
    % p.sw.foldcheck = 1;
    % p.sw.bifcheck = 2;
    % p.sw.verb = 2;
    % p0 = p;
    % t1 = tic;
    % p = cont(p, 280);
    % t1 = toc(t1); % cont with default settings
    % 
    % p = swibra('p', 'bpt8', 'b8', 0.02);
    % p.nc.dsmax = 0.15; 
    % p.sw.spcalc = 1;
    % p.sw.foldcheck = 1;
    % p.sw.bifcheck = 2;
    % p.sw.verb = 2;
    % p0 = p;
    % t1 = tic;
    % p = cont(p, 291);
    % t1 = toc(t1); % cont with default settings
    % 
    % p = swibra('p', 'bpt9', 'b9', 0.02);
    % p.nc.dsmax = 0.15; 
    % p.sw.spcalc = 1;
    % p.sw.foldcheck = 1;
    % p.sw.bifcheck = 2;
    % p.sw.verb = 2;
    % p0 = p;
    % t1 = tic;
    % p = cont(p, 296);
    % t1 = toc(t1); % cont with default settings


    % p = p0;
    % bw = 0;
    % beltol = 1e-4;
    % belmaxit = 5;
    % p = setbel(p, bw, beltol, belmaxit, @lss); % lssbel 
    % t2 = tic;
    % p = cont(p, 50);
    % t2 = toc(t2); 
    % fprintf('t1 = %g, t2 = %g\n', t1, t2);
    % plotsol(p, 1, 1, 1); 

%% plots    
    lw = 4;
    lwst = 4;
    lwun = 4;
    
    ps = 100;
    fs = 40;
    figure(3)
    plotbra('p', 'pt351', 3, 0, 'tyun', ':k', 'tyst', '-k', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    plotbra('b46', 'pt43', 3, 0, 'tyun', ':r', 'tyst', '-r', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    % plot(L, L2norm, 'LineWidth', lw, 'Color', [0.93, 0.69, 0.13])
    xlabel('$\beta$', 'Interpreter', 'latex')
    ylabel('$||\nabla u||_{L^2}$', 'Interpreter', 'latex')
    set(gca, 'fontsize', fs)