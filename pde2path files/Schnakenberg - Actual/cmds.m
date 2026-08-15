%% 1 - initialising the problem
changedir = 1;

if changedir==1
    sleep()
end

    close all;
    clc
    keep pphome;
    p = [];
    % lambda, sigma, d, p
    par = [1.05, - 1.5, 1, 2.1]; 
    h = 1e-2;
    p = schnakinit(p, h, par); % also use nper=80, 120, 160, 200
    p.sw.jac = 0;
    p.nc.bisecmax = 50;
    p.nc.dmax = 1e-3;

%% 2 - continue trivial branch to find BP 
    tic;
    p = findbif(p, 8);
    toc
    % p = loadp('p', 'pt31');
    % p.sw.bifcheck = 2;
    % p.nc.dsmax = 5e-4;
    % p.nc.bisecmax = 1000;
    % p.nc.dsminbis = 1e-20;
    % tic
    % p = findbif(p, 3);
    % toc
    % p = loadp('p', 'pt100');
    % p.sw.bifcheck = 1;
    % p.nc.dsmax = 5e-4;
    % p.nc.nsteps = 2000;
    % tic
    % p = findbif(p, 2);
    % toc

    % bif point: p = 2.20512

    % [Gua, Gun] = jaccheck(p);
    % 
    % max(max(abs(Gun - Gua)));

    % pause

%% 3 - switch to periodic branch and continue. For comparison of \ and 
    % lssbel, switch off stuff not related to lss 
    istep = - 1e-2;

    contset = 1;

    for counter=contset
        p = swibra('p', char(strcat('bpt', string(counter))), char(strcat('b', string(counter))), istep);
        p.nc.dsmax = 1e-2;
        p.nc.dsmin = 1e-10;
        p.nc.bisecmax = 50;
        p.nc.dsminbis = 1e-9;
        p.sw.spcalc = 1;
        p.sw.foldcheck = 1;
        p.sw.bifcheck = 2;
        p.sw.verb = 2;
        p0 = p;
        t1 = tic;
        if counter==3
            p = cont(p, 105);
        else
            p = cont(p, 200);
        end
        t1 = toc(t1); % cont with default settings
    end

    p = swibra('b1', 'bpt1', 'snake', - 1e-3);
    p.nc.dsmin = 1e-10;
    p.nc.dsmax = 1e-2;
    p.nc.bisecmax = 20;
    p.sw.verb = 2;
    p0 = p;
    t1 = tic;
    p = cont(p, 140);
    t1 = toc(t1); % cont with default settings

    % p = swibra('snake', 'bpt1', 'extra', - 1e-2);
    % p.nc.dsmax = 5e-2;
    % p.nc.bisecmax = 20;
    % p.sw.verb = 2;
    % p0 = p;
    % t1 = tic;
    % p = findbif(p, 2);
    % t1 = toc(t1); % cont with default settings
    % p.sol.ds = 1e-2;
    % p = cont(p, 15);
    % 
    % p = loadp('extra', 'pt11');
    % p.branch = [];
    % p.file.count  = 0;
    % p.file.bcount = 0;
    % p.file.fcount = 0;
    % p.file.dir = 'actual1';
    % p.nc.ds = - 1e-2;
    % p = cont(p, 150);
    % p = loadp('extra', 'pt11');
    % p.branch = [];
    % p.file.count  = 0;
    % p.file.bcount = 0;
    % p.file.fcount = 0;
    % p.file.dir = 'actual2';
    % p.nc.ds = - 1e-2;
    % p = cont(p, 150);


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

    contset = 1;
    
    ps = 100;
    fs = 40;
    figure(3)
    plotbra('p', 3, 0, 'tyun', ':k', 'tyst', '-k', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    for counter=contset
        plotbra(char(strcat('b', string(counter))), 3, 0, 'tyun', '-b', 'tyst', '-r', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    end
    plotbra('snake', 3, 0, 'tyun', '-b', 'tyst', '-r', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    
    xlabel('$p$', 'Interpreter', 'latex')
    ylabel('$||u||_{L^2}$', 'Interpreter', 'latex')
    set(gca, 'fontsize', fs)
    load('snake/pt5.mat')
    mesh = getpte(p);
    length = 2*mesh(end);
    yt = yticks;
    yticklabels(yt*sqrt(length));