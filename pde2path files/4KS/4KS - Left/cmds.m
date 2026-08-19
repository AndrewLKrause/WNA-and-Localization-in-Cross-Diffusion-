%% 1 - initialising the problem
changedir = 1;

if changedir==1
    sleep()
end

    close all;
    clc
    keep pphome;
    p = [];
    % a, b, c, D, r, r1, r2
    par = [0.9, 0.0, 4.4, 1.0, 0.12, 0.2, 0.5];
    % par = [0.9, 0.8, 2.9, 1.0, 0.12, 0.2, 0.5];
    h = 2e-2;
    p = fourKSinit(p, h, par); % also use nper=80, 120, 160, 200
    p.sw.jac = 0;
    p.nc.dsmin = 1e-20;
    p.nc.dsmax = 5e-3;

%% 2 - continue trivial branch to find BP 
    tic;
    p = findbif(p, 7);
    toc

    % bif point 1: c = 4.29791
    % bif point 2: c = 2.90548

%% 3 - switch to periodic branch and continue. For comparison of \ and 
    % lssbel, switch off stuff not related to lss 
    istep = - 1e-2;

    contset = 1:7;

    for counter=contset
        p = swibra('p', char(strcat('bpt', string(counter))), char(strcat('b', string(counter))), istep);
        p.nc.dsmax = 1e-2;
        p.nc.dsmin = 1e-20;
        p.sw.spcalc = 1;
        p.sw.foldcheck = 1;
        p.sw.bifcheck = 2;
        p.sw.verb = 2;
        p0 = p;
        t1 = tic;
        p = cont(p, 50);
        t1 = toc(t1); % cont with default settings
    end

    p = swibra('b1', 'bpt1', 'snake', - 1e-3);
    p.nc.dsmax = 2e-3;
    p.nc.dsmin = 1e-20;
    p.sw.verb = 2;
    p0 = p;
    t1 = tic;
    p = cont(p, 335);
    t1 = toc(t1); % cont with default settings

%% plots    
    lw = 4;
    lwst = 4;
    lwun = 4;

    contset = 1:7;
    
    ps = 100;
    fs = 40;
    figure(3)
    plotbra('p', 3, 0, 'tyun', ':k', 'tyst', '-k', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    for counter=contset
        plotbra(char(strcat('b', string(counter))), 3, 0, 'tyun', '-b', 'tyst', '-r', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    end
    plotbra('snake', 3, 0, 'tyun', '-b', 'tyst', '-r', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    
    xlabel('$b$', 'Interpreter', 'latex')
    ylabel('$||u||_{L^2}$', 'Interpreter', 'latex')
    set(gca, 'fontsize', fs)
    load('snake/pt5.mat')
    mesh = getpte(p);
    length = 2*mesh(end);
    yt = yticks;
    yticklabels(yt*sqrt(length));