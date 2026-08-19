%% 1 - initialising the problem
changedir = 1;

if changedir==1
    sleep()
end

    close all;
    clc
    keep pphome;
    p = [];
    % a, b, chi, r, K, D1, D2
    par = [1.0, 20.08, 23.0, 1.0, 1.0, 1.0, 0.1]; 
    h = 1e-2;
    p = KSinit(p, h, par); % also use nper=80, 120, 160, 200
    p.nc.dsmax = 5e-4;
    p.sw.jac = 0;

%% 2 - continue trivial branch to find BP 
    tic;
    p = findbif(p, 20);
    toc;

    % bif point: b = 20.0668

%% 3 - switch to periodic branch and continue. For comparison of \ and 
    % lssbel, switch off stuff not related to lss 
    istep = 1e-2;

    contset = 1:13;

    for counter=contset
        p = swibra('p', char(strcat('bpt', string(counter))), char(strcat('b', string(counter))), istep);
        p.nc.dsmax = 1e-2;
        p.nc.dsmin = 1e-10;
        p.sw.spcalc = 1;
        p.sw.foldcheck = 1;
        p.sw.bifcheck = 2;
        p.sw.verb = 2;
        p0 = p;
        t1 = tic;
        if counter>=11
            p.nc.bisecmax = 20;
        end
        p = cont(p, 160);
        t1 = toc(t1); % cont with default settings
    end

    p = swibra('b1', 'bpt1', 'snake', 1e-2);
    p.nc.dsmin = 1e-10;
    p.sw.verb = 2;
    p0 = p;
    t1 = tic;
    p = cont(p, 560);
    t1 = toc(t1); % cont with default settings

%% plots    
    lw = 4;
    lwst = 4;
    lwun = 4;

    contset = 1:13;
    
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