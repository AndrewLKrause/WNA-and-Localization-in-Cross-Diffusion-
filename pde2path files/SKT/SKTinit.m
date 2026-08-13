function p=schnakinit(p, h, par)
    p = stanparam(p);
    screenlayout(p);
    p = setfn(p, 'p');
    p.nc.lammin = 0.0;
    p.nc.lammax = 3.0;
    p.nc.bisecmax = 10;
    p.nc.neq = 2;
    p.sw.sfem = - 1;
    p.fuha.sG = @sG;
    p.fuha.sGjac = @sGjac; 
    lx = 4*8.224672274330526; % wavenumer of the critical mode
    p.pdeo = stanpdeo1D(lx, h);
    p.np = p.pdeo.grid.nPoints;
    p.nu = p.np*p.nc.neq; 
    p = setfemops(p);
    p.nc.ilam = 4;
    p.sol.xi = 1/p.nu;
    p.sol.ds = 1e-2;
    p.nc.dsmax = 1e-2; 
    u = (par(4) - par(5))/(par(3)*par(4) - par(5)*par(6))*ones(p.np, 1);
    v = (par(3) - par(6))/(par(3)*par(4) - par(5)*par(6))*ones(p.np, 1); % hom.soln 
    p.u = [u; v; par'];
    p.nc.nsteps = 200;
    p.plot.pmod = 1;
    p.file.smod = 5;
end