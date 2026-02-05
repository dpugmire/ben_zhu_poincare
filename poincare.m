%% collect data then generate Poincare plot
%
% last updated by B.Zhu (12/2023)

%% Collect puncture points info and generate Poincare plot
%
%  by B.Zhu (10/2023)

addpath('~/Documents/MATLAB/colormaps/');
addpath('~/Documents/MATLAB/utilities/');

clear all;
close all;

set(groot,'DefaultLineLinewidth',2);
set(groot,'DefaultAxesFontSize',20);
set(groot,'DefaultTextFontSize',20);
set(groot,'DefaultLineMarkerSize',2);

% step 0: equilibrium and grid resolution
g = read_eqdsk('g030306.007850_kin_2','../../../GPEC/30306_7850/');
gridfile =  '../kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc';
nx = 260; ny = 128; nz = 256; zperiod = 1; 
nlines = 256;
direction = 1;

% read in grid info
fid = netcdf.open(gridfile, 'nc_nowrite');
vid = netcdf.inqVarID(fid, 'Rxy');
rxy = netcdf.getVar(fid, vid); rxy = double(rxy); rxy = permute(rxy, [2 1]);
vid = netcdf.inqVarID(fid, 'Zxy');
zxy = netcdf.getVar(fid, vid); zxy = double(zxy); zxy = permute(zxy, [2 1]);
vid = netcdf.inqVarID(fid, 'psixy');  
psixy = netcdf.getVar(fid, vid); psixy = double(psixy); psixy = permute(psixy, [2 1]);
vid = netcdf.inqVarID(fid, 'hthe');  
hthe = netcdf.getVar(fid, vid); hthe = permute(hthe, [2 1]);
vid = netcdf.inqVarID(fid, 'ixseps1');  
nxsep = netcdf.getVar(fid, vid) + 1;
vid = netcdf.inqVarID(fid, 'psi_bndry');
psi_bndry = netcdf.getVar(fid, vid);
vid = netcdf.inqVarID(fid, 'psi_axis');
psi_axis = netcdf.getVar(fid, vid);
psin=(psixy(:,55)-psi_axis)/(psi_bndry-psi_axis);
netcdf.close(fid);

nypf=16;
corex = rxy(1,nypf+1:ny-nypf); corey = zxy(1,nypf+1:ny-nypf);
solx = rxy(1,1:nypf); soly = zxy(1,1:nypf);
solx(nypf+1:2*nypf) = rxy(1,ny-nypf+1:end); soly(nypf+1:2*nypf) = zxy(1,ny-nypf+1:end);
solx(2*nypf:2*nypf+nx-1) = rxy(:,end); soly(2*nypf:2*nypf+nx-1) = zxy(:,end);
solx(2*nypf+nx-1:2*nypf+nx+ny-2) = rxy(end,ny:-1:1); soly(2*nypf+nx-1:2*nypf+nx+ny-2) = zxy(end,ny:-1:1);
solx(2*nypf+nx+ny-2:2*nypf+2*nx+ny-3) = rxy(nx:-1:1,1); soly(2*nypf+nx+ny-2:2*nypf+2*nx+ny-3) = zxy(nx:-1:1,1);
tmp = [1:nypf ny-nypf:-1:nypf+1 ny-nypf+1:ny];
sepx = 0.5*(rxy(nxsep-1,tmp)+rxy(nxsep,tmp)); sepy = 0.5*(zxy(nxsep-1,tmp)+zxy(nxsep,tmp)); 

% step 1: generate poincare plot
cm=jet(nlines);

figure(1)
set(gcf,'Position',[100 100 1600 800])

for iline=1:nlines
    fprintf('\tradial index: %i \n',iline);

    if (direction == 1)
        filename = strcat('./mat_pp/x',num2str(iline),'y55z1_v3lc-01-250p.mat');
    elseif (direction == -1)
        filename = strcat('./mat_pp/x',num2str(iline),'y55z1_v3lc-01-250m.mat');
    end
    
    if isfile(filename)
        load(filename);
        subplot(1,2,1)
        hold on
        plot(v2,v3,'.','color',cm(iline,:),'MarkerSize',2);
        subplot(1,2,2)
        hold on
        plot(v5,v4,'.','color',cm(iline,:),'MarkerSize',2);
    end

end

subplot(1,2,1)
hold on
plot(sepx,sepy,'--k');
plot(g.R_limits,g.Z_limits,'k'); 
axis equal; xlim([1.2,2.3]);
xlabel('R(m)'); ylabel('Z(m)')
subplot(1,2,2)
hold on
plot([psixy(195,55),psixy(195,55)],[0,2],'--k');
xlabel('$\psi$','Interpreter','latex'); ylabel('$\theta/\pi$','Interpreter','latex');