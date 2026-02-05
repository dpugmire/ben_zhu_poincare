% prepare and calcuate apar and its derivatives for field-line tracing
% for shift-circular configuration
%
% note that psi is defined at CELL_YLOW, while B field is at CELL_CENTER,
% so interpolation/shifting along y is needed UNLESS true_apar is true
function [apar,dapardx,dapardy,dapardz] = ...
    get_apar_sc(psi,bxy,psixy,zs,sa,sinty,dy0,dz, ...
    zperiod,interp_opt,deriv_opt,true_apar)

    [nx,ny,nz]=size(psi);
    
    apar=zeros(nx,ny,nz);    dapardx=zeros(nx,ny,nz);
    dapardy=zeros(nx,ny,nz); dapardz=zeros(nx,ny,nz);
    ny_cfr=ny+1; iy_cfr=[1:ny];
    zarray = (0:nz-1)*dz; zarrayp=(0:nz)*dz; zmax = dz*nz;
    
    if (~true_apar) % if not true_apar but psi as in most BOUT++ output

    % step 1, interpolation of psi back to to CELL_CENTER
    psis=psi;
    % shift first point across the branch-cut -- hence becomes ny+1 cell
    for i=1:nx
        zarray_shift=mod(zarray+sa(i),zmax);
        psis_tmp=squeeze(psi(i,1,:)); psis_tmp(end+1)=psis_tmp(1);
        psisp1(i,1,:)=spline(zarrayp,psis_tmp,zarray_shift);
    end
    
    if (interp_opt == 0)
    % linear interpolation
        fprintf('Linear interpolation of psi back to cell center.\n');
        for i=1:nx
            % closed flux surface region
            psi(i,iy_cfr(1:end-1),:)=0.5*(psis(i,iy_cfr(1:end-1),:)+psis(i,iy_cfr(2:end),:));
            psi(i,ny,:)=0.5*(psisp1(i,1,:)+psis(i,ny,:));
        end
    
    elseif (interp_opt == 1)
    % cubic spline interpolation
        fprintf('Cubic spline interpolation of psi back to cell center.\n');
        for i=1:nx
            for k=1:nz
            % closed flux surface region
                psis_cfr=squeeze(psis(i,iy_cfr,k),k);
                psis_cfr(end+1)=psisp1(i,1,k);
                psi(i,1:ny,k)=spline([1:ny_cfr],psis_cfr,[1:ny_cfr-1]+0.5);
            end
        end

    else
        fprintf('Unknown interpolation method for psi!\n');
    end

    % step2, apar=psi*bxy unless psi is truly apar
	for k=1:nz
        apar(:,:,k)=psi(:,:,k).*bxy;
    end

    else % if the input is a true apar in the cell center
        apar=psi;
    end

    % step 3, get apar derivatives
    
    % first shift to flux coordinate for d/dpsi
    apars=zeros(nx,ny,nz); dapardpsi=zeros(nx,ny,nz);
    kz=[0:nz/2 -nz/2+1:1:-1]*zperiod; kz=kz';
    ci=sqrt(-1);
    for i=1:nx
        for j=1:ny
            apars(i,j,:)=real(ifft(fft(squeeze(apar(i,j,:))).*exp(-ci*zs(i,j)*kz)));
        end
    end
    
    if (deriv_opt == 0)
    % numerical differentiate, less accurate
        fprintf('Central finite difference for Apar.\n');
        dpsi=psixy(101,39)-psixy(100,39);
        dapardpsi(2:nx-1,:,:)=0.5*(apars(3:nx,:,:)-apars(1:nx-2,:,:))/dpsi;

        dapardz(:,:,1)=0.5*(apar(:,:,2)-apar(:,:,nz))/dz;
        dapardz(:,:,2:nz-1)=0.5*(apar(:,:,3:nz)-apar(:,:,1:nz-2))/dz;
        dapardz(:,:,nz)=0.5*(apar(:,:,1)-apar(:,:,nz-1))/dz;
    
        for k=1:nz
            dapardx(:,:,k)=dapardpsi(:,:,k)+sinty.*dapardz(:,:,k);
        end
    
        % d/dy is a little tricky
        for i=1:nx
            dapardy(i,iy_cfr(2:end-1),:)=0.5*(apar(i,iy_cfr(3:end),:)-apar(i,iy_cfr(1:end-2),:))/dy0;
            % reverse shift last point to be added before first point -- becomes -1 cell
            zarray_shift=mod(zarray-sa(i),zmax);
            apar_tmp=squeeze(apar(i,ny,:)); apar_tmp(end+1)=apar_tmp(1);
            aparm1(1,1,:)=spline(zarrayp,apar_tmp,zarray_shift);
            dapardy(i,5,:)=0.5*(apar(i,2,:)-aparm1)/dy0;
            % then forward shift first point to be added to last point -- becomes ny+1 cell
            zarray_shift=mod(zarray+sa(i),zmax);
            apar_tmp=squeeze(apar(i,1,:)); apar_tmp(end+1)=apar_tmp(1);
            aparp1(1,1,:)=spline(zarrayp,apar_tmp,zarray_shift);           
            dapardy(i,ny,:)=0.5*(aparp1-apar(i,ny-1,:))/dy0;
        end
    
        % or 4th order cfd (to be implemented)
        % dapardx(3:nx-2,:,:)=(-apar(5:nx,:,:)+8.*apar(4:nx-1,:,:)-8.*apar(2:nx-3)+apar(1:nx-4,:,:))/12./dx;
    
    elseif (deriv_opt == 1)
    % analytical differentiate based on cubic spline interpolation that
    % ensures first and second order derivatives are continous across the
    % cell boundaries
        fprintf('Cubic spline fit then differentiation of Apar.\n');
    
        for k=1:nz
            for j=1:ny
                apar_tmp=apars(:,j,k);
                apar_tmp(end+1)=apar_tmp(end); % pathing last point
                xarray=psixy(:,j);xarray(end+1)=2.*xarray(end)-xarray(end-1);
                pp=spline(xarray,apar_tmp);
                dapardpsi(:,j,k)=pp.coefs(:,3);
            end
        end    
    
        % d/dy again is a little tricky
        yiarray_sol=double(1:ny+1);yiarray_cfr=double(1:ny_cfr);
        for k=1:nz
            for i=1:nx
                % closed flux surface region
                % forward shift first point to be added to last point -- becomes ny+1 cell
                zarray_shift=mod(zarray(k)+sa(i),zmax);
                apar_tmp=squeeze(apar(i,1,:)); apar_tmp(end+1)=apar_tmp(1);
                aparp1=spline(zarrayp,apar_tmp,zarray_shift);
                apar_tmp=squeeze(apar(i,iy_cfr,k)); apar_tmp(end+1)=aparp1; % reuse apar_tmp, caution
                pp=spline(yiarray_cfr,apar_tmp);
                dapardy(i,iy_cfr,k)=pp.coefs(:,3);
            end
        end
        dapardy=dapardy/dy0;
    
        for j=1:ny
            for i=1:nx
                apar_tmp=squeeze(apar(i,j,:));
                apar_tmp(end+1)=apar_tmp(1);
                pp=spline(zarrayp,apar_tmp);
                dapardz(i,j,:)=pp.coefs(:,3);
            end
        end
    
        for k=1:nz
            dapardx(:,:,k)=dapardpsi(:,:,k)+sinty.*dapardz(:,:,k);
        end
    
    end
    
end
