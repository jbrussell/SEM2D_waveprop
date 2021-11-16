function [ out, OUTmovie ] = run_SEM2D( LX, LY, NELX, NELY, P, Fx, Fy, avd_out, Ff0, Fstr, VS, OUTxseis, OUTyseis)
% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% with absorbing boundary conditions,
% zero initial conditions
% and a time dependent force source,
% in a structured undeformed grid.
%
% Version 1a:	domain = rectangular
%            	medium = general (heterogeneous)
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%
% JBR 4/2020: Add absorbing boundary conditions
%

% is_movieout : [0 or 1] % Save output movie
%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

%**** Set here wavefield output : ****
% avd_out : a=accel, v=vel, d=disp
%**** Set here the parameters of the square box domain and mesh : ****
% LX : x-size of the box
% LY : y-size of the box
% NELX : number of elements in x
% NELY : number of elements in y
% P :  polynomial degree (spatial resolution)
%********
dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element
[iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);
nglob = length(x); % (NELX*(NGLL-1)+1)*(NELY*(NGLL-1)+1)

%**** Set here the parameters of the time solver : ****
% Tduration : duration of output seismograms (seconds)
Tduration = max([x;y])./min(VS(:));
CFL   = sqrt(2)/pi; %0.6; %0.6/2; %0.6;       % stability number (time resolution)

%**** Set here receiver locations : ****
% OUTyseis : x-coordinates of stations
% OUTxseis : y-coordinates of stations
OUTnseis = length(OUTyseis);

%**** Set here properties of medium : ****
% VS : velocity model (km/s)

out = [];
%**** Set here the source location : ****
% Fx : x-coordinates of source
% Fy : y-coordinates of source

%**** Set here source properties : ****
% Fstr : 'ricker', 'gaussian' input source type
% Ff0 : % fundamental frequency of source
Ft0 = 1.5/Ff0; % delay
%------------------------------------------
% STEP 2: INITIALIZATION


[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
wgll2 = wgll * wgll' ;

W     = zeros(NGLL,NGLL,NEL);	% for internal forces
M     = zeros(nglob,1);		% global mass matrix, diagonal
rho   = zeros(NGLL,NGLL);	% density will not be stored
mu    = zeros(NGLL,NGLL);	% shear modulus will not be stored
% vs    = zeros(NGLL,NGLL);	% shear velocity
% vs_mat = zeros(NELX,NELY);
% jbr: initialize global arrays of structural properties for BCs
mu_glob = zeros(nglob,1);
rho_glob = zeros(nglob,1);
vs_glob = zeros(nglob,1);

%********

dt = inf; % will be set later

% For this simple mesh the global-local coordinate map (x,y)->(xi,eta)
% is linear, its jacobian is constant
dx_dxi  = 0.5*dxe;
dy_deta = 0.5*dye;
jac = dx_dxi*dy_deta;
coefint1 = jac/dx_dxi^2 ;
coefint2 = jac/dy_deta^2 ;

% FOR EACH ELEMENT ...
% . set physical properties
% . set mass and stiffness matrices
% . set timestep
for ey=1:NELY, 
for ex=1:NELX, 

  e = (ey-1)*NELX+ex;
  ig = iglob(:,:,e);

%**** Set here the physical properties of the heterogeneous medium : ****
 % can be heterogeneous inside the elements
 % and/or discontinuous across elements 
 % example: a low velocity layer
  rho(:,:) = 2.800; % kg*m^-3
  mu(:,:) = rho*VS(ex,ey).^2;
%   vs_mat(ex,ey) = mu(1);
  mu_glob(ig) = mu;
  rho_glob(ig) = rho;
  vs_glob(ig) = sqrt(mu./rho);
%********

 % Diagonal mass matrix
  M(ig) = M(ig) + wgll2 .*rho *jac;

 % Local contributions to the stiffness matrix K
 %  WX(:,:,e) = wgll2 .* mu *jac/dx_dxi^2;
 %  WY(:,:,e) = wgll2 .* mu *jac/dy_deta^2;
  W(:,:,e) = wgll2 .* mu;

 % The timestep dt is set by the stability condition
 %   dt = CFL*min(dx/vs)
  vs = sqrt(mu./rho); 
  if dxe<dye
    vs = max( vs(1:NGLL-1,:), vs(2:NGLL,:) );
    dx = repmat( diff(xgll)*0.5*dxe ,1,NGLL); 
  else
    vs = max( vs(:,1:NGLL-1), vs(:,2:NGLL) );
    dx = repmat( diff(xgll)'*0.5*dye ,NGLL,1); 
  end
  dtloc = dx./vs;
  dt = min( [dt dtloc(1:end)] );

end
end %... of element loop
dt = CFL*dt;
NT = floor(Tduration/dt);

%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

time = (1:NT)'*dt;

%-- SOURCE TERM: point force, time function = Ricker wavelet
%********
[Fx,Fy,Fig] = FindNearestNode(Fx,Fy,x,y);
% Ff0 = 1/10; %0.3; % fundamental frequency
% Ft0 = 1.5/Ff0; % delay
% source time function (at mid-steps)
Ft = src_timef( time-0.5*dt,Fstr, Ff0,Ft0)*1;

%-- initialize data for output seismograms
%********

[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
OUTavd = zeros(OUTnseis,NT);
OUT_trF = zeros(1,NT);

% Points to sample space for traveltime surface calculation
% % xpts = [0:dx_dxi/2:LX];
% % ypts = [0:dy_deta/2:LY];
xpts = [0:dx_dxi:LX];
ypts = [0:dy_deta:LY];
[ypts_mesh,xpts_mesh] = meshgrid(ypts,xpts);
dim_mesh = size(ypts_mesh);
OUTd_pts = zeros(length(xpts)*length(ypts),NT);
[OUTxseis_pts,OUTyseis_pts,OUTiglob_pts,OUTdseis_pts] = FindNearestNode(xpts_mesh(:),ypts_mesh(:),x,y);
xpts_mesh = reshape(OUTxseis_pts,dim_mesh);
ypts_mesh = reshape(OUTyseis_pts,dim_mesh);

%-- initialize data for output snapshots
OUTdt = round(NT/50); %round(NT/24); %50;
OUTit = 0;
OUTindx = Plot2dSnapshot(iglob);


%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark-alpha scheme with
% alpha=1/2, beta=1/2, gamma=1
%
half_dt = 0.5*dt;

%-- Absorbing boundaries (first order): 
[BcLC,iBcL] = BoundaryMatrix(wgll,[NELX NELY],iglob,dy_deta,'L'); % Left
[BcRC,iBcR] = BoundaryMatrix(wgll,[NELX NELY],iglob,dy_deta,'R'); % Right
[BcTC,iBcT] = BoundaryMatrix(wgll,[NELX NELY],iglob,dx_dxi ,'T'); % Top
[BcBC,iBcB] = BoundaryMatrix(wgll,[NELX NELY],iglob,dx_dxi ,'B'); % Bottom
% impedance = 1; % = sqrt(rho*mu)
BcLC = sqrt(rho_glob(iBcL).*mu_glob(iBcL)).*BcLC; %impedance*BcLC;
BcRC = sqrt(rho_glob(iBcR).*mu_glob(iBcR)).*BcRC; %impedance*BcRC;
BcTC = sqrt(rho_glob(iBcT).*mu_glob(iBcT)).*BcTC; %impedance*BcTC;
BcBC = sqrt(rho_glob(iBcB).*mu_glob(iBcB)).*BcBC; %impedance*BcBC;

% The mass matrix needs to be modified at the boundary
% for the implicit treatment of the term C*v.
% Fortunately C is diagonal.
M(iBcL) = M(iBcL) +half_dt*BcLC;
M(iBcR) = M(iBcR) +half_dt*BcRC;
M(iBcT) = M(iBcT) +half_dt*BcTC;
M(iBcB) = M(iBcB) +half_dt*BcBC;

for it=1:NT,

 % prediction of mid-step displacement:
 % d_mid = d_old + 0.5*dt*v_old
  d = d + half_dt*v; 

 % internal forces at mid-step -K*d(t+1/2) :
  a(:) = 0; % store -K*d in a global array
  for e=1:NEL,
   %switch to local (element) representation
    ig = iglob(:,:,e);
    local = d(ig);
   %gradients wrt local variables (xi,eta)
    d_xi  = Ht*local;	
    d_eta = local*H;
   %element contribution to internal forces
   %local = coefint1* H * ( W(:,:,e).*d_xi ) + coefint2* ( W(:,:,e).*d_eta ) *Ht ;
    wloc = W(:,:,e);
    d_xi = wloc.*d_xi;
    d_xi = H * d_xi;
    d_eta = wloc.*d_eta;
    d_eta = d_eta *Ht;
    local = coefint1* d_xi  + coefint2* d_eta ;
   %assemble into global vector
    a(ig) = a(ig) -local;
  end 

 % add external forces
  a(Fig) = a(Fig) + Ft(it);
  
 % absorbing boundaries:
  a(iBcL) = a(iBcL) - BcLC .* v(iBcL);
  a(iBcR) = a(iBcR) - BcRC .* v(iBcR);
  a(iBcT) = a(iBcT) - BcTC .* v(iBcT);
  a(iBcB) = a(iBcB) - BcBC .* v(iBcB);

 % acceleration: a = (-K*d +F)/M
  a = a ./M ;

 % update
 % v_new = v_old + dt*a_new;
 % d_new = d_old + dt*v_old + 0.5*dt^2*a_new
 %       = d_mid + 0.5*dt*v_new
  v = v + dt*a;
  d = d + half_dt*v;


%------------------------------------------
% STEP 4: OUTPUT
  if strcmp(avd_out,'a')
    avd = a;
  elseif strcmp(avd_out,'v')
    avd = v;
  elseif strcmp(avd_out,'d')
    avd = d;
  end
  OUTavd(:,it) = avd(OUTiglob);
  OUT_trF(:,it) = avd(Fig);  % Save trace at source position
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;

    figure(1) % seismograms
    PlotSeisTrace(OUTyseis,time,OUTavd);
%     PlotSeisTrace(OUTyseis,time,OUTd);

    figure(2)
%     Plot2dSnapshot(x,y,v,OUTindx,[-0.5 0.5]);
    Plot2dSnapshot(x,y,avd,OUTindx);
    hold on
    plot(OUTxseis,OUTyseis,'^k',Fx,Fy,'*k','linewidth',2)
    set(gca,'linewidth',1.5,'fontsize',15);
    colormap(bluewhitered)
    hold off
%     xlim([75 225])
   
%    if is_movieout
    rect = get(gcf,'Position'); rect(1:2) = [0 0]; OUTmovie(:,OUTit)=getframe(gcf,rect);
%    end

    drawnow
    

  end
  
  % Save pts in order to calculate traveltime surface
    OUTd_pts(:,it) = avd(OUTiglob_pts);
%     OUTd_pts(:,it) = d;
    
end % ... of time loop
OUTvs_pts = vs_glob(OUTiglob_pts);

%disp('To replay the movie: movie(OUTmovie)')
% if is_movieout
%     vid = VideoWriter('wavefield.avi');
%     vid.FrameRate = 5;
%     open(vid);
%     writeVideo(vid,OUTmovie);
%     close(vid);
% end

figure(3);
Plot2dSnapshot(x,y,vs_glob,OUTindx,[min(vs_glob)*0.99 max(vs_glob)*1.01]);
hold on
plot(OUTxseis,OUTyseis,'^k','MarkerFaceColor','w');
hold off
colormap(flipud(jet));
title('Model')
% xlim([75 225])
% [Y,X] = meshgrid(ynode,xnode);
% contourf(X,Y,mu_mat,'edgecolor','none');
% axis equal 
% axis tight
% title('Model')
% xlabel('X')
% ylabel('Y')

%% Estimate traveltime surface
% Reference trace at source location
tr_ref = OUT_trF;
tt = zeros(size(OUTd_pts,1),1); % traveltime vector
r_all = zeros(size(OUTd_pts,1),1); % straight ray distance from source
% grv_min = 3; % km/s
% grv_max = 5;
for ipt = 1:size(OUTd_pts,1)
    tr = OUTd_pts(ipt,:);
    r = sqrt((Fx-xpts_mesh(ipt)).^2+(Fy-ypts_mesh(ipt)).^2);
    r_all(ipt,1) = r;
%     if r>inf % group velocity window
%         Iwin = time<=r/grv_min & time>=r/grv_max;
%         tr(~Iwin) = 0;
%     end
    
    % Determine phase traveltime by cross-correlation
    [c,lags] = xcorr(tr,tr_ref);
%     [c,lags] = xcorr(abs(hilbert(tr)),abs(hilbert(tr_ref)));
    [~,Imax] = max(c);
    tt(ipt,1) = lags(Imax)*dt;
    
%     [~,I1] = max(tr_ref);
%     [~,I2] = max(tr);
%     tt(ipt,1) = time(I2)-time(I1);
    
    if 0
        figure(99); clf;
        plot(lags*dt,c); hold on;
        plot([tt tt],[0 c(Imax)],'-r');
    end
end
% colorbar('vert')

figure(4); clf;
set(gcf,'color','w');
% Plot2dSnapshot(xpts_mesh(:),ypts_mesh(:),tt,OUTindx,[min(tt) max(tt)]); hold on;
tt_mesh = reshape(tt,dim_mesh);
contourf(xpts_mesh,ypts_mesh,tt_mesh,linspace(min(tt),max(tt),100),'edgecolor','none'); hold on;
contour(xpts_mesh,ypts_mesh,tt_mesh,[0:5:Tduration],'-k','linewidth',2);
plot(OUTxseis,OUTyseis,'^r',Fx,Fy,'*k')
% colormap(flipud(jet));
title('Phase Delays')
axis equal 
axis tight
xlabel('X')
ylabel('Y')
cb = colorbar('vert');
ylabel(cb,'travel time (s)');
% xlim([75 225])
set(gca,'fontsize',12)

if 0
    figure;
%     imagesc(reshape(r_all,dim_mesh)./tt_mesh)
    imagesc((reshape(r_all,dim_mesh)./tt_mesh/4-1)*100)
    colorbar
%     caxis([0.99 1.01]*4)
    caxis(([0.99 1.01]-1)*100)
    colormap(redblue)
end

vs_mesh = reshape(OUTvs_pts,dim_mesh);

% OUTPUTS
out.tt_mesh = tt_mesh;
out.Fx = Fx;
out.Fy = Fy;
out.xpts_mesh = xpts_mesh;
out.ypts_mesh = ypts_mesh;
out.vs_mesh = vs_mesh;
out.xsta = OUTxseis;
out.ysta = OUTyseis;


end

