% Run SEM2D over phase velocity maps for ambient noise
%
% Loop through all stations and simulate a wavefield generated from an impulse
% at each station
%
% % % % % % % % % % % % % % % % % % % % % %
% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% with stress free boundary conditions,
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
clear
setup_parameters_SEM2D;
is_savemat = 1; % Save matfile containing phase delays
is_movieout = 1;
is_overwrite = 0; % overwrite previous runs?
%------------------------------------------
% per_ind = [3 9 12 14 17 19];
% per_ind = [3, 14];
per_ind = [3]; % 12.7764 s

% Get desired frequencies (do not need to be the same as input synthetic maps)
load(['./A01_A03_xsp.mat']);
Tperiods = xspinfo.per_start;

% Load phase velocity maps
load('Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_phvmaps_Line13.mat','map');

% Output path
out_sem2d = './out_SEM2D/';
figpath = './figs/';

station_list = './sta_locs.txt';
% Load station info
[sta.name, sta.lat, sta.lon, sta.dep] = textread(station_list,'%s %f %f %f');


% Convert lat-lon grid to x-y taking into account ellipsoid
lon_gr = map.lon;
lat_gr = map.lat;
olon = min(lon_gr(:));
olat = min(lat_gr(:));
[ x_gr, y_gr, z_gr] = geodetic2enu(lat_gr,lon_gr,0,olat,olon,0,wgs84Ellipsoid);
x_gr = x_gr/1000;
y_gr = y_gr/1000;
z_gr = z_gr/1000;

% Interpolate irregular x-y grid to regular grid
x = linspace(min(x_gr(:)),max(x_gr(:)),size(x_gr,1));
y = linspace(min(y_gr(:)),max(y_gr(:)),size(x_gr,2));
[Y,X] = meshgrid(y,x);
F = scatteredInterpolant(x_gr(:),y_gr(:),z_gr(:));
Z   = F(X,Y);


% % Convert Lat-Lon to X-Y using true ellipsoid
% X = zeros(size(map.lon));
% Y = zeros(size(map.lat));
% for ix=1:size(map.lon,1)
%     for iy=1:size(map.lon,2)
%         [ X(ix,iy), Y(ix,iy) , ~] = lonlat2xy_nomap( olon, olat, map.lon(ix,iy), map.lat(ix,iy) );
%     end
% end

% % Lazy convert Lat-Lon to X-Y assuming great-circle paths on spherical earth...
% X = deg2km(map.lon-olon);
% Y = deg2km(map.lat-olat);
[NELX, NELY] = size(X);
LX = max(X(:));
LY = max(Y(:));

% Set station parameters
[ OUTxseis, OUTyseis, OUTzseis] = geodetic2enu(sta.lat,sta.lon,0,olat,olon,0,wgs84Ellipsoid);
OUTxseis = OUTxseis/1000;
OUTyseis = OUTyseis/1000;
OUTzseis = OUTzseis/1000;
% OUTxseis = deg2km(sta.lon-olon);
% OUTyseis = deg2km(sta.lat-olat);

%**** Set here properties of medium : ****
% Interpolate to desired freqeuncy axis
[ map ] = Tinterp_map( map, Tperiods );

if ~exist(figpath)
    mkdir(figpath);
end
%% Loop over stations and calculate 2D SEM
% for iper = 1:length(Tperiods)
for iper = per_ind
    VS = map.phv(:,:,iper);
    
    % Interpolate VS to new grid
    F = scatteredInterpolant(x_gr(:),y_gr(:),VS(:));
    VS = F(X,Y);
    
    figure(99); clf; box on; hold on;
    set(gcf,'color','w');
    contourf(X,Y,VS,100,'LineStyle','none');
    plot(OUTxseis,OUTyseis,'^k','MarkerFaceColor','w','linewidth',1,'markersize',12);
    cb = colorbar;
    cb.LineWidth = 1.5;
    ylabel(cb,'Phase Velocity (km/s)');
    set(gca,'fontsize',16,'linewidth',1.5,'layer','top');
    colormap(tomo_cmap(100));
    xlabel('X (km)');
    ylabel('Y (km)');
    title([num2str(Tperiods(iper)),' s']);
    save2pdf([figpath,'/PhaseVelocityMap_',num2str(Tperiods(iper)),'s','.pdf'],99,250);
    drawnow
    
    for ista = 1:length(OUTxseis)
        clear run
        disp(['Working on ',num2str(Tperiods(iper)),'s: STA ',sta.name{ista},' ',num2str(ista),'/',num2str(length(OUTxseis))])
        % Set up paths
        out_path = [out_sem2d,'/',num2str(Tperiods(iper)),'/'];
        mat_path = [out_path,'phasedelay_',sta.name{ista},'_',Fstr,'_',num2str(1/Ff0),'s_P',num2str(P),'_',num2str(Tperiods(iper)),'s.mat'];
        if exist(mat_path) && is_overwrite==0
            disp(['Skip finished ',sta.name{ista}]);
            continue
        end
        if is_savemat && ~exist(out_path)
            mkdir(out_path)
        end
        
        %**** Set here the source location : ****
        Fx = OUTxseis(ista); Fy = OUTyseis(ista);

    %   Run SEM-2D
        tic
        [out, OUTmovie] = run_SEM2D( LX, LY, NELX, NELY, P, Fx, Fy, avd_out, Ff0, Fstr, VS, OUTxseis, OUTyseis);
        toc
        run.out = out;
        
        % Convert output back to Lat-Lon mesh
        tt_mesh = out.tt_mesh;
        vs_mesh = out.vs_mesh;
        xpts_mesh = out.xpts_mesh;
        ypts_mesh = out.ypts_mesh;
        zpts_mesh = interp2(Y,X,Z,ypts_mesh,xpts_mesh);
        out.zpts_mesh = zpts_mesh;
        [lat_gr,lon_gr,~] = enu2geodetic(xpts_mesh*1000,ypts_mesh*1000,zpts_mesh*1000,olat,olon,0,wgs84Ellipsoid);
        % Interpolate irregular lat-lon grid to regular grid
        lon = linspace(min(lon_gr(:)),max(lon_gr(:)),size(lon_gr,1));
        lat = linspace(min(lat_gr(:)),max(lat_gr(:)),size(lon_gr,2));
        [lat_mesh,lon_mesh] = meshgrid(lat,lon);
        F = scatteredInterpolant(lat_gr(:),lon_gr(:),vs_mesh(:));
        lalo.vs_mesh = F(lat_mesh,lon_mesh);
        F = scatteredInterpolant(lat_gr(:),lon_gr(:),tt_mesh(:));
        lalo.tt_mesh = F(lat_mesh,lon_mesh);
        lalo.lat_mesh = lat_mesh;
        lalo.lon_mesh = lon_mesh;
        [lalo.F_lat, lalo.F_lon,~] = enu2geodetic(Fx*1000,Fy*1000,OUTzseis(ista)*1000,olat,olon,0,wgs84Ellipsoid);
        
        run.lalo = lalo;
        
        if 1
            figure(5); clf;
            set(gcf,'color','w');
            
            subplot(1,2,1);
            contourf(out.xpts_mesh,out.ypts_mesh,out.tt_mesh,linspace(min(out.tt_mesh(:)),max(out.tt_mesh(:)),100),'edgecolor','none'); hold on;
            contour( out.xpts_mesh,out.ypts_mesh,out.tt_mesh,[0:5:max(out.tt_mesh(:))],'-k','linewidth',2);
            plot(OUTxseis,OUTyseis,'^r',Fx,Fy,'*k')
            title('Phase Delays')
            axis equal 
            axis tight
            xlabel('X')
            ylabel('Y')
            cb1 = colorbar('vert');
            ylabel(cb1,'travel time (s)');
            set(gca,'fontsize',12);
            
            subplot(1,2,2);
            contourf(lalo.lon_mesh,lalo.lat_mesh,lalo.tt_mesh,linspace(min(out.tt_mesh(:)),max(out.tt_mesh(:)),100),'edgecolor','none'); hold on;
            contour( lalo.lon_mesh,lalo.lat_mesh,lalo.tt_mesh,[0:5:max(out.tt_mesh(:))],'-k','linewidth',2);
            [OUTlatseis, OUTlonseis,~] = enu2geodetic(OUTxseis*1000,OUTyseis*1000,OUTzseis*1000,olat,olon,0,wgs84Ellipsoid);
            plot(OUTlonseis,OUTlatseis,'^r',lalo.F_lon,lalo.F_lat,'*k')
            title('Phase Delays')
            axis equal 
            axis tight
            xlabel('Lon')
            ylabel('Lat')
            cb2 = colorbar('vert');
            ylabel(cb2,'travel time (s)');
            set(gca,'fontsize',12);
            
            save2pdf([figpath,'/PhaseDelays_coordtrans_',sta.name{ista},'_',num2str(Tperiods(iper)),'s','.pdf'],5,100);
        end
        
        itstr = [sta.name{ista},'_',num2str(Tperiods(iper)),'s'];
        save2pdf([figpath,'/PhaseDelays_',itstr,'.pdf'],4,100);
        if is_movieout
            vid = VideoWriter([figpath,'/wavefield_',itstr,'.avi']);
            vid.FrameRate = 10;
            open(vid);
            writeVideo(vid,OUTmovie);
            close(vid);
        end
        
        % run(1).periods = Tperiods;
        run(1).period = Tperiods(iper);
        run(1).sta = sta;
        run(1).olon = olon;
        run(1).olat = olat;
        run(1).lon = map.lon;
        run(1).lat = map.lat;
        
        if is_savemat
            % Save to mat
            save(mat_path,'run')
        end
    end
    save2pdf([figpath,'/VS_',num2str(Tperiods(iper)),'s.pdf'],3,100);
end
