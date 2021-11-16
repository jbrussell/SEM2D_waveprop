function [ x, y , z] = lonlat2xy_nomap( olon, olat, lon, lat,h )
% [ x, y, z ] = lonlat2xy_nomap( olon, olat, lon, lat,h )
%
%  Identical usage to lonlat2xy, but for folks without the mapping toolbox!
% 
% Convert lon/lat (deg) to x/y (m) using the reference ellipsoid

if nargin<5 || isempty(h)
    h=0;
end

global ellipsoidGRS80
if isempty(ellipsoidGRS80)
    ellipsoidGRS80 = struct( 'Name','World Geodetic System 1984',...
                               'LengthUnit','meter',...
                               'SemimajorAxis',6378137,...
                               'SemiminorAxis',6356752.314245179,...
                               'InverseFlattening',298.257223563,...
                               'Eccentricity',0.081819190842621,...
                               'Flattening',0.003352810664747,...
                               'ThirdFlattening',0.001679220386384,...
                               'MeanRadius',6.371008771415059e+06,...
                               'SurfaceArea',5.100656217240886e+14,...
                               'Volume',1.083207319801408e+21 );
end

%% If mapping toolbox...
% [x, y, ~] = geodetic2enu(lat, lon, zeros(size(lon)), olat, olon, 0, ellipsoidGRS80);

%% --------- If no mapping toolbox --------
% https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#Geodetic_to/from_ENU_coordinates

% geodetic to ECEF
[Xo,Yo,Zo] = geodetic2ecef_ze(olat,olon,h,ellipsoidGRS80);
[X,Y,Z] = geodetic2ecef_ze(lat,lon,h,ellipsoidGRS80);
[x,y,z] = ecef2nu_ze(olat,olon,Xo,Yo,Zo,X,Y,Z);

% [ x, y ] = dontuse_lonlat2xy( olon, olat, lon, lat ); % test old version
end

function [X,Y,Z] = geodetic2ecef_ze(lat,lon,h,LIPSD)
a = LIPSD.SemimajorAxis;
b = LIPSD.SemiminorAxis;
e = LIPSD.Eccentricity;
N_phi = a./sqrt(1 - e^2*sind(lat).^2);
X = (N_phi + h).*cosd(lat).*cosd(lon);
Y = (N_phi + h).*cosd(lat).*sind(lon);
Z = (b*b*N_phi/a/a + h).*sind(lat);
end

function [x,y,z] = ecef2nu_ze(olat,olon,Xo,Yo,Zo,X,Y,Z)
G = [ -sind(olon)             cosd(olon)            0
      -sind(olat)*cosd(olon) -sind(olat)*sind(olon) cosd(olat)
       cosd(olat)*cosd(olon)  cosd(olat)*sind(olon) sind(olat)];

dX = X(:)-Xo;
dY = Y(:)-Yo;
dZ = Z(:)-Zo;

r = G*[dX'; dY'; dZ'];
x = r(1,:)';
y = r(2,:)';
z = r(3,:)';
end
