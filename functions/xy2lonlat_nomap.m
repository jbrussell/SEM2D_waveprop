function [ lon, lat ] = xy2lonlat_nomap( olon, olat, x, y,az )
%[ lon, lat ] = xy2lonlat( olon, olat, x, y,az )
%
%  Identical usage to xy2lonlat, but for folks without the mapping toolbox!
% 
% Convert x/y (m) to lat/lon (degrees) using the reference ellipsoid
% option to rotate by azimuth "az" (defined as clockwise from N)

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

if nargin < 5
    az = 0;
end

if length(x)==1 && length(y)>1
    x = x.*ones(size(y));
end
if length(y)==1 && length(x)>1
    y = y.*ones(size(x));
end

if az~=0
    x0 = x;
    y0 = y;
    x = x0.*cosd(az) - y0.*sind(az);
    y = x0.*sind(az) + y0.*cosd(az);
end

%% If mapping toolbox...
% [lat, lon, ~] = enu2geodetic(x, y, zeros(size(x)), olat, olon, 0, ellipsoidGRS80);

%% --------- If no mapping toolbox --------
% https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#Geodetic_to/from_ENU_coordinates

% geodetic to ECEF
[Xo,Yo,Zo] = geodetic2ecef_ze(olat,olon,0,ellipsoidGRS80);
[X,Y,Z] = enu2ecef_ze(olat,olon,x,y,0,Xo,Yo,Zo);
[lon,lat,~] = ecef2geodetic_ze(X,Y,Z,ellipsoidGRS80);

% [ lon, lat ] = dontuse_xy2lonlat( olon, olat, x, y ); % test old version
end

function [lon,lat,h] = ecef2geodetic_ze(X,Y,Z,LIPSD)
a = LIPSD.SemimajorAxis;
b = LIPSD.SemiminorAxis;
e = LIPSD.Eccentricity;

% Kaplan 15 step procedure
r = sqrt(X.^2 + Y.^2);
ee2 = (a.^2 - b.^2)./b.^2;
e2 = e.^2;
e4 = e.^4;
E2 = a.^2 - b.^2;
F = 54*b*b*Z.*Z;
G = r.^2 + (1 - e2)*Z.^2 - e2*E2;
C = e4.*F.*r.^2./G.^3;
S = (1 + C + sqrt(C.^2 + 2*C)).^(1/3);
P = F./(3*(S + 1./S + 1).^2.*G.*G);
Q = sqrt(1 + 2*e4*P);
r01 = -(P.*e2.*r)./(1+Q);
r02b = (P.*(1 - e2).*Z.*Z)./(Q.*(1+Q));
r02 = sqrt(0.5*a*a*(1 + 1./Q) - r02b - 0.5*P.*r.*r);
r0 = r01 + r02;
U = sqrt((r - e2.*r0).^2 + Z.^2);
V = sqrt((r - e2.*r0).^2 + (1 - e2).*Z.^2);
Z0 = (b.*b.*Z)./(a.*V);
h = U.*(1 - (b.^2)./(a.*V));
lat = atand((Z + ee2*Z0)./r);
lon = atan2d(Y,X);
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

function [X,Y,Z] = enu2ecef_ze(olat,olon,x,y,z,Xo,Yo,Zo)
if length(z)==1, z = ones(size(x))*z; end

G = [ -sind(olon)             cosd(olon)            0
      -sind(olat)*cosd(olon) -sind(olat)*sind(olon) cosd(olat)
       cosd(olat)*cosd(olon)  cosd(olat)*sind(olon) sind(olat)];

R = G'*[x(:)'; y(:)'; z(:)'];
X = R(1,:)' + Xo;
Y = R(2,:)' + Yo;
Z = R(3,:)' + Zo;
end


