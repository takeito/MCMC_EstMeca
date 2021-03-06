%% PLTXY TRANSFORMS (ALAT,ALONG) TO (X,Y)
% INPUT
% ALAT: LATTITUDE (degree)
% ALON: LONGITUDE (degree)
% ALAT0: ORIGIN OF LATTITUDE (degree)
% ALON0: ORIGIN OF LONGITUDE (degree)
% OUTPUT
% X: EAST DIRECTION AT LOCAL CARTESIAN  (km)
% Y: NORTH DIRECTION AT LOCAL CARTESIAN (km)
%
function [X,Y]=PLTXY(ALAT,ALON,ALAT0,ALON0)
A  =6.378160e3;
E2 =6.6944541e-3;
E12=6.7395719e-3;
D=180/pi;
RD=1.0/D;
RLAT = RD.*ALAT;
SLAT = sin(RLAT);
CLAT = cos(RLAT);
AL   = ALON-ALON0;
PH1  = ALAT + ((1.0 + E12.*CLAT.^2).*AL.^2.*SLAT.*CLAT)./(2.0*D);
AN   = A./sqrt(1.0-E2.*sin(PH1.*RD).^2);
Y    = (A.*(1.0-E2)./sqrt((1.0-E2.*sin((PH1+ALAT0).*0.5.*RD).^2).^3)).*(PH1-ALAT0).*RD;
X    = (AL.*CLAT)./(D./AN)+(AL.^3.*CLAT.*cos(2.0.*RLAT))./(6.0./AN.*D.^3);
end