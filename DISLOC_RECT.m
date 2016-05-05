%% SURFACE DISPLACEMENT DUE TO RECTANGULAR FAULT
%-------------------
%	OLAT, OLON  : Coordinates of observation points (Latitude, Longitude)
%	ALAT0,ALON0 : Coordinates of the fault centroid (Latitude, Longitude)
%	DEP         : Depth of the fault top (DEP > 0)
%	PHAI        : Strike-angle from Norh clock wise (in degrees)
%                 fault dips to the right side of the trace
%	DIP         : Dip-angle  (in degrees, must be scalar)
%	RAK         : Rake-angle (in degrees, must be scalar)
%	AL          : Fault length in strike direction (LEN > 0)
%	AW          : Fault width in dip direction (WIDTH > 0)
%	NU          : 0.25 isotropic Poisson's ratio
%-------------------
% coded by T.Ito AUG 2013
% modified by T.Ito Feb. 2016
% modified by T.Ito June 2016
function [U]=DISLOC_RECT(OLAT,OLON,ALAT0,ALON0,DEP,PHAI,DIP,RAK,AL,AW,NU)
RAD=pi./180;
nsite=length(OLAT);
%
SR=sin(RAK.*RAD);
CR=cos(RAK.*RAD);
%
SD=sin(DIP.*RAD);
CD=cos(DIP.*RAD);
%
ST=sin((PHAI-90).*RAD);
CT=cos((PHAI-90).*RAD);
%
CD(abs(CD)<=1.0e-3)=0;
%
SD(CD==0 & SD > 0)=1;
SD(CD==0 & SD < 0)=-1;
%
[Xobs,Yobs]=PLTXY(OLAT,OLON,ALAT0,ALON0);
%
%	X,Y : Coordinates of observation points (relative to fault center at top)
%
Xobs=Xobs.*1000;  % Unit km --> m
Yobs=Yobs.*1000;  % Unit km --> m
DEP=DEP*1000;  % Unit km --> m
AL=AL*1000;    % Unit km --> m
AW=AW*1000;    % Unit km --> m
%
% Converts fault coordinates (E,N,DEPTH) relative to centroid
% into Okada's reference system (X,Y,DEP)
% DEP is fault's top edge
EC = Xobs + CT.*CD.*AW/2;
NC = Yobs - ST.*CD.*AW/2;
X  = CT.*NC - ST.*EC + AL/2;
Y  = ST.*NC + CT.*EC + CD.*AW;
% Variable substitution (independent from xi and eta)
P = Y.*CD + DEP.*SD.*ones(nsite,1);
Q = Y.*SD - DEP.*CD.*ones(nsite,1);
%
% XI,ET,Q : FAULT COORDINATE
%           ORIGIN AT CENTER OF LOWER SIDE OF FAULT
%
Uxs=zeros(nsite,1); Uys=zeros(nsite,1); Uzs=zeros(nsite,1);
Uxd=zeros(nsite,1); Uyd=zeros(nsite,1); Uzd=zeros(nsite,1);
Uxt=zeros(nsite,1); Uyt=zeros(nsite,1); Uzt=zeros(nsite,1);
%
[uxs,uys,uzs,uxd,uyd,uzd,uxt,uyt,uzt]=STATIC(NU,X,P,Q,SD,CD);
Uxs=Uxs+uxs; Uys=Uys+uys; Uzs=Uzs+uzs;
Uxd=Uxd+uxd; Uyd=Uyd+uyd; Uzd=Uzd+uzd;
Uxt=Uxt+uxt; Uyt=Uyt+uyt; Uzt=Uzt+uzt;
%
[uxs,uys,uzs,uxd,uyd,uzd,uxt,uyt,uzt]=STATIC(NU,X,P-AW,Q,SD,CD);
Uxs=Uxs-uxs; Uys=Uys-uys; Uzs=Uzs-uzs;
Uxd=Uxd-uxd; Uyd=Uyd-uyd; Uzd=Uzd-uzd;
Uxt=Uxt-uxt; Uyt=Uyt-uyt; Uzt=Uzt-uzt;
%
[uxs,uys,uzs,uxd,uyd,uzd,uxt,uyt,uzt]=STATIC(NU,X-AL,P,Q,SD,CD);
Uxs=Uxs-uxs; Uys=Uys-uys; Uzs=Uzs-uzs;
Uxd=Uxd-uxd; Uyd=Uyd-uyd; Uzd=Uzd-uzd;
Uxt=Uxt-uxt; Uyt=Uyt-uyt; Uzt=Uzt-uzt;
%
[uxs,uys,uzs,uxd,uyd,uzd,uxt,uyt,uzt]=STATIC(NU,X-AL,P-AW,Q,SD,CD);
Uxs=Uxs+uxs; Uys=Uys+uys; Uzs=Uzs+uzs;
Uxd=Uxd+uxd; Uyd=Uyd+uyd; Uzd=Uzd+uzd;
Uxt=Uxt+uxt; Uyt=Uyt+uyt; Uzt=Uzt+uzt;
%
U.E=(-CR.*( Uxs(:).*CT+Uys(:).*ST)-SR.*( Uxd(:).*CT+Uyd(:).*ST))./(2*pi); %E
U.N=(-CR.*(-Uxs(:).*ST+Uys(:).*CT)-SR.*(-Uxd(:).*ST+Uyd(:).*CT))./(2*pi); %N
U.U=(-CR.*( Uzs(:)               )-SR.*( Uzd(:)               ))./(2*pi); %D
U.Et=( Uxt(:).*CT+Uyt(:).*ST)./(2*pi); %E
U.Nt=(-Uxt(:).*ST+Uyt(:).*CT)./(2*pi); %N
U.Ut=(Uzt(:)               )./(2*pi); %D
%
end
%====================================================
function [Uxs,Uys,Uzs,Uxd,Uyd,Uzd,Uxt,Uyt,Uzt]=STATIC(nu,xi,eta,q,sd,cd)
%-------------------
% Subroutine STATIC
% coded by T.Ito AUG 2013
% modified by T.Ito Feb. 2016
%-------------------
q2  = q.^2;
X   = sqrt(xi.^2 + q2);
R   = sqrt(xi.^2 + eta.^2 + q2);
yb  = eta.*cd + q.*sd;
db  = eta.*sd - q.*cd;
Rb  = R + db;
cR  = cd.*Rb;
Re  = R + eta;
lR  = log(Re);
qre = q./(R.*Re);
qrx = q./(R.*(R + xi));
NU12=1-2*nu;
if cd == 0
  I5 =  NU12 * 2./cd .* atan((eta.*(X + q.*cd) + X.*(R + X).*sd)) ./(xi.*(R + X).*cd);
  I5(xi==0) = 0;
  I4 =  NU12 * 1./cd * (log(Rb) - sd.*lR);
  I3 =  NU12 * ( yb./cR) + sd./cd.*I4 - lR ;
  I1 =  NU12 * (-xi./cR) - sd./cd.*I5;
else
  I5 = -NU12   * xi.*sd./Rb;
  I4 = -NU12   * q./Rb;
  I3 =  NU12/2 * (eta./Rb + yb.*q./Rb.^2 - lR);
  I1 = -NU12/2 * xi.*q./Rb.^2;
end
I2   =  NU12   * (-lR) - I3;
%
Uxs  =  xi.*qre + I1.*sd;
Uys  =  yb.*qre + q.*cd./Re  + I2.*sd;
Uzs  =  db.*qre + q.*sd./Re  + I4.*sd;
Uxd  =   q./R                - I3.*sd.*cd;
Uyd  =  yb.*qrx - I1.*sd.*cd;
Uzd  =  db.*qrx - I5.*sd.*cd;
Uxt  =   q.*qre               - I3.*sd.*sd;
Uyt  = -db.*qrx - sd.*xi.*qre - I1.*sd.*sd;
Uzt  =  yb.*qrx + cd.*xi.*qre - I5.*sd.*sd;
%
k    = find(q~=0);
ax   = atan(xi(k).*eta(k)./(q(k).*R(k)));
Uxs(k) = Uxs(k) +     ax;
Uyd(k) = Uyd(k) + cd.*ax;
Uzd(k) = Uzd(k) + sd.*ax;
Uyt(k) = Uyt(k) + sd.*ax;
Uzt(k) = Uzt(k) - cd.*ax;
end