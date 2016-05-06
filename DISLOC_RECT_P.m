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
function [U]=DISLOC_RECT_P(OLAT,OLON,ALAT0,ALON0,DEP,PHAI,DIP,RAK,AL,AW,NU)
RAD=pi./180;
nsite=length(OLAT);
nmodel=length(ALAT0);
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
[Xobs,Yobs]=PLTXY_P(OLAT,OLON,ALAT0,ALON0);
%
%	X,Y : Coordinates of observation points (relative to fault center at top)
%
Xobs=Xobs.*1000;  % Unit km --> m
Yobs=Yobs.*1000;  % Unit km --> m
DEP=DEP.*1000;  % Unit km --> m
AL=AL.*1000;    % Unit km --> m
AW=AW.*1000;    % Unit km --> m
AL2=AL/2;
AW2=AW/2;
%
% Converts fault coordinates (E,N,DEPTH) relative to centroid
% into Okada's reference system (X,Y,DEP)
% DEP is fault's top edge
%EC = bsxfun(@plus, Xobs,CT.*CD.*AW2);
%NC = bsxfun(@minus,Yobs,ST.*CD.*AW2);
%X =  bsxfun(@times,EC,CT) - bsxfun(@plus,bsxfun(@times,NC,ST),AL./2);
%Y =  bsxfun(@times,EC,ST) + bsxfun(@plus,bsxfun(@times,NC,CT),CD.*AW);
X =  bsxfun(@times,Xobs,CT) - bsxfun(@plus,bsxfun(@times,Yobs,ST),AL./2);
Y =  bsxfun(@times,Xobs,ST) + bsxfun(@plus,bsxfun(@times,Yobs,CT),CD.*AW);
% Variable substitution (independent from xi and eta)
P = bsxfun(@plus, bsxfun(@times,Y,CD),(DEP+SD.*AW).*SD);
Q = bsxfun(@minus,bsxfun(@times,Y,SD),(DEP+SD.*AW).*CD);
%
% XI,ET,Q : FAULT COORDINATE
%           ORIGIN AT CENTER OF LOWER SIDE OF FAULT
%
Uxs=zeros(nsite,nmodel); Uys=zeros(nsite,nmodel); Uzs=zeros(nsite,nmodel);
Uxd=zeros(nsite,nmodel); Uyd=zeros(nsite,nmodel); Uzd=zeros(nsite,nmodel);
Uxt=zeros(nsite,nmodel); Uyt=zeros(nsite,nmodel); Uzt=zeros(nsite,nmodel);
%
[uxs,uys,uzs,uxd,uyd,uzd,uxt,uyt,uzt]=STATIC_P(NU,bsxfun(@plus,X,AL2),P,Q,SD,CD);
Uxs=Uxs+uxs; Uys=Uys+uys; Uzs=Uzs+uzs;
Uxd=Uxd+uxd; Uyd=Uyd+uyd; Uzd=Uzd+uzd;
Uxt=Uxt+uxt; Uyt=Uyt+uyt; Uzt=Uzt+uzt;
%
[uxs,uys,uzs,uxd,uyd,uzd,uxt,uyt,uzt]=STATIC_P(NU,bsxfun(@plus,X,AL2),bsxfun(@minus,P,AW),Q,SD,CD);
Uxs=Uxs-uxs; Uys=Uys-uys; Uzs=Uzs-uzs;
Uxd=Uxd-uxd; Uyd=Uyd-uyd; Uzd=Uzd-uzd;
Uxt=Uxt-uxt; Uyt=Uyt-uyt; Uzt=Uzt-uzt;
%
[uxs,uys,uzs,uxd,uyd,uzd,uxt,uyt,uzt]=STATIC_P(NU,bsxfun(@minus,X,AL2),P,Q,SD,CD);
Uxs=Uxs-uxs; Uys=Uys-uys; Uzs=Uzs-uzs;
Uxd=Uxd-uxd; Uyd=Uyd-uyd; Uzd=Uzd-uzd;
Uxt=Uxt-uxt; Uyt=Uyt-uyt; Uzt=Uzt-uzt;
%
[uxs,uys,uzs,uxd,uyd,uzd,uxt,uyt,uzt]=STATIC_P(NU,bsxfun(@minus,X,AL2),bsxfun(@minus,P,AW),Q,SD,CD);
Uxs=Uxs+uxs; Uys=Uys+uys; Uzs=Uzs+uzs;
Uxd=Uxd+uxd; Uyd=Uyd+uyd; Uzd=Uzd+uzd;
Uxt=Uxt+uxt; Uyt=Uyt+uyt; Uzt=Uzt+uzt;
%
U.E = (-bsxfun(@times,(bsxfun(@times, Uxs,CT)+bsxfun(@times,Uys,ST)),CR)...
       -bsxfun(@times,(bsxfun(@times, Uxd,CT)+bsxfun(@times,Uyd,ST)),SR))./(2*pi);
U.N = (-bsxfun(@times,(bsxfun(@times,-Uxs,ST)+bsxfun(@times,Uys,CT)),CR)...
       -bsxfun(@times,(bsxfun(@times,-Uxd,ST)+bsxfun(@times,Uyd,CT)),SR))./(2*pi);
U.U = (-bsxfun(@times,Uzs,CR)-bsxfun(@times,Uzd,SR))./(2*pi);
U.Et= ( bsxfun(@times,Uxt,CT)+bsxfun(@times,Uyt,ST))./(2*pi);
U.Nt= (-bsxfun(@times,Uxt,ST)+bsxfun(@times,Uyt,CT))./(2*pi);
U.Ut= Uzt./(2*pi);
%
end
%====================================================
function [Uxs,Uys,Uzs,Uxd,Uyd,Uzd,Uxt,Uyt,Uzt]=STATIC_P(nu,xi,eta,q,sd,cd)
%-------------------
% Subroutine STATIC_P
% coded by T.Ito AUG 2013
% modified by T.Ito June 2016
%-------------------
X = sqrt(xi.^2 + q.^2);
R  =sqrt(xi.^2 + eta.^2 + q.^2);
yb = bsxfun(@times,eta,cd) + bsxfun(@times,q,sd);
db = bsxfun(@times,eta,sd) - bsxfun(@times,q,cd);
Rb = R + db;
cR = bsxfun(@times,Rb,cd);
Re = R + eta;
lR = log(Re);
qre = q./(R.*Re);
qrx = q./(R.*(R + xi));
sdsd=sd.*sd;
sdcd=sd.*cd;
xiqre=xi.*qre;
NU12=1-2*nu;
%
I5 = -NU12   * bsxfun(@times,xi,sd)./Rb;
I4 = -NU12   * q./Rb;
I3 =  NU12/2 * (eta./Rb + yb.*q./Rb.^2 - lR);
I2 =  NU12   * (-lR) - I3;
I1 = -NU12/2 * xi.*q./Rb.^2;
%
IND=cd==0;
I5(:,IND) = bsxfun(@times,...
            atan((eta(:,IND).*(X(:,IND)+bsxfun(@times,q(:,IND),cd(IND))+...
                    X(:,IND).* bsxfun(@times,R(:,IND)+X(:,IND),sd(IND)))))./...
                  (xi(:,IND).*(bsxfun(@times,R(:,IND)+X(:,IND),cd(IND)))),NU12*2./cd(IND));
I5(xi==0) = 0;
I4(:,IND) = bsxfun(@times,log(Rb(:,IND)) - bsxfun(@times,lR(:,IND),sd(IND)),NU12./cd(IND));
I3(:,IND) = NU12*( yb(:,IND)./cR(:,IND)) + bsxfun(@times,I4(:,IND),sd(IND)./cd(IND)) - lR(:,IND) ;
I1(:,IND) = NU12*(-xi(:,IND)./cR(:,IND)) - bsxfun(@times,I5(:,IND),sd(IND)./cd(IND));
%
Uxs  =    xiqre + bsxfun(@times,I1,sd);
Uys  =  yb.*qre + q.*bsxfun(@ldivide,Re,cd) + bsxfun(@times,I2,sd);
Uzs  =  db.*qre + q.*bsxfun(@ldivide,Re,cd) + bsxfun(@times,I4,sd);
Uxd  =     q./R - bsxfun(@times,I3,sdcd);
Uyd  =  yb.*qrx - bsxfun(@times,I1,sdcd);
Uzd  =  db.*qrx - bsxfun(@times,I5,sdcd);
Uxt  =   q.*qre - bsxfun(@times,I3,sdsd);
Uyt  = -db.*qrx - bsxfun(@times,xiqre,sd) - bsxfun(@times,I1,sdsd);
Uzt  =  yb.*qrx + bsxfun(@times,xiqre,cd) - bsxfun(@times,I5,sdsd);
%
IND=find(q==0);
ax  =atan(xi.*eta./(q.*R));
axcd=bsxfun(@times,ax,cd);
axsd=bsxfun(@times,ax,sd);
ax(IND)=0;
axcd(IND)=0;
axsd(IND)=0;
%
Uxs = Uxs + ax;
Uyd = Uyd + axcd;
Uzd = Uzd + axsd;
Uyt = Uyt + axsd;
Uzt = Uzt - axcd;
end
%% PLTXY TRANSFORMS (ALAT,ALONG) TO (X,Y) IN PALLALEL
% INPUT
% ALAT: LATTITUDE (degree)
% ALON: LONGITUDE (degree)
% ALAT0: ORIGIN OF LATTITUDE (degree)
% ALON0: ORIGIN OF LONGITUDE (degree)
% OUTPUT
% X: EAST DIRECTION AT LOCAL CARTESIAN  (km)
% Y: NORTH DIRECTION AT LOCAL CARTESIAN (km)
%-------------------
% coded by T.Ito June 2016
function [X,Y]=PLTXY_P(ALAT,ALON,ALAT0,ALON0)
A  =6.378160e3;
E2 =6.6944541e-3;
E12=6.7395719e-3;
D=180/pi;
RD=1.0/D;
RLAT = RD.*ALAT;
SLAT = sin(RLAT);
CLAT = cos(RLAT);
AL   = bsxfun(@minus,ALON,ALON0);
PH1  = bsxfun(@plus,bsxfun(@times,AL.^2,((1.0 + E12.*CLAT.^2).*SLAT.*CLAT)./(2.0*D)),ALAT);
AN   = A./sqrt(1.0-E2.*sin(PH1.*RD).^2);
PAR  = bsxfun(@plus,PH1,ALAT0).*RD;
Y    = A.*(1.0-E2)./sqrt((1.0-E2.*sin(PAR.*0.5).^2).^3).*PAR;
X    = bsxfun(@ldivide,bsxfun(@minus,AL,CLAT),D./AN)+bsxfun(@ldivide,bsxfun(@times,AL.^3,CLAT.*cos(2.0.*RLAT)),(6.0./AN.*D.^3));
end