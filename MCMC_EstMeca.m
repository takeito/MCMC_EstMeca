function MCMC_EstMeca(infile)
%% READ PARAMETERS
PRM=READ_PARM_EstMeca(infile);
%% Estimate Fault Parameters
for NUM=1:10
  est_MCMC_CPU(PRM,NUM);
end
READ_KEEP(infile)
end
%% ESTIMATE PARAMETER BASED ON MCMC 
function est_MCMC_CPU(PRM,NUM)
% coded by T.Ito Nov.20 2013
[LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,RES,DRP,MW,dPr]=INT_PARM_MCMC(PRM);
LDIM=PRM.NPL.*PRM.KEEP;
PRM.NERN=1./PRM.ERRN;
PRM.NERE=1./PRM.ERRE;
PRM.NERU=1./PRM.ERRU;
Rtime=0;
DPP=1;
RWD=1;
tic;
%%
d2r=pi./180;
SLIP_APR1=repmat(SLIP_VEC(d2r*PRM.ASTR(1),d2r*PRM.ADIP(1),d2r*PRM.ARAK(1)),1,PRM.NPL);
SLIP_APR2=repmat(SLIP_VEC(d2r*PRM.ASTR(2),d2r*PRM.ADIP(2),d2r*PRM.ARAK(2)),1,PRM.NPL);
%%
RT=0;
COUNT=0;
RR=sum(PRM.NERE.*PRM.OBSE.^2)+sum(PRM.NERN.*PRM.OBSN.^2)+sum(PRM.NERU.*PRM.OBSU.^2);
while not(COUNT==2)
  RT=RT+1;
  NACC=0;tic
  for iT=1:PRM.CHIN
    [LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM]=RWR_PARM_MCMC(LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,RWD,PRM);
    if PRM.PR==4; SLP.SMP=100.*10^(PRM.AMW.*1.5+9.1)./(PRM.MU.*LEN.SMP.*WID.SMP.*1000.*1000); end
    if DPP==1; DRP.SMP=(2/pi).*PRM.MU.*((SLP.SMP.*0.01)./(WID.SMP*1000)); end % strike slip for surface break dip slip is DRP*4/3
    MW.SMP=(log10(PRM.MU.*SLP.SMP.*LEN.SMP.*WID.SMP.*1000.*1000.*0.01)-9.1)./1.5;
    Rtime=Rtime+1;
    U=DISLOC_RECT_P(PRM.OLAT,PRM.OLON,LAT.SMP(1,:),LON.SMP(1,:),DEP.SMP(1,:),STR.SMP(1,:),DIP.SMP(1,:),RAK.SMP(1,:),LEN.SMP(1,:),WID.SMP(1,:),0.25);
    RES.SMP(1,:)=sum(bsxfun(@times,bsxfun(@plus,-(bsxfun(@times,U.E,SLP.SMP(1,:))+bsxfun(@times,U.Et,TNS.SMP(1,:))),PRM.OBSE).^2,PRM.NERE))+...
                 sum(bsxfun(@times,bsxfun(@plus,-(bsxfun(@times,U.N,SLP.SMP(1,:))+bsxfun(@times,U.Nt,TNS.SMP(1,:))),PRM.OBSN).^2,PRM.NERN))+...
                 sum(bsxfun(@times,bsxfun(@plus,-(bsxfun(@times,U.U,SLP.SMP(1,:))+bsxfun(@times,U.Ut,TNS.SMP(1,:))),PRM.OBSU).^2,PRM.NERU));
%{    
    for NN=1:PRM.NPL
      U=DISLOC_RECT(PRM.OLAT,PRM.OLON,LAT.SMP(1,NN),LON.SMP(1,NN),DEP.SMP(1,NN),STR.SMP(1,NN),DIP.SMP(1,NN),RAK.SMP(1,NN),LEN.SMP(1,NN),WID.SMP(1,NN),0.25);
      RES.SMP(1,NN)=sum(PRM.NERE.*(PRM.OBSE-SLP.SMP(1,NN).*U.E-TNS.SMP(1,NN).*U.Et).^2)+...
                    sum(PRM.NERN.*(PRM.OBSN-SLP.SMP(1,NN).*U.N-TNS.SMP(1,NN).*U.Nt).^2)+...
                    sum(PRM.NERU.*(PRM.OBSU-SLP.SMP(1,NN).*U.U-TNS.SMP(1,NN).*U.Ut).^2);
    end
%}
    if PRM.PR==0
      dPr.SMP=ones(size(RES.SMP));
      PDF=expm1(0.5.*(-RES.SMP+RES.OLD))+1;
    else
      if PRM.PR==1;
        dPr.SMP=(MW.SMP-PRM.MW).^2;
      elseif PRM.PR==2 || PRM.PR==4
        SLIP_V=SLIP_VEC(d2r*STR.SMP,d2r*DIP.SMP,d2r*RAK.SMP);
        SLIP_D1=SLIP_V-SLIP_APR1;
        SLIP_D2=SLIP_V-SLIP_APR2;
        dPr.SMP=min([SLIP_D1'*SLIP_D1;SLIP_D2'*SLIP_D2]);
      elseif PRM.PR==3
        SLIP_V=SLIP_VEC(d2r*STR.SMP,d2r*DIP.SMP,d2r*RAK.SMP);
        SLIP_D1=bsxfun(@times,SLIP_V,MW.SMP)-PRM.AMW.*SLIP_APR1;
        SLIP_D2=bsxfun(@times,SLIP_V,MW.SMP)-PRM.AMW.*SLIP_APR2;
        dPr.SMP=min([SLIP_D1'*SLIP_D1;SLIP_D2'*SLIP_D2]);
      elseif PRM.PR==5;
        [X,Y]=PLTXY(LAT.SMP,LON.SMP,PRM.PLAT(4),PRM.PLON(4));
        dPr.SMP=sqrt(X.^2+Y.^2);
      elseif PRM.PR==6;
        [X,Y]=PLTXY(LAT.SMP,LON.SMP,PRM.PLAT(4),PRM.PLON(4));
        dPr.SMP=sqrt(X.^2+Y.^2+(DEP.SMP-PRM.PDEP(4)).^2);
      elseif PRM.PR==7;
        [X,Y]=PLTXY(LAT.SMP,LON.SMP,PRM.PLAT(4),PRM.PLON(4));
        dPr.SMP=sqrt(X.^2+Y.^2+(DEP.SMP-PRM.PDEP(4)).^2+LEN.SMP.^2+WID.SMP.^2);
      end
%
%   Pdf_post:(1/sqrt(2*pi*exp(L    ))*1/sqrt(2*pi)*exp(-Re2    /2)*exp(-M2    /(2*exp(L    ))))/
%            (1/sqrt(2*pi*exp(L_old))*1/sqrt(2*pi)*exp(-Re2_old/2)*exp(-M2_old/(2*exp(L_old))))
%
      PDF=expm1(0.5.*...
         ((-RES.SMP-LAM.SMP-exp(-LAM.SMP).*dPr.SMP)...
          +(RES.OLD+LAM.OLD+exp(-LAM.OLD).*dPr.OLD)))+1;
    end
%
    IND_M=PDF > rand(1,PRM.NPL,'single');
%
    [RES,dPr,LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,DRP,MW]=UPD_PARM_MCMC(RES,dPr,LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,DRP,MW,IND_M);
    if iT > PRM.CHIN-PRM.KEEP
      [RES,dPr,LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,DRP,MW]=CHA_KEEP_MCMC(iT,PRM,RES,dPr,LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,DRP,MW);
      NACC=NACC+sum(IND_M);
    end
  end
  AJR=NACC./LDIM;
  fprintf('T=%3d MaxRes=%6.3f MinRes=%6.3f Accept=%5.1f RWD=%5.2f Time=%5.1fsec\n',...
           RT,1-max(RES.OLD)./RR,1-min(RES.OLD)./RR,100*AJR,RWD,toc)
  [LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM]=UPD_STPM_MCMC(LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,PRM);
  fprintf('Lat.   Lon.    Depth Str.   Dip   Rake    Len.   Wid.   Slip   Tens.  LAM    Mw    Re\n') 
  fprintf('%6.3f %7.3f %5.2f %6.2f %5.2f %6.2f %6.2f %6.2f %6.1f %6.1f %6.1f %4.2f %6.1f \n',...
  mean(LAT.OLD),mean(LON.OLD),mean(DEP.OLD),mean(STR.OLD),mean(DIP.OLD),mean(RAK.OLD),...
  mean(LEN.OLD),mean(WID.OLD),mean(SLP.OLD),mean(TNS.OLD),mean(LAM.OLD),mean(MW.OLD),min(RES.OLD)) 
  %if PRM.SHOW==1;
  %  figure(100)
  %  plot(PRM.OLON,PRM.OLAT,'s','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','g')
  %  quiver(PRM.OLON,PRM.OLAT,PRM.OBSE,PRM.OBSN)
  %  drown now 
  %  quiver(PRM.OLON,PRM.OLAT,PRM.OBSE,PRM.OBSN)
  %end
  if AJR > 0.24
    RWD=RWD*1.1;
  elseif AJR < 0.22
    RWD=RWD*0.7;
  else
    COUNT=COUNT+1;
  end
  if RT > PRM.RTIM; break; end;
end
SAV_KEEP_MCMC(LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,RES,DRP,MW,NUM,PRM)
end
%% READ PARAMETERS
function PRM=READ_PARM_EstMeca(infile)
Fid=fopen(infile,'r');
PRM.NOBS=fscanf(Fid,'%d \n',[1,1]);
PRM.OLAT=zeros(PRM.NOBS,1);
PRM.OLON=zeros(PRM.NOBS,1);
PRM.OBSE=zeros(PRM.NOBS,1);
PRM.OBSN=zeros(PRM.NOBS,1);
PRM.OBSU=zeros(PRM.NOBS,1);
PRM.ERRN=zeros(PRM.NOBS,1);
PRM.ERRE=zeros(PRM.NOBS,1);
PRM.ERRU=zeros(PRM.NOBS,1);
for N=1:PRM.NOBS
  PRM.OLAT(N,1)=fscanf(Fid,'%f ',[1,1]);
  PRM.OLON(N,1)=fscanf(Fid,'%f ',[1,1]);
  PRM.OBSE(N,1)=fscanf(Fid,'%f ',[1,1]);
  PRM.OBSN(N,1)=fscanf(Fid,'%f ',[1,1]);
  PRM.OBSU(N,1)=fscanf(Fid,'%f ',[1,1]);
  PRM.ERRE(N,1)=fscanf(Fid,'%f ',[1,1]);
  PRM.ERRN(N,1)=fscanf(Fid,'%f ',[1,1]);
  PRM.ERRU(N,1)=fscanf(Fid,'%f ',[1,1]);
end
fgetl(Fid);
for N=1:2
  PRM.ASTR(N)=fscanf(Fid,'%f ',[1,1]);
  PRM.ADIP(N)=fscanf(Fid,'%f ',[1,1]);
  PRM.ARAK(N)=fscanf(Fid,'%f ',[1,1]);
end
PRM.AMW=fscanf(Fid,'%f \n',[1,1]);
fgetl(Fid);
PRM.PLAT=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.PLON=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.PDEP=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.PSTR=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.PDIP=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.PRAK=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.PLEN=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.PWID=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.PSLP=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.PTNS=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.PLAM=fscanf(Fid,'%f %f %f %f \n',[1,4]);
PRM.MU=fscanf(Fid,'%f \n',[1,1]);PRM.MU=PRM.MU.*1e9;
PRM.KEEP=fscanf(Fid,'%f \n',[1,1]);
PRM.CHIN=fscanf(Fid,'%f \n',[1,1]);
PRM.RTIM=fscanf(Fid,'%f \n',[1,1]);
PRM.OUTD=fscanf(Fid,'%s \n',[1,1]);
PRM.PR=fscanf(Fid,'%d \n',[1,1]);
PRM.NPL=fscanf(Fid,'%d \n',[1,1]);
PRM.SHOW=fscanf(Fid,'%d \n',[1,1]);
fclose(Fid);
%-------------------
fprintf('==================\nINPUT DATA\n==================\n') 
fprintf('Num. of SITE: %i \n',PRM.NOBS)
fprintf('  Lat    Lon EW(cm) NS(cm) UD(cm) ERR_EW ERR_NS ERR_UD \n') 
for N=1:PRM.NOBS
  fprintf('%5.2f %5.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',...
      [PRM.OLAT(N) PRM.OLON(N) PRM.OBSN(N) PRM.OBSE(N) PRM.OBSU(N) PRM.ERRE(N) PRM.ERRN(N) PRM.ERRU(N)]) 
end
fprintf('==================\nINPUT PARAMETERS\n==================\n')
fprintf('Range of Lat. LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PLAT) 
fprintf('Range of Lon. LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PLON) 
fprintf('Range of Dep. LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PDEP) 
fprintf('Range of Str. LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PSTR)
fprintf('Range of Dip  LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PDIP)
fprintf('Range of RAKE LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PRAK)
fprintf('Range of LEN. LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PLEN)
fprintf('Range of WID. LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PWID)
fprintf('Range of SLIP LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PSLP)
fprintf('Range of TNS. LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PTNS)
fprintf('Range of Lam. LOWER, UPPPER, WD and Init: %6.1f %6.1f %5.3f %6.1f \n',PRM.PLAM)
fprintf('==================\n')
fprintf(' aPriori  Str.  Dip   Rake \n') 
fprintf(' Plane_1 %5.1f %5.1f %5.1f \n Plane_2 %5.1f %5.1f %5.1f \n',[PRM.ASTR ; PRM.ADIP ; PRM.ARAK]) 
fprintf(' aPriori  Mw : %5.1f \n',PRM.AMW) 
fprintf('==================\n')
fprintf('Mu(GPa), Keep, Chain and Itration Sample: %5.1f %3.1e %3.1e %3.1e \n',[PRM.MU./1e9 PRM.KEEP PRM.CHIN PRM.RTIM]) 
fprintf('OUT_DIR: %s \n',PRM.OUTD) 
fprintf('index_Pr: %i \n',PRM.PR) 
fprintf('Number of Parallel: %i \n',PRM.NPL)
if PRM.SHOW==1;
  figure(100)
  plot(PRM.OLON,PRM.OLAT,'s','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','g')
  quiver(PRM.OLON,PRM.OLAT,PRM.OBSE,PRM.OBSN)
end
end
%% Init PARAMETERS
function [LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,RES,DRP,MW,dPr]=INT_PARM_MCMC(PRM)
LAT.SMP=ones(1,PRM.NPL,'single');
LAT.STD=ones(1,PRM.NPL,'single').*PRM.PLAT(3);
LAT.OLD=ones(1,PRM.NPL,'single').*PRM.PLAT(4);
LAT.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
LON.SMP=ones(1,PRM.NPL,'single');
LON.STD=ones(1,PRM.NPL,'single').*PRM.PLON(3);
LON.OLD=ones(1,PRM.NPL,'single').*PRM.PLON(4);
LON.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
DEP.SMP=ones(1,PRM.NPL,'single');
DEP.STD=ones(1,PRM.NPL,'single').*PRM.PDEP(3);
DEP.OLD=ones(1,PRM.NPL,'single').*PRM.PDEP(4);
DEP.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
STR.SMP=ones(1,PRM.NPL,'single');
STR.STD=ones(1,PRM.NPL,'single').*PRM.PSTR(3);
STR.OLD=ones(1,PRM.NPL,'single').*PRM.PSTR(4);
STR.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
DIP.SMP=ones(1,PRM.NPL,'single');
DIP.STD=ones(1,PRM.NPL,'single').*PRM.PDIP(3);
DIP.OLD=ones(1,PRM.NPL,'single').*PRM.PDIP(4);
DIP.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
RAK.SMP=ones(1,PRM.NPL,'single');
RAK.STD=ones(1,PRM.NPL,'single').*PRM.PRAK(3);
RAK.OLD=ones(1,PRM.NPL,'single').*PRM.PRAK(4);
RAK.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
LEN.SMP=ones(1,PRM.NPL,'single');
LEN.STD=ones(1,PRM.NPL,'single').*PRM.PLEN(3);
LEN.OLD=ones(1,PRM.NPL,'single').*PRM.PLEN(4);
LEN.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
WID.SMP=ones(1,PRM.NPL,'single');
WID.STD=ones(1,PRM.NPL,'single').*PRM.PWID(3);
WID.OLD=ones(1,PRM.NPL,'single').*PRM.PWID(4);
WID.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
SLP.SMP=ones(1,PRM.NPL,'single');
SLP.STD=ones(1,PRM.NPL,'single').*PRM.PSLP(3);
SLP.OLD=ones(1,PRM.NPL,'single').*PRM.PSLP(4);
SLP.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
TNS.SMP=ones(1,PRM.NPL,'single');
TNS.STD=ones(1,PRM.NPL,'single').*PRM.PTNS(3);
TNS.OLD=ones(1,PRM.NPL,'single').*PRM.PTNS(4);
TNS.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
LAM.SMP=ones(1,PRM.NPL,'single');
LAM.STD=ones(1,PRM.NPL,'single').*PRM.PLAM(3);
LAM.OLD=ones(1,PRM.NPL,'single').*PRM.PLAM(4);
LAM.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
RES.SMP=ones(1,PRM.NPL,'single');
RES.OLD= inf(1,PRM.NPL,'single');
RES.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
DRP.OLD= inf(1,PRM.NPL,'single');
DRP.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
dPr.OLD= inf(1,PRM.NPL,'single');
MW.OLD= inf(1,PRM.NPL,'single');
MW.CHA=ones(1,PRM.NPL*PRM.KEEP,'single');
end
%% SAVE PARAMETERS
function SAV_KEEP_MCMC(LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,RES,DRP,MW,NUM,PRM)
KEEP.LLAT=PRM.PLAT(1,1);
KEEP.LLON=PRM.PLON(1,1);
KEEP.LDEP=PRM.PDEP(1,1);
KEEP.LSTR=PRM.PSTR(1,1);
KEEP.LDIP=PRM.PDIP(1,1);
KEEP.LRAK=PRM.PRAK(1,1);
KEEP.LLEN=PRM.PLEN(1,1);
KEEP.LWID=PRM.PWID(1,1);
KEEP.LSLP=PRM.PSLP(1,1);
KEEP.LTNS=PRM.PTNS(1,1);
KEEP.LLAM=PRM.PLAM(1,1);
KEEP.LRES=min(RES.CHA);
KEEP.LDRP=min(DRP.CHA);
KEEP.ILAT=(PRM.PLAT(1,2)-KEEP.LLAT)./double((intmax('uint8')-1));
KEEP.ILON=(PRM.PLON(1,2)-KEEP.LLON)./double((intmax('uint8')-1));
KEEP.IDEP=(PRM.PDEP(1,2)-KEEP.LDEP)./double((intmax('uint8')-1));
KEEP.ISTR=(PRM.PSTR(1,2)-KEEP.LSTR)./double((intmax('uint8')-1));
KEEP.IDIP=(PRM.PDIP(1,2)-KEEP.LDIP)./double((intmax('uint8')-1));
KEEP.IRAK=(PRM.PRAK(1,2)-KEEP.LRAK)./double((intmax('uint8')-1));
KEEP.ILEN=(PRM.PLEN(1,2)-KEEP.LLEN)./double((intmax('uint8')-1));
KEEP.IWID=(PRM.PWID(1,2)-KEEP.LWID)./double((intmax('uint8')-1));
KEEP.ISLP=(PRM.PSLP(1,2)-KEEP.LSLP)./double((intmax('uint8')-1));
KEEP.ITNS=(PRM.PTNS(1,2)-KEEP.LTNS)./double((intmax('uint8')-1));
KEEP.ILAM=(PRM.PLAM(1,2)-KEEP.LLAM)./double((intmax('uint8')-1));
KEEP.IRES=(max(RES.CHA) -KEEP.LRES)./double((intmax('uint8')-1));
KEEP.IDRP=(max(DRP.CHA) -KEEP.LDRP)./double((intmax('uint8')-1));
KEEP.LAT=uint8((LAT.CHA-KEEP.LLAT)./KEEP.ILAT);
KEEP.LON=uint8((LON.CHA-KEEP.LLON)./KEEP.ILON);
KEEP.DEP=uint8((DEP.CHA-KEEP.LDEP)./KEEP.IDEP);
KEEP.STR=uint8((STR.CHA-KEEP.LSTR)./KEEP.ISTR);
KEEP.DIP=uint8((DIP.CHA-KEEP.LDIP)./KEEP.IDIP);
KEEP.RAK=uint8((RAK.CHA-KEEP.LRAK)./KEEP.IRAK);
KEEP.LEN=uint8((LEN.CHA-KEEP.LLEN)./KEEP.ILEN);
KEEP.WID=uint8((WID.CHA-KEEP.LWID)./KEEP.IWID);
KEEP.SLP=uint8((SLP.CHA-KEEP.LSLP)./KEEP.ISLP);
KEEP.TNS=uint8((TNS.CHA-KEEP.LTNS)./KEEP.ITNS);
KEEP.LAM=uint8((LAM.CHA-KEEP.LLAM)./KEEP.ILAM);
KEEP.RES=uint8((RES.CHA-KEEP.LRES)./KEEP.IRES);
KEEP.DRP=uint8((DRP.CHA-KEEP.LDRP)./KEEP.IDRP); %Unit is Pa
KEEP.MW=uint8(MW.CHA.*25);
save_f=fullfile(PRM.OUTD,strcat('KEEP_',num2str(NUM,'%02d')));
save(save_f,'KEEP');
end
%% RANDAM WARK PARAMETERS
function [LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM]=RWR_PARM_MCMC(LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,RWD,PRM)
LAT.SMP=SMP_PRM(PRM.PLAT(1,1),PRM.PLAT(1,2),LAT.OLD,LAT.STD,RWD,PRM);
LON.SMP=SMP_PRM(PRM.PLON(1,1),PRM.PLON(1,2),LON.OLD,LON.STD,RWD,PRM);
DEP.SMP=SMP_PRM(PRM.PDEP(1,1),PRM.PDEP(1,2),DEP.OLD,DEP.STD,RWD,PRM);
STR.SMP=SMP_PRM(PRM.PSTR(1,1),PRM.PSTR(1,2),STR.OLD,STR.STD,RWD,PRM);
DIP.SMP=SMP_PRM(PRM.PDIP(1,1),PRM.PDIP(1,2),DIP.OLD,DIP.STD,RWD,PRM);
RAK.SMP=SMP_PRM(PRM.PRAK(1,1),PRM.PRAK(1,2),RAK.OLD,RAK.STD,RWD,PRM);
LEN.SMP=SMP_PRM(PRM.PLEN(1,1),PRM.PLEN(1,2),LEN.OLD,LEN.STD,RWD,PRM);
WID.SMP=SMP_PRM(PRM.PWID(1,1),PRM.PWID(1,2),WID.OLD,WID.STD,RWD,PRM);
SLP.SMP=SMP_PRM(PRM.PSLP(1,1),PRM.PSLP(1,2),SLP.OLD,SLP.STD,RWD,PRM);
TNS.SMP=SMP_PRM(PRM.PTNS(1,1),PRM.PTNS(1,2),TNS.OLD,TNS.STD,RWD,PRM);
LAM.SMP=SMP_PRM(PRM.PLAM(1,1),PRM.PLAM(1,2),LAM.OLD,LAM.STD,RWD,PRM);
%
IND=STR.SMP>360; STR.SMP(IND)=STR.SMP(IND)-360;
IND=STR.SMP<0;   STR.SMP(IND)=STR.SMP(IND)+360;
%
end
%% RANDAM WARK PARAMETERS SUB.
function SMP=SMP_PRM(LO,UP,OLD,STD,RWD,PRM)
if LO~=UP
  SMP=OLD+RWD.*STD.*(rand(1,PRM.NPL,'single')-0.5);
  IND_S=find(SMP<LO | SMP>UP);
  count=0;
  while isempty(IND_S)==0
    SMP(IND_S)=OLD(IND_S)+RWD.*STD(IND_S).*(rand(1,length(IND_S),'single')-0.5);
    IND_S=find(SMP<LO | SMP>UP);
    if isempty(IND_S)==1; break; end
    if count > 10; 
      SMP(SMP<LO)=LO;
      SMP(SMP>UP)=UP;
      break;
    end
    count=count+1;
  end
else
  SMP=LO.*ones(1,PRM.NPL,'single');
end
end
%% UPDATE PARAMETERS
function [RES,dPr,LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,DRP,MW]=UPD_PARM_MCMC(RES,dPr,LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,DRP,MW,IND)
RES.OLD(IND)=RES.SMP(IND);
dPr.OLD(IND)=dPr.SMP(IND);
LAT.OLD(IND)=LAT.SMP(IND);
LON.OLD(IND)=LON.SMP(IND);
DEP.OLD(IND)=DEP.SMP(IND);
STR.OLD(IND)=STR.SMP(IND);
DIP.OLD(IND)=DIP.SMP(IND);
RAK.OLD(IND)=RAK.SMP(IND);
LEN.OLD(IND)=LEN.SMP(IND);
WID.OLD(IND)=WID.SMP(IND);
SLP.OLD(IND)=SLP.SMP(IND);
TNS.OLD(IND)=TNS.SMP(IND);
LAM.OLD(IND)=LAM.SMP(IND);
DRP.OLD(IND)=DRP.SMP(IND);
MW.OLD(IND) =MW.SMP(IND);
end
%% UPDATE CHAIN PARAMETERS
function [RES,dPr,LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,DRP,MW]=CHA_KEEP_MCMC(iT,PRM,RES,dPr,LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,DRP,MW)
SN=(iT-(PRM.CHIN-PRM.KEEP)-1)*PRM.NPL+1;
EN=(iT-(PRM.CHIN-PRM.KEEP))*PRM.NPL;
RES.CHA(:,SN:EN)=RES.SMP;
dPr.CHA(:,SN:EN)=dPr.SMP;
LAT.CHA(:,SN:EN)=LAT.SMP;
LON.CHA(:,SN:EN)=LON.SMP;
DEP.CHA(:,SN:EN)=DEP.SMP;
STR.CHA(:,SN:EN)=STR.SMP;
DIP.CHA(:,SN:EN)=DIP.SMP;
RAK.CHA(:,SN:EN)=RAK.SMP;
LEN.CHA(:,SN:EN)=LEN.SMP;
WID.CHA(:,SN:EN)=WID.SMP;
SLP.CHA(:,SN:EN)=SLP.SMP;
TNS.CHA(:,SN:EN)=TNS.SMP;
LAM.CHA(:,SN:EN)=LAM.SMP;
DRP.CHA(:,SN:EN)=DRP.SMP;
MW.CHA(:,SN:EN)=MW.SMP;
end
%% UPDATE STD AND PARAMETERS
function [LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM]=UPD_STPM_MCMC(LAT,LON,DEP,STR,DIP,RAK,LEN,WID,SLP,TNS,LAM,PRM)
LAT.STD=repmat(std(LAT.CHA,1,2),1,PRM.NPL);
LON.STD=repmat(std(LON.CHA,1,2),1,PRM.NPL);
DEP.STD=repmat(std(DEP.CHA,1,2),1,PRM.NPL);
STR.STD=repmat(std(STR.CHA,1,2),1,PRM.NPL);
DIP.STD=repmat(std(DIP.CHA,1,2),1,PRM.NPL);
RAK.STD=repmat(std(RAK.CHA,1,2),1,PRM.NPL);
LEN.STD=repmat(std(LEN.CHA,1,2),1,PRM.NPL);
WID.STD=repmat(std(WID.CHA,1,2),1,PRM.NPL);
SLP.STD=repmat(std(SLP.CHA,1,2),1,PRM.NPL);
TNS.STD=repmat(std(TNS.CHA,1,2),1,PRM.NPL);
LAM.OLD=repmat(median(LAM.CHA,2),1,PRM.NPL);
LAT.OLD=repmat(median(LAT.CHA,2),1,PRM.NPL);
LON.OLD=repmat(median(LON.CHA,2),1,PRM.NPL);
DEP.OLD=repmat(median(DEP.CHA,2),1,PRM.NPL);
STR.OLD=repmat(median(STR.CHA,2),1,PRM.NPL);
DIP.OLD=repmat(median(DIP.CHA,2),1,PRM.NPL);
RAK.OLD=repmat(median(RAK.CHA,2),1,PRM.NPL);
LEN.OLD=repmat(median(LEN.CHA,2),1,PRM.NPL);
WID.OLD=repmat(median(WID.CHA,2),1,PRM.NPL);
SLP.OLD=repmat(median(SLP.CHA,2),1,PRM.NPL);
TNS.OLD=repmat(median(TNS.CHA,2),1,PRM.NPL);
LAM.OLD=repmat(median(LAM.CHA,2),1,PRM.NPL);
end
%% MAKE FIGURES
function MAKE_FIGS(KEEP)
%
global OLOC Lo_Pa Up_Pa OBS
%
I_Pa=2*(Up_Pa-Lo_Pa)./single(intmax('uint8')-1);
max_lat=max(save_Eq_Pa(1,:));
max_lon=max(save_Eq_Pa(2,:));
min_lat=min(save_Eq_Pa(1,:));
min_lon=min(save_Eq_Pa(2,:));
int_lat=I_Pa(1);
n_lat=Lo_Pa(1):I_Pa(1):Up_Pa(1);
n_lon=Lo_Pa(2):I_Pa(2):Up_Pa(2);
%
disp('Plotting MAP')
[slon,slat]=meshgrid(n_lon,n_lat);
count=zeros(size(slat));
%
for nlat=1:length(n_lat);
  index=(save_Eq_Pa(1,:) > n_lat(nlat)-int_lat/2 & save_Eq_Pa(1,:) <= n_lat(nlat)+int_lat/2 );
  count(nlat,:)=histc(save_Eq_Pa(2,index),n_lon);
end
tcount=sum(sum(count));
mcount=max(max(count));
%
Fid=fopen('out.txt','w');
fprintf(Fid,'%8.3f %7.3f %6.4e \n',[slon(:) slat(:) count(:)./mcount]');
fclose(Fid);
%
figure('name','Source parameters');
subplot(2,2,1);
contourf(slon,slat,count);
hold on
plot(OLOC(:,2),OLOC(:,1),'s','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','g')
quiver(OLOC(:,2),OLOC(:,1),OBS(:,1),OBS(:,2))
title('Map view');
hold off
%
clear slat slon count
%
subplot(4,4,3);
hist(save_Eq_Pa(1,:),Lo_Pa(1):I_Pa(1):Up_Pa(1));
xlim([min_lat max_lat])
title('Latitude [deg]');
%
subplot(4,4,4);
hist(save_Eq_Pa(2,:),Lo_Pa(2):I_Pa(2):Up_Pa(2));
xlim([min_lon max_lon])
title('Longitude [deg]');
%
disp('Plotting DEPTH')
subplot(4,4,7);
hist(save_Eq_Pa(3,:),Lo_Pa(3):I_Pa(3):Up_Pa(3));
xlim([Lo_Pa(3) Up_Pa(3)])
title('Depth [km]');
%
subplot(4,4,8);
ang=-save_Eq_Pa(4,:).*pi/180+pi/2;
ang(ang>2.*pi)=ang(ang>2.*pi)-2.*pi;
ang(ang<0)=ang(ang<0)+2.*pi;
rose(ang,60);
title('STRIKE');
%
disp('Plotting DIP')
subplot(4,4,9);
hist(save_Eq_Pa(5,:),Lo_Pa(5):I_Pa(5):Up_Pa(5));
xlim([Lo_Pa(5) Up_Pa(5)])
title('DIP [deg]');
%
disp('Plotting RAKE')
subplot(4,4,10);
hist(save_Eq_Pa(6,:),Lo_Pa(6):I_Pa(6):Up_Pa(6));
xlim([Lo_Pa(6) Up_Pa(6)])
title('RAKE [deg]');
%
disp('Plotting LENGTH')
subplot(4,4,11);
hist(save_Eq_Pa(7,:),Lo_Pa(7):I_Pa(7):Up_Pa(7));
xlim([Lo_Pa(7) Up_Pa(7)])
title('LENGTH [km]');
%
disp('Plotting WIDTH')
subplot(4,4,12);
hist(save_Eq_Pa(8,:),Lo_Pa(8):I_Pa(8):Up_Pa(8));
xlim([Lo_Pa(8) Up_Pa(8)])
title('WIDTH [km]');
%
disp('Plotting SLIP')
subplot(4,4,13);
hist(save_Eq_Pa(9,:),Lo_Pa(9):I_Pa(9):Up_Pa(9));
xlim([Lo_Pa(9) Up_Pa(9)])
title('SLIP [cm]');
%
subplot(4,4,14);
hist(save_Eq_Pa(10,:),Lo_Pa(10):I_Pa(10):Up_Pa(10));
xlim([Lo_Pa(10) Up_Pa(10)])
title('log(Alpha) for APRIORI');
%
subplot(4,4,15);
%{
max_Mw=max(save_Eq_Pa(11,:));
min_Mw=min(save_Eq_Pa(11,:));
%max_Mw=7;
%min_Mw=6;
hist(save_Eq_Pa(11,:),min_Mw:1/25:max_Mw);
xlim([min_Mw max_Mw])
title('Mw');
%}
max_Dp=5;
hist(save_Eq_Pa(13,:)/1e6,0:0.01*max_Dp:max_Dp+0.01*max_Dp);
xlim([0 max_Dp])
title('Stress Drop (MPa)');
%
subplot(4,4,16);
nlen=length(save_Eq_Pa);
nint=fix(nlen./1000);
plot(1:nint:nlen,save_Eq_Pa(12,1:nint:nlen),'.');
min_Re=min(save_Eq_Pa(12,:));
xlim([1 nlen])
ylim([min_Re median(save_Eq_Pa(12,:))])
title('Residual^2');
%
end
%% READ KEEP
function READ_KEEP(infile)
PRM=READ_PARM_EstMeca(infile);
%
I_Pa=single(PRM.Up-PRM.Lo)./single(intmax('uint8')-1);
Max_Nf=9;
dat_dir='./dat_save';
save_Eq_Pa=zeros(13,Keep.*Max_Nf,'single');
IM_Pa=repmat(I_Pa',1,Keep);
LoM_Pa=single(repmat(Lo_Pa',1,Keep));
for Nf=1:Max_Nf
  load_f=fullfile(PRM.OUTD,strcat('KEEP_',num2str(NUM,'%02d')));
  disp(['Reading ',load_f])
  load(load_f);
  save_Eq_Pa(:,((Nf-1)*Keep+1):Nf*Keep)=[IM_Pa.*single(chain_Pa)+LoM_Pa-IM_Pa/10; single(chain_Mw)./25; single(chain_Re)./10; single(chain_Dp)];
end
save('save_dat.mat','save_Eq_Pa')
MAKE_FIGS(save_Eq_Pa);
end