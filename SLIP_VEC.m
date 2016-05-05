%% CONVERT SLIP VECTER FROM FAULT PARAMETERS
% INPUT 
% str : STRIKE (RADIAN) from Norh clock wise
% dip : DIP (RARIAN)
% rake: RAKE (RADIAN)
% OUTPUT
% slip_v : SLIP VECTOR
%
% code by T.Ito 2013/11/29
% modified by T.Ito 2015/05/03
function slip_v=SLIP_VEC(str,dip,rake)
cos_s=cos(str);
sin_s=sin(str);
cos_d=cos(dip);
sin_d=sin(dip);
cos_r=cos(rake);
sin_r=sin(rake);
slip_v=[ cos_r.*cos_s+sin_r.*cos_d.*sin_s;...
         cos_r.*sin_s-sin_r.*cos_d.*cos_s;...
         -sin_r.*sin_d];
end