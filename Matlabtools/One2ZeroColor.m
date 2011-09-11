function [cmin,cmax] = One2ZeroColor(cmin,cmax,cwhite)

% This script modifies the jet colormap in order to go through white when
% the density is one. This is to make sure that most of the colorplot is
% white rather than colored.

if(abs(cmin-cmax)<0.1);
    cmin=cmin-0.05;
    cmax=cmax+0.05;
end

hh=jet(128);
mid=64;


hh(mid,:)=[1 1 1];
delta=24;
for k=1:3
    hhm=hh(mid-delta,k);
    hhp=hh(mid+delta,k);
    hh(mid-delta:mid,k)=linspace(hhm,hh(mid,k),delta+1);
    hh(mid:mid+delta,k)=linspace(hh(mid,k),hhp,delta+1);
end

neutr=round(128*(cwhite-cmin)/double(cmax-cmin));

for k=1:neutr
    hnew(k,:)=hh(round(k*mid/neutr),:);
end
for k=neutr+1:128
    hnew(k,:)=hh(round(mid+mid*(k-neutr)/(128-neutr)),:);
end


colormap(hnew);

end