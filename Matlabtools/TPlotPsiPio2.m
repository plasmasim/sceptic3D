function [] = TPlotPsiPio2(filename)

% Script plotting the temperature tensor color plots in either cartesian or
% spherical coordinates.

opengl neverselect


short=false;readforce=false;cont=1;
readoutput();

% Actual center of the first and last theta cells
if(tcc(1)==1)
    tcc(1)=0.25*(3+tcc(2));
end
if(tcc(nthused)==-1)
    tcc(nthused)=0.25*(-3+tcc(nthused-1));
end

for i=1:npsiused
    pcc(i)=0.+double(i-1)*2*pi/double(npsiused);
end

ptopm=ceil((npsiused-1)/4)+1;
pbotm=ceil(3*(npsiused-1)/4)+1;

% Build the temperature tensor in spherical coordinates
Tr=(vr2sum.*psum-power(vrsum,2))./power(psum,2);
Tt=(vt2sum.*psum-power(vtsum,2))./power(psum,2);
Tp=(vp2sum.*psum-power(vpsum,2))./power(psum,2);
Trt=(vrtsum.*psum-vrsum.*vtsum)./power(psum,2);
Trp=(vrpsum.*psum-vrsum.*vpsum)./power(psum,2);
Ttp=(vtpsum.*psum-vtsum.*vpsum)./power(psum,2);

% Get the temperature terms on the major cross-section
TrT=double(npsiused)/(2*pi)*(Tr(:,:,ptopm)*(pcc(ptopm+1)-pi/2)+Tr(:,:,ptopm+1)*(pi/2-pcc(ptopm)));
TrB=double(npsiused)/(2*pi)*(Tr(:,:,pbotm)*(pcc(pbotm+1)-3*pi/2)+Tr(:,:,pbotm+1)*(3*pi/2-pcc(pbotm)));
TtT=double(npsiused)/(2*pi)*(Tt(:,:,ptopm)*(pcc(ptopm+1)-pi/2)+Tt(:,:,ptopm+1)*(pi/2-pcc(ptopm)));
TtB=double(npsiused)/(2*pi)*(Tt(:,:,pbotm)*(pcc(pbotm+1)-3*pi/2)+Tt(:,:,pbotm+1)*(3*pi/2-pcc(pbotm)));
TpT=double(npsiused)/(2*pi)*(Tp(:,:,ptopm)*(pcc(ptopm+1)-pi/2)+Tp(:,:,ptopm+1)*(pi/2-pcc(ptopm)));
TpB=double(npsiused)/(2*pi)*(Tp(:,:,pbotm)*(pcc(pbotm+1)-3*pi/2)+Tp(:,:,pbotm+1)*(3*pi/2-pcc(pbotm)));

TrtT=double(npsiused)/(2*pi)*(Trt(:,:,ptopm)*(pcc(ptopm+1)-pi/2)+Trt(:,:,ptopm+1)*(pi/2-pcc(ptopm)));
TrtB=double(npsiused)/(2*pi)*(Trt(:,:,pbotm)*(pcc(pbotm+1)-3*pi/2)+Trt(:,:,pbotm+1)*(3*pi/2-pcc(pbotm)));
TrpT=double(npsiused)/(2*pi)*(Trp(:,:,ptopm)*(pcc(ptopm+1)-pi/2)+Trp(:,:,ptopm+1)*(pi/2-pcc(ptopm)));
TrpB=double(npsiused)/(2*pi)*(Trp(:,:,pbotm)*(pcc(pbotm+1)-3*pi/2)+Trp(:,:,pbotm+1)*(3*pi/2-pcc(pbotm)));
TtpT=double(npsiused)/(2*pi)*(Ttp(:,:,ptopm)*(pcc(ptopm+1)-pi/2)+Ttp(:,:,ptopm+1)*(pi/2-pcc(ptopm)));
TtpB=double(npsiused)/(2*pi)*(Ttp(:,:,pbotm)*(pcc(pbotm+1)-3*pi/2)+Ttp(:,:,pbotm+1)*(3*pi/2-pcc(pbotm)));

% liberate some memory
clear *sum

Tcc=[1;tcc;-1];
[R,T]=meshgrid(rcc,Tcc);

XR=R.*T;
YR=R.*sqrt(1-T.*T);
Dtheta1=acos(0.25*(3+tcc(2)));
Dtheta2=acos(tcc(2));
Dr=(rcc(nrused)-1)./(nrused-1);

prec=1;
[r,the]=meshgrid(1:Dr*prec:rcc(nrused)+Dr/nrused,0:pi/round(pi/(prec*Dtheta1)):pi);
xr=r.*cos(the);
yr=r.*sin(the);

% Average the temperature on axis
TrT=[0.5*(TrT(:,1)+TrB(:,1)) TrT 0.5*(TrT(:,nthused)+TrB(:,nthused))];
TrB=[TrT(:,1) TrB TrT(:,nthused+2)];
TtT=[0.5*(TtT(:,1)+TtB(:,1)) TtT 0.5*(TtT(:,nthused)+TtB(:,nthused))];
TtB=[TtT(:,1) TtB TtT(:,nthused+2)];
TpT=[0.5*(TpT(:,1)+TpB(:,1)) TpT 0.5*(TpT(:,nthused)+TpB(:,nthused))];
TpB=[TpT(:,1) TpB TpT(:,nthused+2)];

% The cross-field temperatures involving angles are averaged in order to
% make sense on the upper portion of the domain, hence the off-diagonal
% terms involving angles need a minus sign.
TrtT=[0.5*(TrtT(:,1)-TrtB(:,1)) TrtT 0.5*(TrtT(:,nthused)-TrtB(:,nthused))];
TrtB=[TrtT(:,1) -TrtB TrtT(:,nthused+2)];
TrpT=[0.5*(TrpT(:,1)-TrpB(:,1)) TrpT 0.5*(TrpT(:,nthused)-TrpB(:,nthused))];
TrpB=[TrpT(:,1) -TrpB TrpT(:,nthused+2)];
TtpT=[0.5*(TtpT(:,1)+TtpB(:,1)) TtpT 0.5*(TtpT(:,nthused)+TtpB(:,nthused))];
TtpB=[TtpT(:,1) TtpB TtpT(:,nthused+2)];

% Interpolate the temperature on a new grid with equal theta spacing, to
% occupy less memory
TrTr=griddata(XR,YR,TrT',xr,yr);
TrBr=griddata(XR,YR,TrB',xr,yr);
TtTr=griddata(XR,YR,TtT',xr,yr);
TtBr=griddata(XR,YR,TtB',xr,yr);
TpTr=griddata(XR,YR,TpT',xr,yr);
TpBr=griddata(XR,YR,TpB',xr,yr);

TrtTr=griddata(XR,YR,TrtT',xr,yr);
TrtBr=griddata(XR,YR,TrtB',xr,yr);
TrpTr=griddata(XR,YR,TrpT',xr,yr);
TrpBr=griddata(XR,YR,TrpB',xr,yr);
TtpTr=griddata(XR,YR,TtpT',xr,yr);
TtpBr=griddata(XR,YR,TtpB',xr,yr);

% Rotation of the Temperature tensor to have it on the xy axis
    nthe=size(the,1);
    nrr=size(the,2);
 %   for j=1:nthe
 %       R=[cos(the(j,1)) sin(the(j,1));-sin(the(j,1)) cos(the(j,1))];
 %       for k=1:nrr
 %           TpolT=[TrTr(j,k) TrtTr(j,k);TrtTr(j,k) TtTr(j,k)];
 %           TpolB=[TrBr(j,k) TrtBr(j,k);TrtBr(j,k) TtBr(j,k)];
 %           TcarT=R*TpolT*R^-1;
 %           TcarB=R*TpolB*R^-1;
 %           TzTr(j,k)=TcarT(1,1);
 %           TyTr(j,k)=TcarT(2,2);
 %           TzyTr(j,k)=TcarT(2,1);
 %           TzBr(j,k)=TcarB(1,1);
 %           TyBr(j,k)=TcarB(2,2);
 %           TzyBr(j,k)=TcarB(2,1);
 %       end
 %   end
    
     for j=1:nthe
        R=[cos(the(j,1)) sin(the(j,1)) 0;-sin(the(j,1)) cos(the(j,1)) 0;0 0 1];
        for k=1:nrr
            TpolT=[TrTr(j,k) TrtTr(j,k) TrpTr(j,k);TrtTr(j,k) TtTr(j,k) TtpTr(j,k);TrpTr(j,k) TtpTr(j,k) TpTr(j,k)];
            TpolB=[TrBr(j,k) TrtBr(j,k) TrpBr(j,k);TrtBr(j,k) TtBr(j,k) TtpBr(j,k);TrpBr(j,k) TtpBr(j,k) TpBr(j,k)];
            TcarT=R^-1*TpolT*R;
            TcarB=R*TpolB*R^-1;
            
            
            TzTr(j,k)=TcarT(1,1);
            TyTr(j,k)=TcarT(2,2);
            TxTr(j,k)=TcarT(3,3);
            TzyTr(j,k)=TcarT(2,1);
            % e_x is antiparallel to e_psi in the upper portion of the
            % domain
            TzxTr(j,k)=-TcarT(3,1);
            TyxTr(j,k)=-TcarT(3,2);
             
            TzBr(j,k)=TcarB(1,1);
            TyBr(j,k)=TcarB(2,2);
            TxBr(j,k)=TcarB(3,3);
            TzyBr(j,k)=TcarB(2,1);
            TzxBr(j,k)=-TcarB(3,1);
            TyxBr(j,k)=-TcarB(3,2);
        end
    end
    
    
sinus=sin([0:pi/32:2*pi]);
cosinus=cos([0:pi/32:2*pi]);


%%%%%%%%%%%%%%%%%%%%%
% Start contouring %%
%%%%%%%%%%%%%%%%%%%%%

% Cartesian
if(1)

    
% Tz
figure;hold all;box on

plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
%pcolor(xr,yr,TzTr)
%pcolor(xr,-yr,TzBr)

[C,h]=contour(xr,yr,TzTr,Ti*[.5 .7 .9 1.1 1.5 2 3],'LineWidth',1,'LabelSpacing',2000);
text_handle=clabel(C,h);
set(text_handle,'FontSize',13);
[C,h]=contour(xr,-yr,TzBr,Ti*[.5 .7 .9 1.1 1.5 2 3],'LineWidth',1,'LabelSpacing',2000);
text_handle=clabel(C,h);
set(text_handle,'FontSize',13);

if(0)
shading interp
cmin=0.1*round(10*min(min(min(TzTr)),min(min(TzBr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TzTr)),max(max(TzBr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,Ti);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Tz','FontSize',22);

% Tz Contours
if(cont)
%[C,h]=contour(xr,yr,TzTr,Ti*[0.1 0.2 0.3 0.5 0.7 0.9 1.1 1.5 2 3],'LineWidth',1,'LabelSpacing',200);
%text_handle=clabel(C,h);
%set(text_handle,'FontSize',13);
if(1)
C=contour(xr,yr,TzTr,Ti*[0.1 0.2 0.3 0.5 0.7 0.9 1.1 1.5 2 3]);
Clength=size(C,2);
 temp1=1;
 while(true)
 	temp2=C(2,temp1)+temp1;
  	if(temp2==Clength)
        plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k-','LineWidth',1)
        break;
    end
  	plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k-','LineWidth',1)
 	temp1=temp2+1;
 end
end
%[C,h]=contour(xr,-yr,TzBr,Ti*[0.1 0.2 0.3 0.5 0.7 0.9 1.1 1.5 2 3],'LineWidth',1,'LabelSpacing',200);
%text_handle=clabel(C,h);
%set(text_handle,'FontSize',13);
if(1)
C=contour(xr,-yr,TzBr,Ti*[0.1 0.2 0.3 0.5 0.7 0.9 1.1 1.5 2 3]);
Clength=size(C,2);
temp1=1;
while(true)
    temp2=C(2,temp1)+temp1;
    if(temp2==Clength)
        plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k-','LineWidth',1)
        break;
    end
    plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k-','LineWidth',1)
    temp1=temp2+1;
end
end
end
end

% Ty
figure;hold all;box on

plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
%pcolor(xr,yr,TyTr)
%pcolor(xr,-yr,TyBr)

[C,h]=contour(xr,yr,TyTr,Ti*[.5 .7 .9 1.05  1.2 1.5 2 3],'LineWidth',1,'LabelSpacing',2000);
text_handle=clabel(C,h);
set(text_handle,'FontSize',13);
[C,h]=contour(xr,-yr,TyBr,Ti*[.5 .7 .9 1.05  1.2 1.5 2 3],'LineWidth',1,'LabelSpacing',2000);
text_handle=clabel(C,h);
set(text_handle,'FontSize',13);
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Ty','FontSize',22);

if(0)
shading interp
cmin=0.1*round(10*min(min(min(TyTr)),min(min(TyBr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TyTr)),max(max(TyBr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,Ti);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Ty','FontSize',22);

% Ty Contours
if(cont)
[C,h]=contour(xr,yr,TyTr,Ti*[0.1 0.2 0.3 0.5 0.7 0.9 1.1 1.5 2],'LineWidth',1,'LabelSpacing',200);
text_handle=clabel(C,h);
set(text_handle,'FontSize',13);
 C=contour(xr,yr,TyTr,Ti*[0.1 0.2 0.3 0.5 0.7 0.9 1.1 1.5 2]);
 Clength=size(C,2);
 temp1=1;
 while(true)
 	temp2=C(2,temp1)+temp1;
  	if(temp2==Clength)
        plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k-','LineWidth',1)
        break;
    end
  	plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k-','LineWidth',1)
 	temp1=temp2+1;
 end
[C,h]=contour(xr,-yr,TyBr,Ti*[0.1 0.2 0.3 0.5 0.7 0.9 1.1 1.5 2],'LineWidth',1,'LabelSpacing',200);
text_handle=clabel(C,h);
set(text_handle,'FontSize',13);
C=contour(xr,-yr,TyBr,Ti*[0.1 0.2 0.3 0.5 0.7 0.9 1.1 1.5 2]);
Clength=size(C,2);
temp1=1;
while(true)
    temp2=C(2,temp1)+temp1;
    if(temp2==Clength)
        plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k-','LineWidth',1)
        break;
    end
    plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k-','LineWidth',1)
    temp1=temp2+1;
end
end
end


% Tx
figure;hold all;box on

plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
%pcolor(xr,yr,TxTr)
%pcolor(xr,-yr,TxBr)


[C,h]=contour(xr,yr,TxTr,Ti*[.5 .7 .9 1.05 1.2 1.5 2 3],'LineWidth',1,'LabelSpacing',2000);
text_handle=clabel(C,h);
set(text_handle,'FontSize',13);
[C,h]=contour(xr,-yr,TxBr,Ti*[.5 .7 .9 1.05 1.2 1.5 2 3],'LineWidth',1,'LabelSpacing',2000);
text_handle=clabel(C,h);
set(text_handle,'FontSize',13);
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Tx','FontSize',22);

if(0)
shading interp
cmin=0.1*round(10*min(min(min(TxTr)),min(min(TxBr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TxTr)),max(max(TxBr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,Ti);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Tx','FontSize',22);

% Tzy
figure;hold all;box on

plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,TzyTr)
pcolor(xr,-yr,TzyBr)

%[C,h]=contour(xr,yr,TzyTr,[0.1*Ti,0.3*Ti,Ti,3*Ti,10*Ti],'LineWidth',1,'LabelSpacing',2000);
%text_handle=clabel(C,h);
%set(text_handle,'FontSize',13);
%[C,h]=contour(xr,-yr,TzyBr,[0.1*Ti,0.3*Ti,Ti,3*Ti,10*Ti],'LineWidth',1,'LabelSpacing',2000);
%text_handle=clabel(C,h);
%set(text_handle,'FontSize',13);

shading interp
cmin=0.1*round(10*min(min(min(TzyTr)),min(min(TzyBr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TzyTr)),max(max(TzyBr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,0);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Tzy','FontSize',22);
end

% Tzx
figure;hold all;box on

plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
pcolor(xr,yr,TzxTr)
pcolor(xr,-yr,TzxBr)

shading interp
cmin=0.1*round(10*min(min(min(TzxTr)),min(min(TzxBr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TzxTr)),max(max(TzxBr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,0);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Tzx','FontSize',22);

% Tyx
figure;hold all;box on

plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
pcolor(xr,yr,TyxTr)
pcolor(xr,-yr,TyxBr)

shading interp
cmin=0.1*round(10*min(min(min(TyxTr)),min(min(TyxBr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TyxTr)),max(max(TyxBr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,0);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Tyx','FontSize',22);

end

% Spherical
if(0)

figure;hold all;box on

%contour(xr,yr,TrTr);
%contour(xr,-yr,TrBr);

plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,TrTr)
pcolor(xr,-yr,TrBr)
shading interp
cmin=0.1*round(10*min(min(min(TrTr)),min(min(TrBr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TrTr)),max(max(TrBr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,4*Ti);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')

axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Tr','FontSize',22);

figure;hold all;box on

%contour(xr,yr,TtTr);
%contour(xr,-yr,TtBr);

plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,TtTr)
pcolor(xr,-yr,TtBr)

shading interp
cmin=0.1*round(10*min(min(min(TtTr)),min(min(TtBr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TtTr)),max(max(TtBr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,Ti);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')

axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Tt','FontSize',22);

figure;hold all;box on

%contour(xr,yr,TrtTr);
%contour(xr,-yr,TrtBr);

plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,TrtTr)
pcolor(xr,-yr,TrtBr)
shading interp
cmin=0.1*round(10*min(min(min(TrtTr)),min(min(TrtBr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TrtTr)),max(max(TrtBr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,Ti);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')

axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);
title('Trt','FontSize',22);


end




