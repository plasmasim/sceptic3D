function [] = TPlotTh0(filename)

% Script plotting the temperature tensor color plots in either cartesian or
% spherical coordinates.

opengl neverselect


short=false;
readforce=false;
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


ttopm=floor(nthused/2);

% Build the temperature tensor in spherical coordinates
Tr=(vr2sum.*psum-power(vrsum,2))./power(psum,2);
Tt=(vt2sum.*psum-power(vtsum,2))./power(psum,2);
Tp=(vp2sum.*psum-power(vpsum,2))./power(psum,2);
Trt=(vrtsum.*psum-vrsum.*vtsum)./power(psum,2);
Trp=(vrpsum.*psum-vrsum.*vpsum)./power(psum,2);
Ttp=(vtpsum.*psum-vtsum.*vpsum)./power(psum,2);

% Get the temperature terms on the major cross-section
TrT=double(nthused-1)/(2)*(Tr(:,ttopm,:)*(0-tcc(ttopm+1))+Tr(:,ttopm+1,:)*(tcc(ttopm)-0));
TtT=double(nthused-1)/(2)*(Tt(:,ttopm,:)*(0-tcc(ttopm+1))+Tt(:,ttopm+1,:)*(tcc(ttopm)-0));
TpT=double(nthused-1)/(2)*(Tp(:,ttopm,:)*(0-tcc(ttopm+1))+Tp(:,ttopm+1,:)*(tcc(ttopm)-0));

TrtT=double(nthused-1)/(2)*(Trt(:,ttopm,:)*(0-tcc(ttopm+1))+Trt(:,ttopm+1,:)*(tcc(ttopm)-0));
TrpT=double(nthused-1)/(2)*(Trp(:,ttopm,:)*(0-tcc(ttopm+1))+Trp(:,ttopm+1,:)*(tcc(ttopm)-0));
TtpT=double(nthused-1)/(2)*(Ttp(:,ttopm,:)*(0-tcc(ttopm+1))+Ttp(:,ttopm+1,:)*(tcc(ttopm)-0));

Tr=reshape(TrT,nrused,npsiused);
Tt=reshape(TtT,nrused,npsiused);
Tp=reshape(TpT,nrused,npsiused);
Trt=reshape(TrtT,nrused,npsiused);
Trp=reshape(TrpT,nrused,npsiused);
Ttp=reshape(TtpT,nrused,npsiused);

Tr=[Tr(:,:) Tr(:,1)];
Tt=[Tt(:,:) Tt(:,1)];
Tp=[Tp(:,:) Tp(:,1)];
Trt=[Trt(:,:) Trt(:,1)];
Trp=[Trp(:,:) Trp(:,1)];
Ttp=[Ttp(:,:) Ttp(:,1)];

% liberate some memory
clear *sum

prec=1;
Pcc=[pcc;0];
[R,T]=meshgrid(rcc,Pcc);
XR=R.*cos(T);
YR=R.*sin(T);

[r,the]=meshgrid(linspace(1,rcc(nrused),80),linspace(0,2*pi,40));
xr=r.*cos(the);
yr=r.*sin(the);

TrTr=griddata(XR,YR,Tr',xr,yr);
TtTr=griddata(XR,YR,Tt',xr,yr);
TpTr=griddata(XR,YR,Tp',xr,yr);

TrtTr=griddata(XR,YR,Trt',xr,yr);
TrpTr=griddata(XR,YR,Trp',xr,yr);
TtpTr=griddata(XR,YR,Ttp',xr,yr);

% Rotation of the Temperature tensor to have it on the xyz axis
     for j=1:40
        R=[cos(the(j,1)) sin(the(j,1)) 0;-sin(the(j,1)) cos(the(j,1)) 0;0 0 1];
        for k=1:80
            TpolT=[TrTr(j,k) TrpTr(j,k) TrtTr(j,k);TrpTr(j,k) TpTr(j,k) TtpTr(j,k);TrtTr(j,k) TtpTr(j,k) TtTr(j,k)];
            TcarT=R^-1*TpolT*R;
            
            TxTr(j,k)=TcarT(1,1);
            TyTr(j,k)=TcarT(2,2);
            TzTr(j,k)=TcarT(3,3);
            TyxTr(j,k)=TcarT(2,1);
            TzxTr(j,k)=-TcarT(3,1);
            TzyTr(j,k)=-TcarT(3,2);
        end
    end


%%%%%%%%%%%%%%%%%%%%%
% Start contouring %%
%%%%%%%%%%%%%%%%%%%%%
   
sinus=sin([0:pi/32:2*pi]);
cosinus=cos([0:pi/32:2*pi]);

% Tz
figure;hold all; box on;
plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,TzTr)
hold all
shading interp
cmin=0.1*round(10*min(min(min(TzTr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TzTr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,Ti);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('NorthOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('x','FontSize',22);
title('Tz','FontSize',22);

% Ty
figure;hold all; box on;
plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,TyTr)
hold all
shading interp
cmin=0.1*round(10*min(min(min(TyTr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TyTr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,Ti);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('NorthOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('x','FontSize',22);
title('Ty','FontSize',22);


% Tx
figure;hold all; box on;
plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,TxTr)
hold all
shading interp
cmin=0.1*round(10*min(min(min(TxTr))));
cmax=min(5*Ti,0.1*round(10*max(max(max(TxTr)))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,Ti);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('NorthOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('x','FontSize',22);
title('Tx','FontSize',22);

% Tzy
figure;hold all; box on;
plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,TzyTr)
hold all
shading interp
cmin=0.1*round(10*min(min(min(TzyTr))));
cmax=0.1*round(10*max(max(max(TzyTr))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,0);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('NorthOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('x','FontSize',22);
title('Tzy','FontSize',22);


% Tzx
figure;hold all; box on;
plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,TzxTr)
hold all
shading interp
cmin=0.1*round(10*min(min(min(TzxTr))));
cmax=0.1*round(10*max(max(max(TzxTr))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,0);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('NorthOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('x','FontSize',22);
title('Tzx','FontSize',22);

% Tyx
figure;hold all; box on;
plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,TyxTr)
hold all
shading interp
cmin=0.1*round(10*min(min(min(TyxTr))));
cmax=0.1*round(10*max(max(max(TyxTr))));
%colormap(hot);
[cmin,cmax]=One2ZeroColor(cmin,cmax,0);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('NorthOutside')
axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('x','FontSize',22);
title('Tyx','FontSize',22);
    
end




