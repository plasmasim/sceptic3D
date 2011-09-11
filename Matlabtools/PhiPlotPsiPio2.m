function [] = PhiPlotPsiPio2(filename)


%function [] = PhiPlotPsiPio2(filename)
%
% Color-plot of potential contours on the plane {0,e_z,e_y}. Also
% shows the fluid strealines.
%
% In order to have white color when the potential
% is unperturbed (phi=0), I call the function function
% One2ZeroColor to change the colormap.
short=false;readforce=false;
opengl neverselect

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read SCEPTIC output file
%%%%%%%%%%%%%%%%%%%%%%%%%%

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

PhiT=double(npsiused)/(2*pi)*(phi(:,:,ptopm)*(pcc(ptopm+1)-pi/2)...
    +phi(:,:,ptopm+1)*(pi/2-pcc(ptopm)));
PhiB=double(npsiused)/(2*pi)*(phi(:,:,pbotm)*(pcc(pbotm+1)-3*pi/2)...
    +phi(:,:,pbotm+1)*(3*pi/2-pcc(pbotm)));

vrT=double(npsiused)/(2*pi)*(vrsum(:,:,ptopm)*(pcc(ptopm+1)-pi/2)...
    +vrsum(:,:,ptopm+1)*(pi/2-pcc(ptopm)));
vrB=double(npsiused)/(2*pi)*(vrsum(:,:,pbotm)*(pcc(pbotm+1)-3*pi/2)...
    +vrsum(:,:,pbotm+1)*(3*pi/2-pcc(pbotm)));

vtT=double(npsiused)/(2*pi)*(vtsum(:,:,ptopm)*(pcc(ptopm+1)-pi/2)...
    +vtsum(:,:,ptopm+1)*(pi/2-pcc(ptopm)));
vtB=double(npsiused)/(2*pi)*(vtsum(:,:,pbotm)*(pcc(pbotm+1)-3*pi/2)...
    +vtsum(:,:,pbotm+1)*(3*pi/2-pcc(pbotm)));

psumT=double(npsiused)/(2*pi)*(psum(:,:,ptopm)*(pcc(ptopm+1)-pi/2)...
    +psum(:,:,ptopm+1)*(pi/2-pcc(ptopm)));
psumB=double(npsiused)/(2*pi)*(psum(:,:,pbotm)*(pcc(pbotm+1)-3*pi/2)...
    +psum(:,:,pbotm+1)*(3*pi/2-pcc(pbotm)));

% Interpolate values on axis

PhiT=[0.5*(PhiT(:,1)+PhiB(:,1)) PhiT 0.5*(PhiT(:,nthused)+PhiB(:,nthused))];
PhiB=[PhiT(:,1) PhiB PhiT(:,nthused+2)];

psumT=[0.5*(psumT(:,1)+psumB(:,1)) psumT 0.5*(psumT(:,nthused)+psumB(:,nthused))];
psumB=[psumT(:,1) psumB psumT(:,nthused+2)];


vrT=[0.5*(vrT(:,1)+vrB(:,1)) vrT 0.5*(vrT(:,nthused)+vrB(:,nthused))];
vrB=[vrT(:,1) vrB vrT(:,nthused+2)];

% There was a sign error when averaging the theta velocities
vtT=[0.5*(vtT(:,1)-vtB(:,1)) vtT 0.5*(vtT(:,nthused)-vtB(:,nthused))];
vtB=[-vtT(:,1) vtB -vtT(:,nthused+2)];


Tcc=[1;tcc;-1];


[R,T]=meshgrid(rcc,Tcc);

XR=R.*T;
YR=R.*sqrt(1-T.*T);

Dtheta1=acos(0.25*(3+tcc(2)));
Dtheta2=acos(tcc(2));
Dr=(rcc(nrused)-1)./(nrused-1);

%%%%%%%%%%%%%%%%%%%%%
% Start contouring %%
%%%%%%%%%%%%%%%%%%%%%

[r,the]=meshgrid(linspace(1,rcc(nrused),40),linspace(0,pi,20));

%xr=r.*cos(the);
%yr=r.*sin(the);
%PhiTr=griddata(XR,YR,PhiT',xr,yr);
%PhiBr=griddata(XR,YR,PhiB',xr,yr);

PhiTr=griddata(R,T,PhiT',r,cos(the));
PhiBr=griddata(R,T,PhiB',r,cos(the));
xr=r.*cos(the);
yr=r.*sin(the);

BaxR=[-rcc(nrused) 0 rcc(nrused)]*1.1;
sinus=sin([0:pi/32:2*pi]);
cosinus=cos([0:pi/32:2*pi]);

figure;hold all
%plot(BaxR*cB,BaxR*sqrt(1-cB^2),'r-','LineWidth',2)
%fill(cosinus,sinus,[0.9 0.9 0.9]);
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,PhiTr)
hold all
pcolor(xr,-yr,PhiBr)

shading interp;
cmin=0.1*round(10*min(min(min(PhiBr)),min(min(PhiTr))));
cmax=max(0,0.1*round(10*max(max(max(PhiBr)),max(max(PhiTr)))));
%colormap(hot);
One2ZeroColor(cmin,cmax,0);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')

axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('z','FontSize',22);

% Now, retake the dat with more precision for the contours and velocity plots
[r,the]=meshgrid(linspace(1,rcc(nrused),nrused),linspace(0,pi,nthused));


%xr=r.*cos(the);
%yr=r.*sin(the);
%PhiTr=griddata(XR,YR,PhiT',xr,yr);
%PhiBr=griddata(XR,YR,PhiB',xr,yr);

PhiTr=griddata(R,T,PhiT',r,cos(the));
PhiBr=griddata(R,T,PhiB',r,cos(the));
xr=r.*cos(the);
yr=r.*sin(the);


%
% Do the contours by hand, to have them all back
%

if(0)
    figure;hold all;
    
    cosinus=cos([0:pi/32:2*pi]);
    plot(cosinus,sinus,'k')
    plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
    fill(cosinus,sinus,[.85 .85 .85]);
    
    
    [C,h]=contour(xr,yr,PhiTr,[-10 -8 -3  -1 -.3 -.1 -.03],'LineWidth',1,'LabelSpacing',1000);
    text_handle=clabel(C,h);
    set(text_handle,'FontSize',12,'BackgroundColor',[1 1 1],...
    'Edgecolor',[0 0 0]);
 
    %colormap(hot);
    colormap([0 0 1;0 0 1]);


    [C,h]=contour(xr,-yr,PhiBr,[-10 -8 -3  -1 -.3 -.1 -.03],'LineWidth',1,'LabelSpacing',1000);
    text_handle=clabel(C,h);
    set(text_handle,'FontSize',12,'BackgroundColor',[1 1 1],...
    'Edgecolor',[0 0 0]);
    
    %colormap(hot);
    colormap([0 0 1;0 0 1]);
    
    axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
    axis equal
    plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
    ylabel('y','FontSize',22);
    xlabel('z=r cos(\theta)','FontSize',22);

else
    %figure;hold all
    C=contour(xr,yr,PhiTr,[-10 -8 -3  -1 -.3 -.1 -.03]);
    Clength=size(C,2);
    temp1=1;
    while(true)
        temp2=C(2,temp1)+temp1;
        if(temp2==Clength)
            plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k--','LineWidth',1)
            break;
        end
        plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k--','LineWidth',1)
        temp1=temp2+1;
    end

    C=contour(xr,-yr,PhiBr,[-10 -8 -3  -1 -.3 -.1 -.03]);
    Clength=size(C,2);
    temp1=1;
    while(true)
        temp2=C(2,temp1)+temp1;
        if(temp2==Clength)
            plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k--','LineWidth',1)
            break;
        end
        plot(C(1,temp1+1:temp2),C(2,temp1+1:temp2),'k--','LineWidth',1)
        temp1=temp2+1;
    end
    axis equal
    axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
    ylabel('y','FontSize',22);
    xlabel('z','FontSize',22);
end

axx=-rcc(nrused)*0.85+[0 vd]*c_d*rcc(nrused)*.2;
axy=-rcc(nrused)*0.95+[0 vd]*sqrt(1-c_d^2)*rcc(nrused)*.2;
[arrx,arry]=dsxy2figxy(gca, axx, axy);
har=annotation('textarrow',arrx,arry);
%content=sprintf('v_d=%4.2f',vd)
content=sprintf('');
set(har,'String',content,'HeadStyle','vback3','FontSize',14,'HeadLength',3,'HeadWidth',6,'Color','b')
%har=text(axx(1),axy(1),'s');
%content=sprintf('v_d');
%set(har,'String',content,'FontSize',14)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Velocity arrows %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(0)

[r,the]=meshgrid(2:Dr*8*prec:rcc(nrused)+Dr/nrused, ...
    pi/round(pi/(prec*Dtheta1)):pi/round(pi/(prec*Dtheta2)):pi-pi/round(pi/(prec*Dtheta1)));

%[r,the]=meshgrid(2:Dr*8*prec:rcc(nrused)+Dr/nrused, ...
%    pi/2:pi/round(pi/(prec*Dtheta2)):pi-pi/round(pi/(prec*Dtheta1)));

xr=r.*cos(the);
yr=r.*sin(the);

vxT=griddata(XR,YR,(vrT'.*T-vtT'.*sqrt(1-T.^2))./(psumT'+0.001),xr,yr);
vyT=griddata(XR,YR,(vtT'.*T+vrT'.*sqrt(1-T.^2))./(psumT'+0.001),xr,yr);
vxB=griddata(XR,YR,(vrB'.*T-vtB'.*sqrt(1-T.^2))./(psumB'+0.001),xr,yr);
vyB=griddata(XR,YR,(vtB'.*T+vrB'.*sqrt(1-T.^2))./(psumB'+0.001),xr,yr);

figure;hold all
quiver(xr,yr,vxT,vyT,0,'k')
quiver(xr,-yr,vxB,-vyB,0,'k')

end


%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fluid streamlines %%
%%%%%%%%%%%%%%%%%%%%%%%%

if(0)

[xr,yr]=meshgrid(-0.97*rcc(nrused):.2:0.97*rcc(nrused),0.1:.2:0.97*rcc(nrused));

vxT=griddata(XR,YR,(vrT'.*T-vtT'.*sqrt(1-T.^2))./(psumT'+0.001),xr,yr);
vyT=griddata(XR,YR,(vtT'.*T+vrT'.*sqrt(1-T.^2))./(psumT'+0.001),xr,yr);

vxB=griddata(XR,YR,(vrB'.*T-vtB'.*sqrt(1-T.^2))./(psumB'+0.001),xr,yr);
vyB=griddata(XR,YR,(vtB'.*T+vrB'.*sqrt(1-T.^2))./(psumB'+0.001),xr,yr);


yy=-yr(end:-1:1,:);
xx=xr(end:-1:1,:); % This is not necessary
YY=[yy;yr];
XX=[xr;xr];
vy=-vyB(end:-1:1,:);
vx=vxB(end:-1:1,:);
VY=[vy;vyT];
VX=[vx;vxT];

% To avoid plotting stream lines in the probe
RR2=XX.^2+YY.^2;
VX(RR2<1)=NaN;
VY(RR2<1)=NaN;

%figure
r=0.95*rcc(nrused);
cthe=-0.95:0.1:0.95;
sx=r.*cos(acos(cthe)-acos(c_d)+pi/2);
sy=-r.*sin(acos(cthe)-acos(c_d)+pi/2);
%[sx,sy] = meshgrid(-4:0.5:4,-4);
h=streamline(stream2(XX,YY,VX,VY,sx,sy));
set(h,'Color','b','LineStyle','--')
hold all;

end

% This fill causes troubles. Need to do it at the end
fill(cosinus,sinus,[0.9 0.9 0.9]);

%plot(cosinus,sinus,'magenta')
%plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
box on

end
