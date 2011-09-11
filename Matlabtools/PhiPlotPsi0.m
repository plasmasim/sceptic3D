function [] = PhiPlotPsi0(filename)


%function [] = PhiPlotPsi0(filename)
%
% Color-plot of potential contours on the plane {0,e_x,e_z}. Also
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

ptopm=1;
pbotm=ceil((npsiused-1)/2)+1;

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
ylabel('x','FontSize',22);
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
    ylabel('x','FontSize',22);
    xlabel('z','FontSize',22);

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
    ylabel('x','FontSize',22);
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


box on


end