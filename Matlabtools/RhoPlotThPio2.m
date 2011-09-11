function [] = RhoPlotTh0(filename)

%function [] = RhoPlotTh0(filename)
%
% Color-plot of ion charge-density contours on the plane {0,e_x,e_y}. Also
% shows the fluid strealines.
%
% In order to have white color when the density
% is unperturbed (n=1), I call the function function
% One2ZeroColor to change the colormap.


short=false;readforce=false;
% There is a problem with opengl on the cmod workstations. Matlab just
% crashes if it calls it
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

ttopm=floor(nthused/2);

Rho=double(nthused-1)/(2)*(rho(:,ttopm,:)*(0-tcc(ttopm+1))...
    +rho(:,ttopm+1,:)*(tcc(ttopm)-0));

vr=double(nthused-1)/(2)*(vrsum(:,ttopm,:)*(0-tcc(ttopm+1))...
    +vrsum(:,ttopm+1,:)*(tcc(ttopm)-0));

vp=double(nthused-1)/(2)*(vpsum(:,ttopm,:)*(0-tcc(ttopm+1))...
    +vpsum(:,ttopm+1,:)*(tcc(ttopm)-0));

ps=double(nthused-1)/(2)*(psum(:,ttopm,:)*(0-tcc(ttopm+1))...
    +psum(:,ttopm+1,:)*(tcc(ttopm)-0));

Rho=reshape(Rho,nrused,npsiused);
vr=reshape(vr,nrused,npsiused);
vp=reshape(vp,nrused,npsiused);
ps=reshape(ps,nrused,npsiused);

Rho=[Rho(:,:) Rho(:,1)];
vr=[vr vr(:,1)];
vp=[vp vp(:,1)];
psum=[ps ps(:,1)];


prec=1;

Pcc=[pcc;0];

[R,T]=meshgrid(rcc,Pcc);

XR=R.*cos(T);
YR=R.*sin(T);

Dtheta1=pcc(2)-pcc(1);
Dr=(rcc(nrused)-1)./(nrused-1);

%%%%%%%%%%%%%%%%%%%%%
% Start contouring %%
%%%%%%%%%%%%%%%%%%%%%

[r,the]=meshgrid(linspace(1,rcc(nrused),40),linspace(0,2*pi,40));

xr=r.*cos(the);
yr=r.*sin(the);
Rhor=griddata(XR,YR,Rho',xr,yr);


sinus=sin([0:pi/32:2*pi]);
cosinus=cos([0:pi/32:2*pi]);

figure;hold all
plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,Rhor)
hold all
%colormap(hot);
shading interp;
cmin=0.1*round(10*min(min(min(Rhor))));
cmax=0.1*round(10*max(max(max(Rhor))));
%colormap(hot);
One2ZeroColor(cmin,cmax,1);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')

axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('x','FontSize',22);

% Now, retake the dat with more precision for the contours and velocity plots
[r,the]=meshgrid(linspace(1,rcc(nrused),nrused),linspace(0,2*pi,npsiused));

xr=r.*cos(the);
yr=r.*sin(the);
Rhor=griddata(XR,YR,Rho',xr,yr);


axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('x','FontSize',22);

%
% Do the contours by hand, to have them all back
%

if(0)
    figure;hold all
    [C,h]=contour(xr,yr,Rhor,[0.5 0.6 0.7 0.8 0.9 0.97 1 1.01 1.02 1.03 1.04 1.05 1.1 1.2],'LineWidth',1,'LabelSpacing',1000);
    text_handle=clabel(C,h);
    set(text_handle,'FontSize',13);

    
    axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
    axis equal
    plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
    ylabel('y','FontSize',22);
    xlabel('x','FontSize',22);

else
    C=contour(xr,yr,Rhor,[0.3 0.4 0.5 0.6 0.8 0.9 0.95 1.1]);
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

%%%%%%%%%%%%%%%%%%%%%%%
% Fluid streamlines
%%%%%%%%%%%%%%%%%%%%%%%

if(1)
    
[xr,yr]=meshgrid(-0.99*rcc(nrused):.2:0.99*rcc(nrused),-0.99*rcc(nrused):.2:0.99*rcc(nrused));


vx=griddata(XR,YR,(vr'.*cos(T)-vp'.*sin(T))./(psum'+0.001),xr,yr);
vy=griddata(XR,YR,(vp'.*cos(T)+vr'.*sin(T))./(psum'+0.001),xr,yr);


% To avoid plotting stream lines in the probe
rr2=yr.^2+xr.^2;
vx(rr2<1)=NaN;
vy(rr2<1)=NaN;


r=0.95*rcc(nrused);
cthe=-0.95:0.1:0.95;
sx=r.*cthe;
sy=-r.*sin(acos(cthe));
%[sx,sy] = meshgrid(-4:0.5:4,-4);
h=streamline(stream2(xr,yr,vx,vy,sx,sy));
set(h,'Color','b','LineStyle','--')
hold all;
plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
box on

end
%%%%%%%%%%%%%%%%%%%%%%%


end
