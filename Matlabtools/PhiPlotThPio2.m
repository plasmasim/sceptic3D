function [] = PhiPlotTh0(filename)


%function [] = PhiPlotTh0(filename)
%
% Color-plot of potential contours on the plane {0,e_x,e_y}. Also
% shows the fluid strealines.
%
% In order to have white color when the potential
% is unperturbed (phi=0), I call the function function
% One2ZeroColor to change the colormap.

short=false;readforce=false;
opengl neverselect


readoutput();

% Actual center of the first and last theta cells
%if(tcc(1)==1)
%    tcc(1)=0.25*(3+tcc(2));
%end
%if(tcc(nthused)==-1)
%    tcc(nthused)=0.25*(-3+tcc(nthused-1));
%end

%for i=1:npsiused
%    pcc(i)=0.+double(i-1)*2*pi/double(npsiused);
%end

ttopm=floor(nthused/2);

Phi=double(nthused)/(2)*(phi(:,ttopm,:)*(0-tcc(ttopm+1))...
    +phi(:,ttopm+1,:)*(tcc(ttopm)-0));

vr=double(nthused-1)/(2)*(vrsum(:,ttopm,:)*(0-tcc(ttopm+1))...
    +vrsum(:,ttopm+1,:)*(tcc(ttopm)-0));

vp=double(nthused-1)/(2)*(vpsum(:,ttopm,:)*(0-tcc(ttopm+1))...
    +vpsum(:,ttopm+1,:)*(tcc(ttopm)-0));

ps=double(nthused-1)/(2)*(psum(:,ttopm,:)*(0-tcc(ttopm+1))...
    +psum(:,ttopm+1,:)*(tcc(ttopm)-0));

%vt=double(nthused)/(2)*(vtsum(:,ttopm,:)*(0-tcc(ttopm+1))...
%    +vtsum(:,ttopm+1,:)*(tcc(ttopm)-0));

for k=1:prod(size(pcc))
    vx(:,k)=-vp(:,k)*sin(pcc(k))+vr(:,k)*cos(pcc(k));
    vy(:,k)=vp(:,k)*cos(pcc(k))+vr(:,k)*sin(pcc(k));
end


Phi=reshape(Phi,nrused,npsiused);
vr=reshape(vr,nrused,npsiused);
vp=reshape(vp,nrused,npsiused);
ps=reshape(ps,nrused,npsiused);

Phi=[Phi(:,:) Phi(:,1)];
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

[r,the]=meshgrid(1:Dr*prec:rcc(nrused)+Dr/nrused,[0:pi/round(pi/(prec*Dtheta1)):2*pi]);
xr=r.*cos(the);
yr=r.*sin(the);

Phir=griddata(XR,YR,Phi',xr,yr);

sinus=sin([0:pi/32:2*pi]);
cosinus=cos([0:pi/32:2*pi]);
figure;hold all
plot(cosinus,sinus,'k')
plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')

pcolor(xr,yr,Phir)
hold all
%colormap(hot);
shading interp;
cmin=0.1*round(10*min(min(min(Phir))));
cmax=max(0,0.1*round(10*max(max(max(Phir)))));
%colormap(hot);
One2ZeroColor(cmin,cmax,0);
caxis([cmin-10^-5,cmax+10^-5]);
colorbar('EastOutside')

axis equal
axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
ylabel('y','FontSize',22);
xlabel('x','FontSize',22);


if(0)
    figure;hold all;
    
    cosinus=cos([0:pi/32:2*pi]);
    plot(cosinus,sinus,'k')
    plot(rcc(nrused)*cosinus,rcc(nrused)*sinus,'k')
    fill(cosinus,sinus,[.85 .85 .85]);
    
    [C,h]=contour(xr,yr,Phir,[-12 -10 -8 -6 -3  -1 -.3],'LineWidth',1,'LabelSpacing',1000);
    text_handle=clabel(C,h);
    set(text_handle,'FontSize',12,'BackgroundColor',[1 1 1],...
    'Edgecolor',[0 0 0]);
    cmin=0.1*round(10*min(min(min(Phir))));
    cmax=max(0,0.1*round(10*max(max(max(Phir)))));
    %colormap(jet);
    colormap([0 0 1;0 0 1]);


    axis equal
    axis([-rcc(nrused) rcc(nrused) -rcc(nrused) rcc(nrused)]);
    ylabel('y','FontSize',22);
    xlabel('x','FontSize',22);
else
    C=contour(xr,yr,Phir,[-10 -8 -3 -1 -.3 -.1 -.03]);
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

if(0)

if(0)
    prec=2;
    [r,the]=meshgrid(1:Dr*3*prec:rcc(nrused)+Dr/nrused,[0:pi/round(pi/(prec*Dtheta1)):2*pi]);
    xr=r.*cos(the);
    yr=r.*sin(the);
    vxT=griddata(XR,YR,vx'./(psum'+0.001),xr,yr);
    vyT=griddata(XR,YR,vy'./(psum'+0.001),xr,yr);
    figure;hold all
    quiver(xr,yr,vxT,vyT,1,'k')
else
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

end

end