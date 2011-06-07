function [] = PostprocForceB(folder,cart,Tiunits,outer)

% Plots the ion-drag force components as a function of ion magnetization.
% All parameters must be the same except for the magnetic field
%     folder:  folder where the list of output files to consider is placed
%              ("titles")
%     cart:    Boolean. If "true", plot Fx,Fy,Fz. If "false", plot Fx,FB,FBperp,
%              Fv.
%     Tiunits: Boolean. If "true", normalized the forces in ion thermal units
%     outer:   Boolean. It "true", overimpose the forces from the outer boundary.


    short=false;readforce=true;
    out=strcat(folder,'/titles')
    fid=fopen(out,'r');
    TS=textscan(fid,'%s%s%s');
    ts=TS{2};
    dim=size(ts);
    dim=dim(1);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Cycle through the output files %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for kk=1:dim
        f=ts{kk};
        filename=strcat(folder,'/',f);

        readoutput();
        
        Cang(kk)=c_d;
        B(kk)=Bz;
        
        Charge1(kk)=charge1;
        Charge2(kk)=charge2;
        
        XFconv(kk)=-dbl^2*Charge1(kk)*vd*Bz*(cB*sqrt(1-c_d^2)-sqrt(1-cB^2)*c_d);
        
        XFfield(kk)=dbl^2*Xffield1;
        XFelec(kk)=Xfelec1;
        XFion(kk)=Xfion1;
        XFlor(kk)=Xflor1;
        XFtot(kk)=Xftot1;
        
        XFfield2(kk)=dbl^2*Xffield2;
        XFelec2(kk)=Xfelec2;
        XFion2(kk)=Xfion2;
        XFlor2(kk)=Xflor2;
        XFtot2(kk)=Xftot2;
        
        YFfield(kk)=dbl^2*Yffield1;
        YFelec(kk)=Yfelec1;
        YFion(kk)=Yfion1;
        YFlor(kk)=Yflor1;
        YFtot(kk)=Yftot1;
        
        YFfield2(kk)=dbl^2*Yffield2;
        YFelec2(kk)=Yfelec2;
        YFion2(kk)=Yfion2;
        YFlor2(kk)=Yflor2;
        YFtot2(kk)=Yftot2;
        
        ZFfield(kk)=dbl^2*Zffield1;
        ZFelec(kk)=Zfelec1;
        ZFion(kk)=Zfion1;
        ZFlor(kk)=Zflor1;
        ZFtot(kk)=Zftot1;
        
        ZFfield2(kk)=dbl^2*Zffield2;
        ZFelec2(kk)=Zfelec2;
        ZFion2(kk)=Zfion2;
        ZFlor2(kk)=Zflor2;
        ZFtot2(kk)=Zftot2;
        
        %Ffield(kk)=ZFfield(kk)*c_d+YFfield(kk)*sqrt(1-c_d^2);
        %Felec(kk)=ZFelec(kk)*c_d+YFelec(kk)*sqrt(1-c_d^2);
        %Fion(kk)=ZFion(kk)*c_d+YFion(kk)*sqrt(1-c_d^2);
        %Ftot(kk)=ZFtot(kk)*c_d+YFtot(kk)*sqrt(1-c_d^2);
        
        %Ffield2(kk)=ZFfield2(kk)*c_d+YFfield2(kk)*sqrt(1-c_d^2);
        %Felec2(kk)=ZFelec2(kk)*c_d+YFelec2(kk)*sqrt(1-c_d^2);
        %Fion2(kk)=ZFion2(kk)*c_d+YFion2(kk)*sqrt(1-c_d^2);
        %Ftot2(kk)=ZFtot2(kk)*c_d+YFtot2(kk)*sqrt(1-c_d^2);
        
        
        % Calculate the internal Lorentz force: Only works if B is along
        % the z axis
        [FLx(kk),FLy(kk),Currtot] = InternalCurrent(filename,5);
         
    end
    
    
    ErrX=2*abs((XFtot2-XFtot)./(XFtot2+XFtot));
       
    
    % Sort for increasing Mag field
    XF=[B',XFfield',XFelec',XFion',XFlor',XFtot',YFfield',YFelec',YFion',YFlor',YFtot',ZFfield',...
        ZFelec',ZFion',ZFlor',ZFtot',XFfield2',XFelec2',XFion2',XFlor2',XFtot2',YFfield2',YFelec2',...
        YFion2',YFlor2',YFtot2',ZFfield2',ZFelec2',ZFion2',ZFlor2',ZFtot2',Charge1',Charge2',XFconv',...
        ErrX',FLx',FLy'];XF=sortrows(XF,1);
    XF=XF';
    B=XF(1,:);
    
    
    % Fix force units (B needs to be extracted before this)
    if(Tiunits)
        ytext='Forces / (R_p^2 N_{\infty} T_{i\infty})';
        XF=XF/Ti;
    else
        ytext='Forces / (R_p^2 N_{\infty} T_e)';
    end
    
    XFfield=XF(2,:);XFelec=XF(3,:);XFion=XF(4,:);XFlor=XF(5,:);XFtot=XF(6,:);
    YFfield=XF(7,:);YFelec=XF(8,:);YFion=XF(9,:);YFlor=XF(10,:);YFtot=XF(11,:);
    ZFfield=XF(12,:);ZFelec=XF(13,:);ZFion=XF(14,:);ZFlor=XF(15,:);ZFtot=XF(16,:);
    XFfield2=XF(17,:);XFelec2=XF(18,:);XFion2=XF(19,:);XFlor2=XF(20,:);XFtot2=XF(21,:);
    YFfield2=XF(22,:);YFelec2=XF(23,:);YFion2=XF(24,:);YFlor2=XF(25,:);YFtot2=XF(26,:);
    ZFfield2=XF(27,:);ZFelec2=XF(28,:);ZFion2=XF(29,:);ZFlor2=XF(30,:);ZFtot2=XF(31,:);
    ErrX=XF(35,:);
    Charge1=XF(32,:);
    Charge2=XF(33,:);
    XFconv=XF(34,:);
    FLx=XF(36,:);
    FLy=XF(37,:);
    
    % Go to beta units
    B=B*sqrt(2/(pi*Ti));
    B=B./(1+B);
    %figure
    %plot(B,ErrX)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%% Start plotting %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % capacitance
    %figure
    %Capac=Charge1/(4*pi*Vp)
    %plot(B,(Capac-1).^(-1))
    %return
    
    
     % Additional forces
    figure
    plot(B,FLx,'b^-',B,FLy,'b^--',B,XFconv,'ko-',B,XFtot,'rs-',B,YFtot,'rs--','LineWidth',1)
    legend('F_{J x}','F_{J y}','F_Q','Drag_x','Drag_y')
    xlabel('\beta_i/(1+\beta_i)','FontSize',22)
    ylabel(ytext,'FontSize',22)
    
    % Fx
    figure
    plot(B,XFfield,'b^-',B,XFelec,'k-',B,XFion,'og-',B,XFlor,'vc-',B,XFtot,'rs-',B,XFconv,'ko','LineWidth',1)
    hold all
    if(outer)
        plot(B,XFfield2,'b^--',B,XFelec2,'k--',B,XFion2,'og--',B,XFlor2,'vc--',B,XFtot2,'rs--','MarkerSize',8,'LineWidth',1)
    end
    title('Forces along the x axis','FontSize',22)
    legend('Ffield','Felec','Fion','FMag','Ftot','Fconv')
    xlabel('\beta_i/(1+\beta_i)','FontSize',22)
    ylabel(ytext,'FontSize',22)
    
    if(cart)
    
    % Fy    
    figure
    plot(B,YFfield,'b^-',B,YFelec,'k-',B,YFion,'og-',B,YFlor,'vc-',B,YFtot,'rs-','LineWidth',1)
    hold all
    if(outer)
        plot(B,YFfield2,'b^--',B,YFelec2,'k--',B,YFion2,'og--',B,YFlor2,'vc--',B,YFtot2,'rs--','MarkerSize',8,'LineWidth',1)
    end
    title('Forces along the y axis ','FontSize',22)
    legend('Ffield','Felec','Fion','FMag','Ftot')
    xlabel('\beta_i/(1+\beta_i)','FontSize',22)
    ylabel(ytext,'FontSize',22)
    
    
    % Fz
    figure
    plot(B,ZFfield,'b^-',B,ZFelec,'k-',B,ZFion,'og-',B,ZFlor,'vc-',B,ZFtot,'rs-','LineWidth',1)
    hold all
    if(outer)
        plot(B,ZFfield2,'b^--',B,ZFelec2,'k--',B,ZFion2,'og--',B,ZFlor2,'vc--',B,ZFtot2,'rs--','MarkerSize',8,'LineWidth',1)
    end
    title('Forces along the z axis','FontSize',22)
    legend('Ffield','Felec','Fion','FMag','Ftot')
    xlabel('\beta_i/(1+\beta_i)','FontSize',22)
    ylabel(ytext,'FontSize',22)
    
    else
    % Plot forces in the B, Bperp, v direction

        a=acos(cB);
        R=[cos(a) sin(a);-sin(a) cos(a)];
        for k=1:dim;
            BFfield(:,k)=R*[ZFfield(k);YFfield(k)];
            BFfield2(:,k)=R*[ZFfield2(k);YFfield2(k)];
            BFelec(:,k)=R*[ZFelec(k);YFelec(k)];
            BFelec2(:,k)=R*[ZFelec2(k);YFelec2(k)];
            BFion(:,k)=R*[ZFion(k);YFion(k)];
            BFion2(:,k)=R*[ZFion2(k);YFion2(k)];
            BFlor(:,k)=R*[ZFlor(k);YFlor(k)];
            BFlor2(:,k)=R*[ZFlor2(k);YFlor2(k)];
            BFtot(:,k)=R*[ZFtot(k);YFtot(k)];
            BFtot2(:,k)=R*[ZFtot2(k);YFtot2(k)];
        end
        BFfield1=BFfield(1,:);BFfield12=BFfield2(1,:);
        BFelec1=BFelec(1,:);BFelec12=BFelec2(1,:);
        BFion1=BFion(1,:);BFion12=BFion2(1,:);
        BFlor1=BFlor(1,:);BFlor12=BFlor2(1,:);
        BFtot1=BFtot(1,:);BFtot12=BFtot2(1,:);
        
        PFfield1=BFfield(2,:);PFfield12=BFfield2(2,:);
        PFelec1=BFelec(2,:);PFelec12=BFelec2(2,:);
        PFion1=BFion(2,:);PFion12=BFion2(2,:);
        PFlor1=BFlor(2,:);PFlor12=BFlor2(2,:);
        PFtot1=BFtot(2,:);PFtot12=BFtot2(2,:);
        
        
        % FB
        figure
        plot(B,BFfield1,'b^-',B,BFelec1,'k-',B,BFion1,'og-',B,BFlor1,'vc-',B,BFtot1,'rs-','LineWidth',1)
        hold all
        if(outer)
            plot(B,BFfield12,'b^--',B,BFelec12,'k--',B,BFion12,'og--',B,BFlor12,'vc--',B,BFtot12,'rs--','MarkerSize',8,'LineWidth',1)
        end
        title('Forces along the B axis','FontSize',22)
        legend('Ffield','Felec','Fion','FMag','Ftot')
        xlabel('\beta_i/(1+\beta_i)','FontSize',22)
        ylabel(ytext,'FontSize',22)
        
        
        % FBperp
        figure
        plot(B,PFfield1,'b^-',B,PFelec1,'k-',B,PFion1,'og-',B,PFlor1,'vc-',B,PFtot1,'rs-','LineWidth',1)
        hold all
        if(outer)
            plot(B,PFfield12,'b^--',B,PFelec12,'k--',B,PFion12,'og--',B,PFlor12,'vc--',B,PFtot12,'rs--','MarkerSize',8,'LineWidth',1)
        end
        title('Forces perp to the B axis','FontSize',22)
        legend('Ffield','Felec','Fion','FMag','Ftot')
        xlabel('\beta_i/(1+\beta_i)','FontSize',22)
        ylabel(ytext,'FontSize',22)
        
        a=acos(c_d);
        R=[cos(a) sin(a);-sin(a) cos(a)];
        for k=1:dim;
            VFfield(:,k)=R*[ZFfield(k);YFfield(k)];
            VFfield2(:,k)=R*[ZFfield2(k);YFfield2(k)];
            VFelec(:,k)=R*[ZFelec(k);YFelec(k)];
            VFelec2(:,k)=R*[ZFelec2(k);YFelec2(k)];
            VFion(:,k)=R*[ZFion(k);YFion(k)];
            VFion2(:,k)=R*[ZFion2(k);YFion2(k)];
            VFlor(:,k)=R*[ZFlor(k);YFlor(k)];
            VFlor2(:,k)=R*[ZFlor2(k);YFlor2(k)];
            VFtot(:,k)=R*[ZFtot(k);YFtot(k)];
            VFtot2(:,k)=R*[ZFtot2(k);YFtot2(k)];
        end
        VFfield1=VFfield(1,:);VFfield12=VFfield2(1,:);
        VFelec1=VFelec(1,:);VFelec12=VFelec2(1,:);
        VFion1=VFion(1,:);VFion12=VFion2(1,:);
        VFlor1=VFlor(1,:);VFlor12=VFlor2(1,:);
        VFtot1=VFtot(1,:);VFtot12=VFtot2(1,:);
        
        %Fv
        figure
        plot(B,VFfield1,'b^-',B,VFelec1,'k-',B,VFion1,'og-',B,VFlor1,'vc-',B,VFtot1,'rs-','LineWidth',1)
        hold all
        if(outer)
            plot(B,VFfield12,'b^--',B,VFelec12,'k--',B,VFion12,'og--',B,VFlor12,'vc--',B,VFtot12,'rs--','MarkerSize',8,'LineWidth',1)
        end
        title('Forces along the v axis','FontSize',22)
        legend('Ffield','Felec','Fion','FMag','Ftot')
        xlabel('\beta_i/(1+\beta_i)','FontSize',22)
        ylabel(ytext,'FontSize',22)
        
    end
    
end
  
    
    