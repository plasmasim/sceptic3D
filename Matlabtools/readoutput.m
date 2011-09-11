%   Script able to read SCEPTIC3D outputs, from "filename"
%   I use a script rather than a
%   function because this is called by other diagnostic routines, and I did
%   not want to pass on every variables by hand.

% Open the file
    filename
	fid=fopen(filename,'r');
    
% No collisions for now in SCEPTIC3D
    ncoll=1;

    % Read the header  

    if(ncoll)
        TS=textscan(fid,'%f',13,'headerlines',1);
        TS=TS{1};
    
        dt=TS(1);vd=TS(2);c_d=TS(3);Ti=TS(4);steps=TS(5);rhoinf=TS(6);
        phiinf=TS(7);fave=TS(8);dbl=TS(9);Vp=TS(10);
        Bz=TS(11);cB=TS(12);nrused=TS(13);%coll=TS(13);nrused=TS(14);

    else
        TS=textscan(fid,'%f',15,'headerlines',1);
        TS=TS{1};
    
        dt=TS(1);vd=TS(2);c_d=TS(3),Ti=TS(4);steps=TS(5);rhoinf=TS(6);
        phiinf=TS(7);fave=TS(8);dbl=TS(9);Vp=TS(10);
        Bz=TS(11);cB=TS(12);nrused=TS(13);coll=TS(14);nrused=TS(15)
    end

% Read mesh diag informations
    TS=textscan(fid,'%f%f%f',nrused,'headerlines',1);
    rcc=TS{1};diagphi=TS{2};diagrho=TS{3};
     
% Read flux to probe at each step
    TS=textscan(fid,'%d',steps,'headerlines',3);
    fluxprobe=TS{1};
    
    
% Read particle flux data
    TS=textscan(fid,'%d',3,'headerlines',2);
    TS=TS{1};
    nthused=TS(1);npsiused=TS(2);

    TS=textscan(fid,'%d','headerlines',1);    
    
% Read averaged particle flux data (flux, average velocity and average
% square velocity)
    TS =textscan(fid,'%d',1,'headerlines',1);
    nastep=TS{1};    

    for k=1:npsiused
        TS=textscan(fid,'%f',nthused);
        TS=TS{1};
        if (k==1);
            nincell=TS;
        else
            nincell=horzcat(nincell,TS);
        end
    end

    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        TS=textscan(fid,'%f',nthused);
        TS=TS{1};
        if (k==1);
            vrincell=TS;
        else
            vrincell=horzcat(vrincell,TS);
        end
    end
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        TS=textscan(fid,'%f',nthused);
        TS=TS{1};
        if (k==1);
            vr2incell=TS;
        else
            vr2incell=horzcat(vr2incell,TS);
        end
    end


% Read Mesh potential
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    
    for k=1:nrused
        for l=1:npsiused
            TS=textscan(fid,'%f',nthused);
            TS=TS{1};
            if (l==1);
                phiA=TS;
            else
                phiA=horzcat(phiA,TS);
            end
        end
        if(k==1)
            phi=phiA;
        else
            phi=cat(3,phi,phiA);
        end
    end
    phi=permute(phi,[3 1 2]);
    clear phiA;
    
% Read Mesh density (normalised to
% rhoinf)
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:nrused
        for l=1:npsiused
            TS=textscan(fid,'%f',nthused);
            TS=TS{1};
            if (l==1);
                rhoA=TS;
            else
                rhoA=horzcat(rhoA,TS);
            end
        end
        if(k==1)
            rho=rhoA;
        else
            rho=cat(3,rho,rhoA);
        end
    end
    rho=permute(rho,[3 1 2]);
    clear rhoA; 

% Read inverse volume of cells
    TS=textscan(fid,'%f','headerlines',2);
    volinv=TS{1}; 
     
    TS=textscan(fid,'%f',8,'headerlines',2);
    TS=TS{1};
    dt=TS(1);vd=TS(2);Ti=TS(3);steps=TS(4);rmax=TS(5);rhoinf=TS(6);
    dbl=TS(7);Vp=TS(8);

% Read psum
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',2);
    for k=1:npsiused
        for l=1:nthused
            TS=textscan(fid,'%f',nrused);
            TS=TS{1};
            if (l==1);
                psumA=TS;
            else
                psumA=horzcat(psumA,TS);
            end
        end
        if(k==1)
            psum=psumA;
        else
            psum=cat(3,psum,psumA);
        end
    end
    clear psumA;
        
 
% short enables to skip the following if not needed, in order to save time
   if(short==false) 
    
        
% Read vrsum
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        for l=1:nthused
            TS=textscan(fid,'%f',nrused);
            TS=TS{1};
            if (l==1);
                vrsumA=TS;
            else
                vrsumA=horzcat(vrsumA,TS);
            end
        end
        if(k==1)
            vrsum=vrsumA;
        else
            vrsum=cat(3,vrsum,vrsumA);
        end
    end
    clear vrsumA;

    
% Read vtsum
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        for l=1:nthused
            TS=textscan(fid,'%f',nrused);
            TS=TS{1};
            if (l==1);
                vtsumA=TS;
            else
                vtsumA=horzcat(vtsumA,TS);
            end
        end
        if(k==1)
            vtsum=vtsumA;
        else
            vtsum=cat(3,vtsum,vtsumA);
        end
    end
    clear vtsumA;
   
% Read vpsum
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        for l=1:nthused
            TS=textscan(fid,'%f',nrused);
            TS=TS{1};
            if (l==1);
                vpsumA=TS;
            else
                vpsumA=horzcat(vpsumA,TS);
            end
        end
        if(k==1)
            vpsum=vpsumA;
        else
            vpsum=cat(3,vpsum,vpsumA);
        end
    end
    clear vpsumA;

 
% Read vr2sum
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        for l=1:nthused
            TS=textscan(fid,'%f',nrused);
            TS=TS{1};
            if (l==1);
                vr2sumA=TS;
            else
                vr2sumA=horzcat(vr2sumA,TS);
            end
        end
        if(k==1)
            vr2sum=vr2sumA;
        else
            vr2sum=cat(3,vr2sum,vr2sumA);
        end
    end
    clear vr2sumA;

% Read vt2sum
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        for l=1:nthused
            TS=textscan(fid,'%f',nrused);
            TS=TS{1};
            if (l==1);
                vt2sumA=TS;
            else
                vt2sumA=horzcat(vt2sumA,TS);
            end
        end
        if(k==1)
            vt2sum=vt2sumA;
        else
            vt2sum=cat(3,vt2sum,vt2sumA);
        end
    end    
    clear vt2sumA;
    
% Read vp2sum
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        for l=1:nthused
            TS=textscan(fid,'%f',nrused);
            TS=TS{1};
            if (l==1);
                vp2sumA=TS;
            else
                vp2sumA=horzcat(vp2sumA,TS);
            end
        end
        if(k==1)
            vp2sum=vp2sumA;
        else
            vp2sum=cat(3,vp2sum,vp2sumA);
        end
    end
    clear vp2sumA;
    
% Read vrtsum
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        for l=1:nthused
            TS=textscan(fid,'%f',nrused);
            TS=TS{1};
            if (l==1);
                vrtsumA=TS;
            else
                vrtsumA=horzcat(vrtsumA,TS);
            end
        end
        if(k==1)
            vrtsum=vrtsumA;
        else
            vrtsum=cat(3,vrtsum,vrtsumA);
        end
    end
    clear vrtsumA;
    
% Read vrpsum
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        for l=1:nthused
            TS=textscan(fid,'%f',nrused);
            TS=TS{1};
            if (l==1);
                vrpsumA=TS;
            else
                vrpsumA=horzcat(vrpsumA,TS);
            end
        end
        if(k==1)
            vrpsum=vrpsumA;
        else
            vrpsum=cat(3,vrpsum,vrpsumA);
        end
    end
    clear vrpsumA;
    
% Read vtpsum
    TS=textscan(fid,'%s',1,'delimiter','\n','headerlines',1);
    for k=1:npsiused
        for l=1:nthused
            TS=textscan(fid,'%f',nrused);
            TS=TS{1};
            if (l==1);
                vtpsumA=TS;
            else
                vtpsumA=horzcat(vtpsumA,TS);
            end
        end
        if(k==1)
            vtpsum=vtpsumA;
        else
            vtpsum=cat(3,vtpsum,vtpsumA);
        end
    end
    clear vtpsumA;

    
% Read crap
    TS=textscan(fid,'%f',nrused,'headerlines',2);
    TS=textscan(fid,'%f','headerlines',2);
    
% Read tcc
    TS=textscan(fid,'%[^t]',1);
    TS=textscan(fid,'%f',nthused,'headerlines',1);
    tcc=TS{1};
% Read pcc
    TS=textscan(fid,'%f',npsiused,'headerlines',2);
    pcc=TS{1};
    
        if(readforce)

% Read various forces ...

        TS=textscan(fid,'%f','headerlines',2);
        TS=TS{1};
        charge1=TS(1);charge2=TS(7);Zffield1=TS(2);Zffield2=TS(8);
        Zfelec1=TS(3);Zfelec2=TS(9);Zfion1=TS(4);Zfion2=TS(10);Zflor1=TS(5);Zflor2=TS(11);
        Zftot1=TS(6);Zftot2=TS(12);
    
        TS=textscan(fid,'%f','headerlines',1);
        TS=TS{1};
        Xffield1=TS(1);Xffield2=TS(6);
        Xfelec1=TS(2);Xfelec2=TS(7);Xfion1=TS(3);Xfion2=TS(8);Xflor1=TS(4);Xflor2=TS(9);
        Xftot1=TS(5);Xftot2=TS(10);
    
        TS=textscan(fid,'%f','headerlines',1);
        TS=TS{1};
        Yffield1=TS(1);Yffield2=TS(6);
        Yfelec1=TS(2);Yfelec2=TS(7);Yfion1=TS(3);Yfion2=TS(8);Yflor1=TS(4);Yflor2=TS(9);
        Yftot1=TS(5);Yftot2=TS(10);
        
% Energy flux to the probe
        TS=textscan(fid,'%f','headerlines',2);
        TS=TS{1};
        encoll=TS(1);

        end
   
   end
   
% close file
    fclose(fid);
        
 %   end