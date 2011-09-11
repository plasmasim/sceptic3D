function [] = PostprocB(folder)
%   Script reading the list of SCEPTIC3D output files from the text file
%   'folder/titles', and plotting the total ion current to the sphere for
%   increasing Mag field, and normalized to the ion thermal current.
    

% We don't need all the infos
    short=true;

    out=strcat(folder,'/titles')

    fid=fopen(out,'r');

    TS=textscan(fid,'%s%s%s');
    ts=TS{2};

    

% Loop over the output files
    dim=size(ts,1);
    
    for kk=1:dim
        f=ts{kk};
        filename=strcat(folder,'/',f);

        % Read the current output file
        readoutput();
        
        Cang(kk)=c_d;
        L(kk)=dbl;
        B(kk)=Bz;
        vt=sqrt(2*Ti);
        
        % Ion thermal flux
        flux0=vt/(2*sqrt(pi));
        
        %Total ion current to the sphere, normalized to the ion thermal
        %current
        fluxtot(kk)=sum(sum(nincell))/(4*pi*rhoinf*dt*double(nastep))/flux0;
        
    end
    
    % Sort for increasing Mag field
    XF=[B',fluxtot'];XF=sortrows(XF,1);
    XF=XF';
    B=XF(1,:);
    fluxtot=XF(2,:);
       
    Beta=B/sqrt(Ti*pi/2);

    plot(Beta./(1+Beta),fluxtot,'ko-','LineWidth',1)       
    xlabel('\beta_i/(1+\beta_i)','FontSize',22);
    ylabel('I_i/I_i^0','FontSize',22);
    
end

    
    