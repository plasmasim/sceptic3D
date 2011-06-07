function [cang,fluxofangle] = Postproc(filename)

% This functions plots the angular ion flux distribution as a funtion of
% cang, cosine of the angle along the drift axis. If v_d//e_z, cang=tcc.
%If B=0, the this should output a plot independent of the angle between v_d and e_z.


% Read output file
short=true;
readoutput();

for j=1:nthused
    tcc(j)=1-2*double(j-1)/double(nthused-1);
end
for i=1:npsiused
    pcc(i)=0.+double(i-1)*2*pi/double(npsiused);
end

% angsize is the size of cang. Better if odd
angsize=min(npsiused,nthused);

for j=1:angsize
    cang(j)=1-2*double(j-1)/double(angsize-1);
end

fluxofangle=zeros(1,angsize);
sheathpot=zeros(1,angsize);
aweight=zeros(1,angsize);
aweightp=zeros(1,angsize);

% Backward interpolation theta/psi->cang
for j=1:nthused
    for k=1:npsiused
% Interpolate position on the drift axis
        cangle=sqrt(1-tcc(j)^2)*sin(pcc(k))*sqrt(1-c_d^2)+tcc(j)*c_d;
        
        for ial=2:angsize
            if(cang(ial)<=cangle)
                ial=ial-1;
                break
            end
        end
        if(ial<angsize)
            af=(cangle-cang(ial))/(cang(ial+1)-cang(ial));
        else
            af=1;ial=angsize-1;
        end
        
        fluxofangle(ial)=fluxofangle(ial)+nincell(j,k)*(1-af);
        fluxofangle(ial+1)=fluxofangle(ial+1)+nincell(j,k)*af;     
        aweight(ial)=aweight(ial)+(1-af);
        aweight(ial+1)=aweight(ial+1)+af;
        if(or(j==1,j==nthused))
            aweight(ial)=aweight(ial)-0.5*(1-af);
            aweight(ial+1)=aweight(ial+1)-0.5*af;
        end  
        
        sheathpot(ial)=sheathpot(ial)+phi(1,j,k)*(1-af);
        sheathpot(ial+1)=sheathpot(ial+1)+phi(1,j,k)*af;
        aweightp(ial)=aweightp(ial)+(1-af);
        aweightp(ial+1)=aweightp(ial+1)+af;

        
    end
end

fluxofangle=fluxofangle.*double(npsiused*(nthused-1))./(aweight+1e-5)/(4*pi*rhoinf*dt*double(nastep));
sheathpot=sheathpot./(aweightp+1e-5);

vt=sqrt(2*Ti);
flux0=vt/(2*sqrt(pi));
fluxofangle=fluxofangle/flux0;
%fluxtot=trapz(-cang,fluxofangle)
rb=rcc(nrused);
npart=sum(sum(sum(psum)));

fluxtot=sum(sum(nincell))/(4*pi*rhoinf*dt*double(nastep))/flux0;

figure(1);
%title('','FontSize',22)
plot(cang,fluxofangle,'LineWidth',0.5)
xlabel('cos \chi','FontSize',22)
ylabel('\Gamma_i/\Gamma_i^0','FontSize',22)

figure(2);
%title('','FontSize',22)
plot(cang,sheathpot,'LineWidth',0.5)
xlabel('cos \chi','FontSize',22)
ylabel('\phi_s','FontSize',22)


end
        