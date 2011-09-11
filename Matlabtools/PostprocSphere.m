function [] = PostprocSphere(filename)

% Read output file

opengl neverselect

short=true;
readoutput();

for j=1:nthused
    tcc(j)=1-2*double(j-1)/double(nthused-1);
end
for i=1:npsiused
    pcc(i)=0.+double(i-1)*2*pi/double(npsiused);
end

% The following does not make much of a difference.
%if(tcc(1)==1)
%    tcc(1)=0.25*(3+tcc(2));
%end
%if(tcc(nthused)==-1)
%    tcc(nthused)=0.25*(-3+tcc(nthused-1));
%end


Csound=sqrt(1+Ti);

fluxofangle=nincell*double(npsiused)*double(nthused-1)/(4*pi*rhoinf*dt*double(nastep))/Csound;
fluxofangle(1,:)=fluxofangle(1,:)*2;
fluxofangle(end,:)=fluxofangle(end,:)*2;

fluxofangle(:,end+1)=fluxofangle(:,1);
pcc(npsiused+1)=2*pi;

% 1st graph

[T,P]=meshgrid(tcc,pcc);
figure;
%pcolor(T',P',fluxofangle./sqrt(1-(cos(P').*sin(acos(T'))).^2))
pcolor(T',P',fluxofangle)
hold all
xlabel('cos\theta','FontSize',22)
ylabel('\psi','FontSize',22)
title('\Gamma_i/N_{\infty}c_{sI}','FontSize',22)
shading interp
colorbar;

% 2d graph
pcc(npsiused+1)=pcc(1);


for j=1:nthused
    for k=1:npsiused+1
        z(j,k)=tcc(j);
        x(j,k)=sqrt(1-tcc(j).^2).*cos(pcc(k));
        y(j,k)=sqrt(1-tcc(j).^2).*sin(pcc(k));
    end
end

figure
%surf(x,y,z,fluxofangle./sqrt(1-(cos(P').*sin(acos(T'))).^2));
surf(x,y,z,fluxofangle);
hold all
colorbar;


axis equal
shading interp
xlabel('x','FontSize',22)
ylabel('y','FontSize',22)
zlabel('z','FontSize',22)



end
        