%Sina Mohammadi - Flow over a flat plate
clc
clear all
close all
%Parameters
%Flow is air, temperature is 25 Centigrade
U=0.75;
L=0.02;
Density=1.184;
D_Viscos=1.849*10^-5;
K_Viscos=D_Viscos/Density;
Reynolds=(U*L)/K_Viscos;
m=950;
n=950;
L_prime=(50*L)/(Reynolds^0.5);
deltax=(L/(m-1));
deltay=(L_prime/(n-1));
X=linspace(0,L,m);
Y=linspace(0,L_prime,n);
Reynoldsnolds=(U.*X)./K_Viscos;
u=zeros(n,m);
v=zeros(n,m);
for j=1:n-1
    unew(j,1)=U;
    u(j,1)=U;
    vnew(j,1)=0;
end
for i=1:m
    unew(1,i)=0;
    vnew(1,i)=0;
    u(n,i)=U;
end
for i=1:m
    unew(n,i)=U;
    u(1,i)=0;
    vnew(n,i)=0;
end
 
e=zeros(n,m);
u0=zeros(n,m);
u0(:,:)=0.01;
 
%% u
for i=1:m-1
    
   e(:,i+1)=1;
   while max(e(:,i+1))>0.000001
       
      for j=2:n-1
          
        a=(K_Viscos*deltax)/(unew(j,i)*(deltay)^2);
        b=(vnew(j,i)*deltax)/(2*unew(j,i)*deltay);
        %unew(j,i+1)=(unew(j,i)-b*(unew(j+1,i)-unew(j-1,i))+a*(unew(j-1,i+1)+unew(j+1,i+1)))/(1+2*a);
       unew(j,i+1)=(unew(j,i)-b*(unew(j+1,i)-unew(j-1,i))+a*(unew(j-1,i+1)+unew(j+1,i+1)))/(1+2*a);
      end
        
      e(:,i+1)=abs(unew(:,i+1)-u0(:,i+1));
      u0(:,i+1)=unew(:,i+1);
 
       
   end
    
   for j=2:n-1
        a=(K_Viscos*deltax)/(unew(j,i)*(deltay)^2);
      
        b=(vnew(j,i)*deltax)/(2*unew(j,i)*deltay);
        %unew(j,i+1)=(unew(j,i)-b*(unew(j+1,i)-unew(j-1,i))+a*(unew(j-1,i+1)+unew(j+1,i+1)))/(1+2*a);
       vnew(j,i+1)=vnew(j-1,i+1)-(deltay/(2*deltax))*(unew(j,i+1)-unew(j,i)+unew(j-1,i+1)-unew(j-1,i));
    
    end
  
end
 
 
%% shear
for i=1:m
    Shear(i,1)=D_Viscos*((unew(2,i)-unew(1,i))/deltay);
end
X1=X';
Y1=Y';
Shear_Blasius=(0.412*Density*U^2)./(Reynoldsnolds.^0.5);
 
Shear_Blasius1=Shear_Blasius';
 
%% delta
Delta_Blasius=(5*X)./(Reynoldsnolds.^0.5);
 
Delta_Blasius1=Delta_Blasius';
for i=1:m
    for j=1:n
        if unew(j,i)>=0.99*U
            delta(i,1)=deltay*(j-1);
            break
        end
    end
end
 
%Figures
%Contour of velocity
figure;
contourf(unew);
colormap(jet);
colorbar
%ylim([0 150])
title('Contour of velocity over a flat plate')
%shear comparison
figure;
plot(Shear_Blasius,'b--','linewidth',1.5);
hold on
plot(Shear,'r','linewidth',1);
title('Comparison of shear stress results from numerical solution and software solution','FontSize',7.5)
legend({'Numerical Solution','Software Solution'},'Location','northeast')
%boundary layer comparison
figure;
plot(Delta_Blasius,'b:','linewidth',1.5);
hold on
plot(delta,'r','linewidth',1);
title('Comparison of the thickness of the boundary layer fluid flow from numerical solution and software solution','FontSize',7.5)
legend({'Numerical Solution','Software Solution'},'Location','southeast')
