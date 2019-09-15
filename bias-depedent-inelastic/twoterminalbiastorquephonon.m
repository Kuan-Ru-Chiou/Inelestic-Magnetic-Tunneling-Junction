clc;
clear all;
%non---equibrium  2017/11/19   check Torque-parallel V -; T-perpendicular not V^2 ;
%%%check ok
%%voltage asymmetry of spin transfer        contact1 is Z up and
%down.   contact 2 is x up and down



%constant and input
hbar = 1.06e-34; %plank constant
q = 1.6e-19;  %free electron charge
m0 = 9.1e-31; %free electon mass
mf =  0.73*m0        %FM effective mass
mO =  0.32*m0        %Oxside effective mass
a = 1e-10;       %MgO unit cell length ( I don't know real value , This is come from some article. Please check)
zplus = 1i*1e-18;
kT = 0.025; %ev in room temp.¡@¡@  % KT = 0.00025  A0 =0.05 hw = 1e-5 exchange =0 for phonon bias  good
KTs = kT;
IE = q*q/(2*pi*hbar); % A/ev
tO = (hbar^2)/(2*mO*(a^2)*q); %Oxside hopping energy (eV)  (oxide a and Fe a using the same value because i dont know real value)
tf = (hbar^2)/(2*mf*(a^2)*q); %FM hopping energy (eV)
exchange =0;  %exchange splitting energy
Ef = 2.25;  %fermi-energy
Ec = 0;  %conduction band minimum
NL =2;  NOx =10;  NR =2;
Np = NL+NOx+NR;
Ub = Ef+0.93 ; 
% Ub = 0 ;
%oxide barrier ev
%pauli matrix
sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1];


%rotation matrix
theta = 0; 
%theta = pi/4;
R = [cos(0.5*theta) -sin(0.5*theta); sin(0.5*theta) cos(0.5*theta)];


%construct hamiltonian

alphaL = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange-0.5*sz*exchange;
alphaL = kron(diag([ones(1,NL) zeros(1,NOx+NR)]),alphaL);


alphaox = [2*tO 0;0 2*tO];
alphaox = kron(diag([zeros(1,NL) ones(1,NOx) zeros(1,NR)]),alphaox);


alphaR = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange-0.5*(sx*sin(theta)+sz*cos(theta))*exchange;
alphaR = kron(diag([zeros(1,NOx+NL) ones(1,NR)]),alphaR);


beta = [ones(1,NL)*(-tf) ones(1,NOx-1)*(-tO) ones(1,NR)*(-tf)];
beta = kron(diag(beta,1),eye(2));

%potential barrier
UBB = [zeros(1,NL) Ub*ones(1,NOx) zeros(1,NR)];
UB = kron(diag(UBB),eye(2));

H = zeros(2*Np,2*Np);
H = H+alphaL+alphaR+alphaox+beta+beta';

%Scattering related data
% hw = 0.0001*tf;
hw = 1e-2;

N =  1/(exp(hw/KTs)-1); %Bose function
  A0 =0;  %A0 =0 ; no interactin
  
  
Dab = N*A0*eye(2*Np);
Dem =(N+1)*A0*eye(2*Np);
  

%    Jz = 0; Jx =0; Jy = 0; %initialize spin current in Jzy Jxy Jyy (transport in y direction)
%    Dos =0;
%    
%    Jz1 = 0; Jx1 =0; Jy1 = 0;
%    Jz2 = 0; Jx2 =0; Jy2= 0;
%    Jz3 = 0; Jx3 =0; Jy3 = 0;
   
%  NE = 600;
  
%define Gn, Gn, G, A
% Gn = zeros(2*Np,2*Np,NE);
% Gp = zeros(2*Np,2*Np,NE);
% G = zeros(2*Np,2*Np,NE);
% A = zeros(2*Np,2*Np,NE);
% sig1 = zeros(2*Np,2*Np,NE);
% sig2 = zeros(2*Np,2*Np,NE);
% gam1 =  zeros(2*Np,2*Np,NE);
% gam2 =  zeros(2*Np,2*Np,NE);
% gam1in = zeros(2*Np,2*Np,NE);
% gam2in = zeros(2*Np,2*Np,NE);
% Ssi = zeros(2*Np,2*Np,NE);
% Ssr = zeros(2*Np,2*Np,NE);
% Ss = zeros(2*Np,2*Np,NE);
% Ssin = zeros(2*Np,2*Np,NE);
% Ssout = zeros(2*Np, 2*Np,NE);
% Ssinnew = zeros(2*Np,2*Np,NE);
% Ssoutnew = zeros(2*Np,2*Np,NE);
%    
% I1 = zeros(NE,1);
% I2 = zeros(NE,1);
% Is = zeros(NE,1);
% sI1 = zeros(Nv,1);
% sI2 = zeros(Nv,1);
% sIs = zeros(Nv,1);
    

%     E = linspace(Ef-0.5*tf*1e-3,Ef+0.5*tf*1e-3,200);
%     dE = E(2)-E(1);
%bias 
Nv = 10;
V = linspace(-1,1,Nv);  


for iv = 1:Nv
    
    
    NE = 600;
    E = linspace(min(Ef-0.5*V(iv)-15*kT*sign(V(iv)),Ef+0.5*V(iv)+15*kT*sign(V(iv))),max(Ef-0.5*V(iv)-15*kT*sign(V(iv)),Ef+0.5*V(iv)+15*kT*sign(V(iv))),NE);
%      E = linspace(min(Ef-0.5*V(iv),Ef+0.5*V(iv)),max(Ef-0.5*V(iv),Ef+0.5*V(iv)),NE);
    dE = E(2)-E(1);
    
    
    
  %define Gn, Gn, G, A
Gn = zeros(2*Np,2*Np,NE);
Gp = zeros(2*Np,2*Np,NE);
G = zeros(2*Np,2*Np,NE);
A = zeros(2*Np,2*Np,NE);
sig1 = zeros(2*Np,2*Np,NE);
sig2 = zeros(2*Np,2*Np,NE);
gam1 =  zeros(2*Np,2*Np,NE);
gam2 =  zeros(2*Np,2*Np,NE);
gam1in = zeros(2*Np,2*Np,NE);
gam2in = zeros(2*Np,2*Np,NE);
Ssi = zeros(2*Np,2*Np,NE);
Ssr = zeros(2*Np,2*Np,NE);
Ss = zeros(2*Np,2*Np,NE);
Ssin = zeros(2*Np,2*Np,NE);
Ssout = zeros(2*Np, 2*Np,NE);
Ssinnew = zeros(2*Np,2*Np,NE);
Ssoutnew = zeros(2*Np,2*Np,NE);
   
% I1 = zeros(NE,1);
% I2 = zeros(NE,1);
% Is = zeros(NE,1);
% sI1 = zeros(Nv,1);
% sI2 = zeros(Nv,1);
% sIs = zeros(Nv,1);
%     
 
%phonon spectrum
hwn = floor(hw/dE);
scllim = hwn+1;
sculim = NE - (hwn+1);
mu1 = Ef+0.5*V(iv);
mu2 = Ef-0.5*V(iv);

    U = [0.5*V(iv)*ones(1,NL) V(iv)*linspace(0.5,-0.5,NOx) -0.5*V(iv)*ones(1,NR)] + UBB;
    U = kron(diag(U),eye(2));
   
 %%%equlibrium part
   for k = 1:NE
       
     fL = 1/(1+exp((E(k)-mu1)/kT));
     fR = 1/(1+exp((E(k)-mu2)/kT));  
     sig1(1,1,k) = selfdatta1D(E(k),U(1,1),tf);
     sig1(2,2,k) = selfdatta1D(E(k),U(1,1)+exchange,tf);
     sig2(2*Np-1:2*Np,2*Np-1:2*Np,k) = R*[selfdatta1D(E(k),U(2*Np,2*Np),tf) 0;0 selfdatta1D(E(k),U(2*Np,2*Np)+exchange,tf)]*R';
     gam1(:,:,k) = 1i*(sig1(:,:,k) - sig1(:,:,k)');  gam2(:,:,k) = 1i*(sig2(:,:,k) - sig2(:,:,k)');
     gam1in(:,:,k) = fL*gam1(:,:,k);   gam2in(:,:,k) = fR*gam2(:,:,k);
     G(:,:,k) = inv((E(k))*eye(2*Np)-sig1(:,:,k)-sig2(:,:,k)-H-U);
     Gn(:,:,k) = G(:,:,k)*(gam1in(:,:,k) + gam2in(:,:,k))*G(:,:,k)';
     A(:,:,k) = 1i*(G(:,:,k)-G(:,:,k)');
     Gp(:,:,k) = A(:,:,k)-Gn(:,:,k);
   end
   
 %%% inelastic phonon loop
    err = 1000;
     
    while err > 1e-7
        err1 = 0;
        err2 = 0;
        for k = scllim:sculim
            Ssinnew = Dem.*diag(diag(Gn(:,:,k+hwn)))+Dab.*diag(diag(Gn(:,:,k-hwn)));
            Ssoutnew =Dab.*diag(diag(Gp(:,:,k+hwn)))+Dem.*diag(diag(Gp(:,:,k-hwn)));
            err1 = err1+sum(sum(abs(Ssout(:,:,k)-Ssoutnew)))/...
                sum(sum(abs(Ssout(:,:,k))+abs(Ssoutnew)));
            err2 = err2+sum(sum(abs(Ssin(:,:,k)-Ssinnew)))/...
                sum(sum(abs(Ssin(:,:,k)+abs(Ssinnew))));
            Ssout(:,:,k) =(Ssout(:,:,k)+Ssoutnew)/2;
            Ssin(:,:,k) = (Ssin(:,:,k)+Ssinnew)/2;
            Ssi(:,:,k) = (-1i/2)*(Ssin(:,:,k)+Ssout(:,:,k));
            G(:,:,k) = inv((E(k)+zplus)*eye(2*Np)-H-U-sig1(:,:,k)-sig2(:,:,k)-Ss(:,:,k));
            Gn(:,:,k) = G(:,:,k)*(gam1in(:,:,k)+gam2in(:,:,k)+Ssin(:,:,k))*G(:,:,k)';
            A(:,:,k) = 1i*(G(:,:,k)-G(:,:,k)');
            Gp(:,:,k) = A(:,:,k)-Gn(:,:,k);
        end
        %Hilbert tranform for real part of Ss
        for ii = 1:2
            for jj = 1:2
                Ssr(ii,jj,:) = -imag(hilbert(imag(Ssi(ii,jj,:))));
            end
        end
        Ss = Ssr+Ssi;
        %Ss = Ssi;
        err = err1 + err2;
        fprintf(1,'Error = %g\n',err);
    end
    
    
    
    %calculate I
    for k = 1:NE
      
    I1(k,iv) = real(trace(gam1in(:,:,k)*A(:,:,k)-gam1(:,:,k)*Gn(:,:,k)));
    I2(k,iv) = real(trace(gam2in(:,:,k)*A(:,:,k)-gam2(:,:,k)*Gn(:,:,k)));
    Is(k,iv) = real(trace(Ssin(:,:,k)*A(:,:,k)-1i*(Ss(:,:,k)-Ss(:,:,k)')*Gn(:,:,k)));
    end
    
    
    
%     HH = H+U;
%     
    
    
end


gllobalcurrent = sum((I1+I2+Is)*dE,1); %1 meaans column sum


fprintf(1,'global current conservation = %g\n',gllobalcurrent) 


figure
plot(V(1:Nv),sum(I1(:,1:Nv),1)*dE,'o-');
xlabel('bias (V)')
ylabel('current')

%     for k = 1:NE
%        
%      for j = 1:Np-1
%      Iin(2:Np,k)=(trace(HH(2*j-1:2*j,2*j+1:2*j+2)*Gn(2*j+1:2*j+2,2*j-1:2*j,k)-HH(2*j+1:2*j+2,2*j-1:2*j)*Gn(2*j-1:2*j,2*j+1:2*j+2,k)));
%      end
%        
%     end
II1 = sum(I1(:,1:Nv),1)*dE;



% Iin(1,:) = I1;
% Iin(Np+1,:) = I2;
% L = 1:Np;
% Lc = 0:Np;
  
 %plot currents

%plot currents
% figure(13)
% hold on;
% grid on;
% box on;
% plot(V,sI1,'b*-');
% plot(V,sI2,'ro-');
% plot(V,sIs,'gx-');
% legend('Contact 1','Contact 2','Scatterer');
% title('Current v Voltage');
% xlabel('Voltage (V)');
% ylabel('Current (A)');
%     shading interp
   
    %%%%using keldish Gless formula to check datta formula (checking ok)
    %       Sigless = 1i*(gamL*fL + gamR*fR);   %%  antihermitian
    %       Gless =  g*Sigless*g';              %%  antihermitian
    
    
    %density(j)=real(trace(gn)); %electron density
%     Tp(k)= real(trace(gamL*g*gamR*g'));  %parallel magnet
%     I = I + 1/(2*pi)*dE*Tp(k)*(fL-fR)  %paralaell total current
%     
    
  


% AAA = diff(sI1)./diff(V');
    
%     HH = H+U;
%     
    %spin current at interface point for calculation torque
    %calculate Np-NR-3 spot's spincurrent
    
    %%%%%%%%%check spin current conserve in n,n-1,n+1 points  spin-conservation
    %%%%%%%%%
    
%     Jz1 = Jz1-trace(0.5*(HH(2*(Np-1)-1:2*(Np-1),2*(Np-1)-3:2*(Np-1)-2)*sz +sz*HH(2*(Np-1)-1:2*(Np-1),2*(Np-1)-3:2*(Np-1)-2))*Gless(2*(Np-1)-3:2*(Np-1)-2,2*(Np-1)-1:2*(Np-1))-...
%         0.5*(HH(2*(Np-1)-3:2*(Np-1)-2,2*(Np-1)-1:2*(Np-1))*sz +sz*HH(2*(Np-1)-3:2*(Np-1)-2,2*(Np-1)-1:2*(Np-1))*Gless(2*(Np-1)-1:2*(Np-1),2*(Np-1)-3:2*(Np-1)-2)));
%     %Jz1 = Jz1 - (tf/(4*pi))*trace(sz*(Gless(2*(Np-2)-1:2*(Np-2),2*(Np-2)-3:2*(Np-2)-2)-Gless(2*(Np-2)-3:2*(Np-2)-2,2*(Np-2)-1:2*(Np-2))));
%     Jzy1(i,k) = Jz1;
%     
%     
%     Jx1 = Jx1-trace(0.5*(HH(2*(Np-1)-1:2*(Np-1),2*(Np-1)-3:2*(Np-1)-2)*sx +sx*HH(2*(Np-1)-1:2*(Np-1),2*(Np-1)-3:2*(Np-1)-2))*Gless(2*(Np-1)-3:2*(Np-1)-2,2*(Np-1)-1:2*(Np-1))-...
%         0.5*(HH(2*(Np-1)-3:2*(Np-1)-2,2*(Np-1)-1:2*(Np-1))*sx +sx*HH(2*(Np-1)-3:2*(Np-1)-2,2*(Np-1)-1:2*(Np-1))*Gless(2*(Np-1)-1:2*(Np-1),2*(Np-1)-3:2*(Np-1)-2)));
%     %Jx1 = Jx1 - (tf/(4*pi))*trace(sx*(Gless(2*(Np-2)-1:2*(Np-2),2*(Np-2)-3:2*(Np-2)-2)-Gless(2*(Np-2)-3:2*(Np-2)-2,2*(Np-2)-1:2*(Np-2))));
%     Jxy1(i,k) = Jx1;
%     
%     
%     
%     Jy1 = Jy1-trace(0.5*(HH(2*(Np-1)-1:2*(Np-1),2*(Np-1)-3:2*(Np-1)-2)*sy +sy*HH(2*(Np-1)-1:2*(Np-1),2*(Np-1)-3:2*(Np-1)-2))*Gless(2*(Np-1)-3:2*(Np-1)-2,2*(Np-1)-1:2*(Np-1))-...
%         0.5*(HH(2*(Np-1)-3:2*(Np-1)-2,2*(Np-1)-1:2*(Np-1))*sy +sy*HH(2*(Np-1)-3:2*(Np-1)-2,2*(Np-1)-1:2*(Np-1))*Gless(2*(Np-1)-1:2*(Np-1),2*(Np-1)-3:2*(Np-1)-2)));
%     %Jy1 = Jy1 - (tf/(4*pi))*trace(sy*(Gless(2*(Np-2)-1:2*(Np-2),2*(Np-2)-3:2*(Np-2)-2)-Gless(2*(Np-2)-3:2*(Np-2)-2,2*(Np-2)-1:2*(Np-2))));
%     Jyy1(i,k) = Jy1;
%     
    


% Ipp(i) = I;
%       


%%every bias spin current at the right oxide FM interface 
% 
%  Jzy= sum(Jzy1,2)*dE;  
%  Jxy = sum(Jxy1,2)*dE;
%  Jyy = sum(Jyy1,2)*dE;
%    
% 
% Jp =sqrt(Jzy.^2 + Jxy.^2); % in plane current
% JO = Jyy;   %out plane current


% figure(1)
% plot(V,JO);
% xlabel('bias')
% ylabel('outplane torque')
% 
% figure(2)
% plot(V,Jp)
%     xlabel('bias')
% ylabel('inplane torque')



%check charge conservation
% for i=1:Np-1
%    J00(i) = (-1)*trace(H(2*i+1:2*i+2,2*i-1:2*i)*Gless(2*i-1:2*i,2*i+1:2*i+2) - H(2*i-1:2*i,2*i+1:2*i+2)*Gless(2*i+1:2*i+2,2*i-1:2*i));
% end
% 

