clc;
clear all;
close all;
;
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
kT = 0.025; %ev in room temp.
KTs = kT;
IE = q*q/(2*pi*hbar); % A/ev
tO = (hbar^2)/(2*mO*(a^2)*q); %Oxside hopping energy (eV)  (oxide a and Fe a using the same value because i dont know real value)
tf = (hbar^2)/(2*mf*(a^2)*q); %FM hopping energy (eV)
exchange =2.15;  %exchange splitting energy
Ef = 2.25;  %fermi-energy

Ec = 0;  %conduction band minimum
NL =2;  NOx =10;  NR =2;
Np = NL+NOx+NR;
Ub =Ef+0.93;     %oxide barrier ev

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


%bias
V = 0.05;


%  %energy grid
% E = linspace(Ef-0.5*0.1,Ef+0.5*0.1,500);
E = linspace(Ef-0.5*V-15*kT,Ef+0.5*V+15*kT,600);

dE = E(2)-E(1);
NE = length(E);
  
  %Scattering related data
  % hw = 0.0001*tf;
  
  %phonon spectrum
  hw = 1e-2;
  hwn = floor(hw/dE);
  scllim = hwn+1;
  sculim = NE - (hwn+1);
  N =  1/(exp(hw/KTs)-1); %Bose function
  A0 =0.03;  %A0 =0 ; no interactin   A0=0.05 is ok

  Dab = N*A0*eye(2*Np);
  Dem =(N+1)*A0*eye(2*Np);


%loop over every energy
   
   mu1 = Ef+0.5*V;  mu2 = Ef-0.5*V;

 
 
   Jz = 0; Jx =0; Jy = 0; %initialize spin current in Jzy Jxy Jyy (transport in y direction)
   Dos =0;
   
   Jz1 = 0; Jx1 =0; Jy1 = 0;
   Jz2 = 0; Jx2 =0; Jy2= 0;
   Jz3 = 0; Jx3 =0; Jy3 = 0;
   
   
     
     U = [0.5*V*ones(1,NL) V*linspace(0.5,-0.5,NOx) -0.5*V*ones(1,NR)] + UBB;
%    U = [0.5*V*ones(1,NL) V*linspace(0.5,-0.5,NOx) -0.5*V*ones(1,NR)] ;
     U = kron(diag(U),eye(2));
   
   
   
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
   
I1 = zeros(NE,1);
I2 = zeros(NE,1);
Is = zeros(NE,1);
Iin = zeros(Np+1,NE);
     
 %%%equlibrium part
   for k = 1:NE
       
     fL(k) = 1/(1+exp((E(k)-mu1)/kT));
     fR(k) = 1/(1+exp((E(k)-mu2)/kT)); 
     

     %%%%%%%%%%%%%%%%
%      ka1 =  acos(1- (E(k)-U(1))/2*tf);
%      ka2 = acos(1- (E(k)-exchange-U(1))/2*tf);
%      ka3 = acos(1- (E(k)-U(2*Np))/2*tf);
%      ka4 = acos(1- (E(k)-exchange-U(2*Np))/2*tf);
%      
%      
%      sig1(1,1,k) = -tf*exp(1i*ka1);
%      sig1(2,2,k) = -tf*exp(1i*ka2);
%      sig2(2*Np-1:2*Np,2*Np-1:2*Np,k) = R*[-tf*exp(1i*ka3) 0;0 -tf*exp(1i*ka4)]*R';
     %%%%%%%%%%%%%%%
     
     sig1(1,1,k) = selfdatta1D(E(k),0.5*V,tf);
     sig1(2,2,k) = selfdatta1D(E(k),0.5*V+exchange,tf);
     sig2(2*Np-1:2*Np,2*Np-1:2*Np,k) = R*[selfdatta1D(E(k),-0.5*V,tf) 0;0 selfdatta1D(E(k),-0.5*V+exchange,tf)]*R';
     
     
     gam1(:,:,k) = 1i*(sig1(:,:,k) - sig1(:,:,k)');  gam2(:,:,k) = 1i*(sig2(:,:,k) - sig2(:,:,k)');
     gam1in(:,:,k) = fL(k)*gam1(:,:,k);   gam2in(:,:,k) = fR(k)*gam2(:,:,k);
     G(:,:,k) = inv((E(k)+zplus)*eye(2*Np)-sig1(:,:,k)-sig2(:,:,k)-H-U);
     Gn(:,:,k) = G(:,:,k)*(gam1in(:,:,k) + gam2in(:,:,k))*G(:,:,k)';
     A(:,:,k) = 1i*(G(:,:,k)-G(:,:,k)');
     Gp(:,:,k) = A(:,:,k)-Gn(:,:,k);
   end
   
  
    err = 1000;
     
    while err > 1e-7
        err1 = 0;
        err2 = 0;
        for k = scllim:sculim
            Ssinnew = Dem.*diag(diag(Gn(:,:,k+hwn)))+Dab.*diag(diag(Gn(:,:,k-hwn)));
            Ssoutnew =Dab.*diag(diag(Gp(:,:,k+hwn)))+Dem.*diag(diag(Gp(:,:,k-hwn)));
            err1 = err1+0.05*(sum(sum(abs(Ssout(:,:,k)-Ssoutnew)))/...
                sum(sum(abs(Ssout(:,:,k))+abs(Ssoutnew))));
            err2 = err2+0.05*(sum(sum(abs(Ssin(:,:,k)-Ssinnew)))/...
                sum(sum(abs(Ssin(:,:,k)+abs(Ssinnew)))));
     

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
%  Ss = Ssi;
        err = err1 + err2;
        fprintf(1,'Error = %g\n',err);
    end
    
    HH = H+U;
   
    
    for k = 1:NE
    I1(k) = real(trace(gam1in(:,:,k)*A(:,:,k)-gam1(:,:,k)*Gn(:,:,k)));
    I2(k) = real(trace(gam2in(:,:,k)*A(:,:,k)-gam2(:,:,k)*Gn(:,:,k)));
    Is(k) = real(trace(Ssin(:,:,k)*A(:,:,k)-1i*(Ss(:,:,k)-Ss(:,:,k)')*Gn(:,:,k)));
    % Iin(2:Np,k) = imag((2*tf*(diag(Gn(:,:,k),1))));
    Dos(k) = 1i*(trace(G(:,:,k)-G(:,:,k)'));
       
      
    end
    
    
   
    for k = 1:NE
        
        
        
        
       
     for j = 1:Np-1
     Iin(j+1,k)=(1i)*(trace(HH(2*j-1:2*j,2*j+1:2*j+2)*Gn(2*j+1:2*j+2,2*j-1:2*j,k)-HH(2*j+1:2*j+2,2*j-1:2*j)*Gn(2*j-1:2*j,2*j+1:2*j+2,k)));
     end

    end
    
    
          
    for k = 1:NE
        
        for j = 1:Np
            dos(1:Np,k)=(1i)*(trace(G(2*j-1:2*j,2*j-1:2*j,k)-G(2*j-1:2*j,2*j-1:2*j,k)'));
        end
        
    end
    
    
    
    for k = 1:NE
        
        for j = 1:Np
            GGGn(1:Np,k)=(trace(Gn(2*j-1:2*j,2*j-1:2*j,k)));
        end
        
    end
    
   
  fprintf(1,'global current conservation = %g\n',sum(I1+I2+Is)*dE)   
    
Iin(1,:) = I1;
Iin(Np+1,:) = -I2;
L = 1:Np;
Lc = 0:Np;
  
 %plot currents
figure;
hold on;
grid on;
box on;
plot(E,I1,'b*-');
plot(E,I2,'ro-');
plot(E,Is,'gx-');
legend('Contact 1','Contact 2','Scatterer');
title('Current vs Energy ');
xlabel('Energy (eV)');
ylabel('Current (1/eV)');

% logI vs position is more good
figure;
contourf(a*Lc,E,real(transpose(log(Iin))),500,'linestyle','none'); 
% contourf(a*Lc,E,real(transpose(Iin)),500,'linestyle','none'); 
colormap(jet);
colorbar; shading interp;
 xlabel('position');
 ylabel('Energy');
 title('energy resolved current 1/eV ');
 
 
figure;
plot(E,Dos)
xlabel('Energy')
ylabel('Total-DOS------arb.unit')
 
 
 
figure;
contourf(a*L,E,real(transpose(dos)),500,'linestyle','none'); 
colormap(jet);
colorbar; shading interp;
 xlabel('position');
 ylabel('Energy');
 title('LDOS 1/eV ');
 

figure;
contourf(a*L,E,real(transpose(GGGn)),500,'linestyle','none'); 
colormap(jet);
colorbar; shading interp;
 xlabel('position');
 ylabel('Energy');
 title('electron density 1/ev ');
 
% logI vs position is more good
figure;
contourf(a*Lc,E,real(transpose(Iin)),500,'linestyle','none'); 
% contourf(a*Lc,E,real(transpose(Iin)),500,'linestyle','none'); 
colormap(jet);
colorbar; shading interp;
 xlabel('position');
 ylabel('Energy');
 title('energy resolved current 1/eV ');
 

