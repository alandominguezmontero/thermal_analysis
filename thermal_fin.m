clear all
close all

% Geometrical data
w=5e-3;
t=1e-3;
L=0.1;
% Temperatures
Tb=100
Tinf=20
% Thermodynamic properties
h=25;
k=385;

%% Calculate area and perimeter

P=2*w+2*t;
Ac=w*t;


%% Solve equation

m=sqrt(h*P/(k*Ac))
m2=m*m;
hmk=h/(m*k);

Nx=10; xa= linspace(0,L,Nx)';

Ta=Tinf+(Tb-Tinf)*(cosh(m*(L-xa)) + ...
    hmk*sinh(m*(L-xa)))/(cosh(m*L)+hmk*sinh(m*L))

M = sqrt(k*Ac*h*P)*(Tb-Tinf);
qf=M.*(sinh(m*L) + hmk.*cos(m*L))./(cosh(m*L)+hmk*sinh(m*L))



%% Using finite differences scheme

N=40;
dx=L/N;
dx2=dx*dx;
xn=dx:dx:L;
xn=xn';
A=zeros(N,N);
b=zeros(N,1);

% Node 1
A(1,1) = -(2+m2*dx2);
A(1,2) = 1;
b(1)= -m2*dx2*Tinf -Tb;

% Interior nodes
for i = 2:N-1
    A(i,i-1) = 1;
    A(i,i)=-(2+m2*dx2);
    A(i,i+1)=1;
    b(i)=-m2*dx2*Tinf;
end

% Node n

A(N,N-2)=1;
A(N,N) = -(2*m2*dx2 + 1 + 2*h*dx/k);
b(N)=-(2*h*dx/k + 2*m2*dx2)*Tinf;


Tn=A\b;

xn=[0;xn];
Tn=[Tb;Tn];

%% Plot results

plot(xa,Ta,'b*')
hold on
plot(xn,Tn,'r-','LineWidth',2)
ylabel('Temperature [ºC]')
xlabel('X Position [m]')
grid on
legend('Eq sol','Finite differences sol')
title('Fin Heat Equation Solution')
