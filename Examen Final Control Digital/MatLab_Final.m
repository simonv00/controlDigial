%% Punto 1
Yss = 35; 
Uss = 10;
Yo = 0;
to = 0;
theta = 0;
tao = 6.25; %63.2
K = (Yss-Yo)/(Uss);
num = K;
den = [tao 1];
Gs = tf (num, den, 'inputdelay', theta);

%Discretización
Ts = 0.5;
Gz = c2d(Gs,Ts)
%step (Gs, Gz)
B = Gz.Numerator{1}
A = Gz.Denominator{1}

%% Punto A
%Ecuaciones Caracteristicas Desadas
tss = 17.5;
Mp = 0.05;
Xid = -log(Mp)/sqrt(log(Mp)^2+pi^2);
Wnd = 4/(tss*Xid);
Mag=exp(-Xid*Wnd*Ts);
th=57.3*Wnd*Ts*sqrt(1-Xid^2);
polos=Mag*[cosd(th)+i*sind(th), cosd(th)-i*sind(th)];
Q=poly(polos)

% Constantes
q0 = (Q(2)-A(2)+1)/B(2)
q1 = (Q(3)+A(2))/B(2)

%% Punto B
% Constantes
R0 = 1
S0 = (Q(2)-A(2)+1)/B(2)
S1 = (Q(3)+A(2))/B(2)
T = S1+S0
%Graficación
T0 = 0;
per(1:100) = 1; 
per(101:200) = 0; %PERTURBACIÓN
nit=200;    %Numero de interacciones
u(1:nit) = 0; 
ym(1:nit) = 0;  
r(1:nit) = T0;

% Referência
r(1:60) = 10; 
r(61:120) = 40; 
r(121:200) = 30; 

%Referencia, desde tiempo inicial a final
%Tambien se debe añadir la verdadera referenacia
%r(10:80) = 45; r(81:nit) = 30; 

for k = 2:nit  %K se define de forma que nunca de un vector en la posicion 0  
    t = 0:Ts:(k-1)*Ts;
    ym=lsim(Gs,u(:,1:k),t,'zoh')+per(k);
    %Ley de Control
    u(k)=(1/R0)*(T*r(k)-S0*ym(k)-S1*ym(k-1)+R0*u(k-1));

    if u(k) > 100
        u(k) = 100;
    else if u(k)<0
            u(k) = 0;
        end
    end
        
end

t = 0:Ts:(nit-1)*Ts;
figure
subplot(2,1,1)
stairs(t,r,'--k','Linewidth',2),hold on
stairs(t,ym,'-r','Linewidth',2)
xlabel('Tiempo (s)');
ylabel('Temperatura (C)');
legend('r','y','Location','SouthEast')
grid on;
hold
subplot(2,1,2)
stairs(t,u,'b','Linewidth',2)
hold on
stairs(t,per,'m','Linewidth',2)
xlabel('Tiempo (s)');
ylabel('Heater (%)');
legend('u')
grid on;




