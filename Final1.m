%% Punto 1
% Datos
Yo = 0;
Yf = 35;
Uo = 0;
Uf = 10;
K = (Yf-Yo)/(Uf-Uo);
To = 0;
theta = 0;
tao = 6.25; %63.2
% Funcion continua
num = [0 K];
den = [tao 1];
Gs = tf (num, den, 'inputdelay', theta)
% Tiempo de muestreo
[Tmin, Tmax] = TiempoMuestreo_Metodo_AnchoBanda(num,den);
% Discretizacion de la señal
Gz = c2d(Gs, Tmin, 'zoh')
% Graficacion de las señales
step(Gs, Gz)
% Extrar informacion de la señal
B = Gz.Numerator{1}
A = Gz.Denominator{1}
%Tmin = 0.5
%% Punto A
% Comportamiento deseado
teq = 0.7*tao;
Mp = 0.05;
% Frecuencia natural y coeficiente de amortiguamiento
Xid = -log(Mp)/sqrt(log(Mp)^2+pi^2);
Wnd = 1/(Xid*teq);
tss = 4/(Wnd*Xid); % 2%   
%5%->3/(Wnd*Xid)
%Wnd = 4/(tss*Xid);
% Magnitud y angulo
Mag=exp(-Xid*Wnd*Tmin);
th=57.3*Wnd*Tmin*sqrt(1-Xid^2);
% Calculo de los polos
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
    t = 0:Tmin:(k-1)*Tmin;
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

t = 0:Tmin:(nit-1)*Tmin;
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
