%% 
% Código base para la programación de un filtro de Kalman para la
% estimación de velocidad a partir de la medición de corriente de un motor
% DC como parte del examen final del curso Control Digital 
% Por: Mario Alejandro Giraldo Vásquez y Juan David Núñez López
% Última edición: 07/11/2021

clear
close all
clc

%% Definición de parámetros motor DC

Ra=99.5833;                 %Resistencia de armadura [Ohmios]
La=0.58378;                 %Inductancia de armadura [Henrios]
Kv=0.09168;                 %Constante de fuerza contraelectromotriz [Voltios/(Rad*Web)]
J=0.001537;                 %Momento de inercia [Kg*m^2]
Beta=0.000752;              %Coeficiente de fricción viscosa [N*m*s/Rad]
Va=12;                      %Voltaje de armadura [Voltios]       
TL=0;                       %torque de carga [N*m]

%% Modelo en espacio de estados

% Los estados son corriente de armadura Ia [A] y velocidad angular Wa [rad/s]
% La entrada es voltaje Va [V]
% La salida medida es la corriente de armadura Ia [A]

A = [-Ra/La -Kv/La;Kv/J -Beta/J]; % Matriz de transición del estado
B = [1/La;0]; % Matriz de entrada
C =[1 0]; % Matriz de salida
D = 0; % Matriz de transferencia directa
%% 
% Validación de modelo 

% Datos de corriente

Motor = load ('current.txt');
tc = Motor(:,1);
Va_t = Motor(:,2);
ut = [0;Va_t];
Ia_t = Motor(:,3);

% Modelo

% Condiciones iniciales de entrada y estados
u=0;
xs=pinv(A)*(-B*u);

% Tiempo de simulacion (ms)
delta=0.00625;            %Período de muestreo para Ia
delta2=0.00125;           %Período de muestreo para wa
tiempo_simulacion=25;
tiempo_simulacion2=5;
n=(tiempo_simulacion/delta);  
n2=(tiempo_simulacion2/delta2);

% Las mediciones de corriente se hacen en 25ms, las mediciones de velocidad
% se hacen en 5ms

% Tiempo de ingreso del escalon (ms)
t_escalon_Va=6.3;
T_escalon=(t_escalon_Va/delta);
Ia(1)=xs(1);        
wa(1)=xs(2);     

for i=1:n 

    cambio_xs = A*xs(:,i) + B*u;
    
    xs(:,i+1) = xs(:,i) + delta*cambio_xs;                    
    
    if i>T_escalon
        u(1) = Va;
     else
        u(1) = 0; 
     end
end

% Gráficas de validacion del modelo

% Comparación para corriente
t=0:delta:n*delta;
t2=0:delta2:n2*delta2;
figure (1);
plot(t,xs(1,:),'r','LineWidth',3);
hold on
plot(tc,Ia_t);
grid on
title('CORRIENTE DE ARMADURA')
xlabel('Tiempo (ms)')
ylabel('Corriente (A)')
legend('Modelo','Datos experimentales')


%% Filtro se Kalman

sysc=ss(A,B,C,0);
lambdac=eig(A); %Polos del sistema en continuo son los valores propios
%Lambda más pequeño es más lento
Taoc=1/abs(real(min(lambdac)));
Tm=Taoc*0.02147;
%%
%Sistema lineal en tiempo discreto
sysd=c2d(sysc,Tm,'zoh');
[Ad,Bd,Cd,Dd]=ssdata(sysd);
lambdad=eig(Ad);
Taod=-0.3/(log(max(abs(lambdad))));


Q = eye(2)*4;
R = 30;
xhat_pos = [0;0];
P_pos = Q;
Ibhat = [];
Wbhat = [];

for i=1:size(tc,1)
    % Etapa de predicción
    xhat_pri = Ad*xhat_pos + Bd*ut(i);
    P_pri = Ad*P_pos*Ad' + Q;
    
    % Etapa de actualización y corrección
    Kk = P_pri*Cd'*(Cd*P_pri*Cd'+R)^-1;
    xhat_pos = xhat_pri + Kk*(Ia_t(i) - Cd*xhat_pri);
    P_pos = (eye(length(Ad))-Kk*Cd)*P_pri;
    
    xhat1 = xhat_pos(1);
    xhat2 = xhat_pos(2);
    
    Ibhat = [Ibhat xhat1];
    Wbhat = [Wbhat xhat2];
end
k = 0:1:4000;
k =k*0.00625; %Pasa a segundos
plot(k,Ibhat,'m','linewidth',1.5)
hold on
plot(tc,Ia_t,'g','linewidth',1.5)

plot(k,Wbhat,'c','linewidth',1.5)

xlim([0 25])


