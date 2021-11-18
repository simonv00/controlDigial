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
B = [1/La 0]; % Matriz de entrada
C =[1 0]; % Matriz de salida
D = 0; % Matriz de transferencia directa
%% 
% Validación de modelo 

% Datos de corriente

Motor = load ('current.txt');
tc = Motor(:,1);
Va_t = Motor(:,2);
Ia_t = Motor(:,3);
% Modelo

% Condiciones iniciales de entrada y estados

u=[0 0]';

xv=pinv(A)*(-B*u);


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
Ia(1)=xv(1);        
wa(1)=xv(2);   

for i=1:n 

    cambio_xv = A*xv(:,i) + B*u;
    
    xv(:,i+1) = xv(:,i) + delta*cambio_xv;                    
    
    if i>T_escalon
        u(1) = Va;
        
     else
        u(1) = 0; 
     end
end

%%inicializamos las variables del modelo a estimar
u=[0 0]';
xs=pinv(A)*(-B*u);

%%

%measurement error matrix
R=cov(Motor(2500:4000,3));

%Error dynamic matrix (model uncertainties)
Q=100*R*eye(2);

%Error covariance matrix
P=1000*eye(2);
%%


for i=1:n 
    
    %cargamos el modelo del sistema
    cambio_xs = A*xs(:,i) + B*u;
    xs(:,i+1) = xs(:,i) + delta*cambio_xs;
    
    if i>T_escalon
        u(1) = Va;
     else
        u(1) = 0;
    end
    
    P=A*P*A'+Q;
    %Kalman gain
    K = P*C'*inv(C*P*C'+R);
    %State estimation update
    xs(:,i+1) = xs(:,i+1) + K*(Motor(i,3)- C*xs(:,i+1));
    %Covariance matrix update
    P = (eye(2) - K*C)*P;
end

%% Model simulation

% Gráficas de validacion del modelo

% Comparación para corriente
t=0:delta:n*delta;
t2=0:delta2:n2*delta2;
figure (1);
plot(t,xs(1,:),'b','LineWidth',3);
hold on
plot(tc,Ia_t);
grid on
title('CORRIENTE DE ARMADURA')
xlabel('Tiempo (ms)')
ylabel('Corriente (A)')
legend('Datos estimados','Datos experimentales')
%comparacion modelo velocidad
figure(2)
plot(t,xs(2,:),'b','LineWidth',4);
hold on
plot(tc,xv(2,:),'lineWidth',2);
grid on
title('VELOCIDAD ANGULAR ')
xlabel('Tiempo (ms)')
ylabel('velocidad rad/s')
legend('Estimado','Modelo')
