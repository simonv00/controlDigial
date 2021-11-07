%% Control de velocidad
% Sistema
L = 300;
s = tf("s");
Tao = 0.31;
thetaOriginal = 0.03;
k = 6.8;
Gs = (k*exp(-0*s))/(Tao*s+1);
step(Gs)
step(feedback(1*Gs,1))
% Sistema discreto
Ts = 0.01;
Gz = c2d(Gs,Ts,'zoh');
step(Gs,Gz);

%% Ziegler y Nichols tunning
% Basic parameters
T = 0.01;
theta = thetaOriginal + (T/2);

% Control P
Kc_P1 = (Tao)/(k*theta);

q0_P1 = Kc_P1;

% Control PI
Kc_PI1 = (9*Tao)/(10*k*theta);
ti_PI1 = (10*theta)/3;

q0_PI1 = Kc_PI1*(1+(T/(2*ti_PI1)));
q1_PI1 = -Kc_PI1*(1-(T/(2*ti_PI1)));

% Control PID
Kc_PID1 = (6*Tao)/(5*k*theta);
ti_PID1 = 2*theta;
td_PID1 = 0.5*theta;

q0_PID1 = Kc_PID1*(1+T/(2*ti_PID1)+td_PID1/T);
q1_PID1 = -Kc_PID1*(1-T/(2*ti_PID1)+td_PID1/T);
q2_PID1 = Kc_PID1*td_PID1/T;

%C_1 = (q0_PID1*z^2+q1_PID1*z+q2_PID1)/(z^2-z)
%feedback(C_1*Gz,1)

%% Ziegler y Nichols limit gain
% Basic parameters
bode(Gs/(Gs+1))
%Gm es la ganancia limite, Pm es phase margin, Wcg Phase crossover
%frequency, Wcp Gain crossover frequency
[Gm,Pm,Wcg,Wcp] = margin(Gs);

%% Asignacion de polos
% Basic parameters
T = 0.01;
Gz = c2d(Gs,T,'zoh');
step(Gs,Gz)
b1 = Gz.Numerator{1}(2);
a1 = Gz.Denominator{1}(2);
ts = 10;
Mp = 0.1;
Xi = -log(Mp)/sqrt(log(Mp)^2 + pi^2);
Wn = 4/(ts*Xi);
Mag = exp(-Xi*Wn*T);
Ang = Wn*T*sqrt(1-Xi^2);
Pol1 = Mag*(cos(Ang)+1i*sin(Ang));
Pol2 = conj(Pol1);
Ecd = conv([1 -Pol1],[1 -Pol2]);

% Control PI
q0_PI3 = (Ecd(2)+1-a1)/b1;
q1_PI3 = (Ecd(3)+a1)/b1;

% Control PID
Pol3 = -0.001;
Ecd2 = conv (Ecd, [1 -Pol3])
q0_PID3 = (Ecd2(2)+1-a1)/b1
q1_PID3 = (Ecd2(3)+a1)/b1
q2_PID3 = Ecd2(4)/b1

%% Metodo de la integral del cuadrado del error (ICE)
% Basic parameters
T = 0.01;
theta = thetaOriginal + (T/2);

%Control P
a = 1.411;
b = -0.917;
Kc_P4 = (a/k)*(theta/Tao)^b;

q0_P4 = Kc_P4;

%Control PI
a = 1.305;
b = -0.959;
Kc_PI4 = (a/k)*(theta/Tao)^b;

a = 0.492;
b = 0.739;
ti_PI4 = (Tao/a)*(theta/Tao)^b;

q0_PI4 = Kc_PI4*(1+(T/(2*ti_PI4)));
q1_PI4 = -Kc_PI4*(1-(T/(2*ti_PI4)));

%Control PID
a = 1.495;
b = -0.945;
Kc_PID4 = (a/k)*(theta/Tao)^b;

a = 1.101;
b = 0.771;
ti_PID4 = (Tao/a)*(theta/Tao)^b;

a = 0.56;
b = 1.006;
td_PID4 = (Tao*a)*(theta/Tao)^b;

q0_PID4 = Kc_PID4*(1+T/(2*ti_PID4)+td_PID4/T);
q1_PID4 = -Kc_PID4*(1-T/(2*ti_PID4)+td_PID4/T);
q2_PID4 = Kc_PID4*td_PID4/T;

%% Metodo de la integral del valor absoluto del error (IAE)
% Basic parameters
T = 0.01;
theta = thetaOriginal + (T/2);

%Control P
a = 0.902;
b = -0.985;
Kc_P5 = (a/k)*(theta/Tao)^b;

q0_P5 = Kc_P5;

%Control PI
a = 0.984;
b = -0.986;
Kc_PI5 = (a/k)*(theta/Tao)^b;

a = 0.608;
b = 0.707;
ti_PI5 = (Tao/a)*(theta/Tao)^b;

q0_PI5 = Kc_PI5*(1+(T/(2*ti_PI5)));
q1_PI5 = -Kc_PI5*(1-(T/(2*ti_PI5)));

%Control PID
a = 1.435;
b = -0.921;
Kc_PID5 = (a/k)*(theta/Tao)^b;

a = 0.878;
b = 0.749;
ti_PID5 = (Tao/a)*(theta/Tao)^b;

a = 0.482;
b = 1.137;
td_PID5 = (Tao*a)*(theta/Tao)^b;

q0_PID5 = Kc_PID5*(1+T/(2*ti_PID5)+td_PID5/T);
q1_PID5 = -Kc_PID5*(1-T/(2*ti_PID5)+td_PID5/T);
q2_PID5 = Kc_PID5*td_PID5/T;

%% Metodo de la integral del error absoluto del error por el tiempo (IAET)
% Basic parameters
T = 0.01;
theta = thetaOriginal + (T/2);

%Control P
a = 0.94;
b = -1.084;
Kc_P6 = (a/k)*(theta/Tao)^b;

q0_P6 = Kc_P6;

%Control PI
a = 0.859;
b = -0.977;
Kc_PI6 = (a/k)*(theta/Tao)^b;

a = 0.674;
b = 0.68;
ti_PI6 = (Tao/a)*(theta/Tao)^b;

q0_PI6 = Kc_PI6*(1+(T/(2*ti_PI6)));
q1_PI6 = -Kc_PI6*(1-(T/(2*ti_PI6)));

%Control PID
a = 1.357;
b = -0.947;
Kc_PID6 = (a/k)*(theta/Tao)^b;

a = 0.842;
b = 0.738;
ti_PID6 = (Tao/a)*(theta/Tao)^b;

a = 0.381;
b = 0.995;
td_PID6 = (Tao*a)*(theta/Tao)^b;

q0_PID6 = Kc_PID6*(1+T/(2*ti_PID6)+td_PID6/T);
q1_PID6 = -Kc_PID6*(1-T/(2*ti_PID6)+td_PID6/T);
q2_PID6 = Kc_PID6*td_PID6/T;

%% Metodo 3c
% Basic parameters
T = 0.01;
theta = thetaOriginal + (T/2);

% Control PI
Kc_PI7 = (0.928/k)*(theta/Tao)^-0.946;
ti_PI7 = 0.928*Tao*(theta/Tao)^0.503;

q0_PI7 = Kc_PI7*(1+(T/(2*ti_PI7)));
q1_PI7 = -Kc_PI7*(1-(T/(2*ti_PI7)));

% Control PID
Kc_PID7 = (1.307/k)*(theta/Tao)^-0.95;
ti_PID7 = 0.74*Tao*(theta/Tao)^0.738;
td_PID7 = 0.365*Tao*(theta/Tao)^0.95;

q0_PID7 = Kc_PID7*(1+T/(2*ti_PID7)+td_PID7/T);
q1_PID7 = -Kc_PID7*(1-T/(2*ti_PID7)+td_PID7/T);
q2_PID7 = Kc_PID7*td_PID7/T;

%% Metodo de Cohen-Coon

% Basic parameters
T = 0.01;
theta = thetaOriginal + (T/2);

% Control P
Kc_P8 = (Tao/(k*theta))*(1+theta/(3*Tao));

q0_P8 = Kc_P8;

% Control PI
Kc_PI8 = (Tao/(theta*k))*(0.9+theta/(12*Tao));
ti_PI8 = (theta*(30*Tao+3*theta))/(9*Tao+20*theta);

q0_PI8 = Kc_PI8*(1+(T/(2*ti_PI8)));
q1_PI8 = -Kc_PI8*(1-(T/(2*ti_PI8)));

% Control PID
Kc_PID8 = (Tao/(theta*k))*((4/3)+theta/(4*Tao));
ti_PID8 = (theta*(32*Tao+6*theta))/(13*Tao+8*theta);
td_PID8 = (4*theta*Tao)/(11*Tao+2*theta);

q0_PID8 = Kc_PID8*(1+T/(2*ti_PID8)+td_PID8/T);
q1_PID8 = -Kc_PID8*(1-T/(2*ti_PID8)+td_PID8/T);
q2_PID8 = Kc_PID8*td_PID8/T;

%% Metodo LGR
rlocus(Gz)

% Control P
% 0.656 < Kc < 9
Kc_P9 = 0.656;

q0_P9 = Kc_P9;

%% Controlador DeadBeat

% De orden normal (m)
B = Gz.Denominator{1};
b1 = B(1);
b2 = B(2);
% revisar la z en el numerador
z = tf('z',Ts);
A = Gz.Numerator{1};
Ap = A*z^-1;

q0_DBM = 1/(b1+b2);

Dz_M10 = (q0_DBM*Ap(2))/(1-q0_DBM*(B(1)+B(2)*z^-1));

% De orden incremental (m+1)

q0_DBI = 1/((1-A(2))*(b1+b2));
alpha_10 = 1/A(2);
Dz_I10 = (q0_DBI*Ap(2)*(1-z^-1/alpha_10))/(1-q0_DBI*(B(1)+B(2)*z^-1)*(1-z^-1/alpha_10));

%% Metodo Dalhin
feedback(Gs,1);
Teq = 0.31/7.8;
% 1.25/Teq <= lambda <= 4/Teq
lambda = 3/Teq;
% N = enterero del numero
N = thetaOriginal/Ts;
Dz_D11 = (1-exp(-lambda*Ts))*z^(-N-1)/(1-exp(-lambda*Ts)*z^(-1)-(1-exp(-lambda*Ts))*z^(-N-1));

%% RST incremental
% Basic parameters
Ts = 0.01;
Gz = c2d(Gs,Ts,'zoh');
step(Gs,Gz)
a0 = 1
b0 = Gz.Numerator{1}(2)
a1 = Gz.Denominator{1}(2)
d = 1

ts = 5;
Mp = 0.1;
Xi = -log(Mp)/sqrt(log(Mp)^2 + pi^2);
Wn = 4/(ts*Xi);
Mag = exp(-Xi*Wn*Ts);
Ang = Wn*Ts*sqrt(1-Xi^2);
Pol1 = Mag*(cos(Ang)+1i*sin(Ang));
Pol2 = conj(Pol1);
Ecd = conv([1 -Pol1],[1 -Pol2]);
B = Gz.numerator{1};
A = Gz.denominator{1};
a1 = A(2);
b1 = B(2);
R0 = 1;
S0 = 0.2932;
S1 = -0.277;
Tt = S0+S1;

% inicializa parametros de Simulacion
T0 = 0;  %Temperatura Ambiente
nit=500;    %Numero de interacciones
u(1:nit) = 0; 
ym(1:nit) = 0;  
r(1:nit) = T0;

% Referência
r(10:nit) = 300; 

for k=4:nit
        t = 0:Ts:(k-1)*Ts;
        ym=lsim(Gs,u(:,1:k),t,'zoh')+T0;
        %Ley de Control
        u(k)=(1/R0)*(Tt*r(k)-S0*ym(k)-S1*ym(k-1)+R0*u(k-1));

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
ylabel('Velocidad (RPM)');
legend('r','y','Location','SouthEast')
grid on;
hold
subplot(2,1,2)
stairs(t,u,'b','Linewidth',2)
xlabel('Tiempo (s)');
ylabel('PWM (%)');   
legend('u')
grid on;

%% Control de angulo
% Sistema
s = tf("s");
Tao = 0.31;
thetaOriginal = 0.01;
k_A = 6.8;
U = 40;
Gs_A = (k_A*exp(-0*s))/(s*(Tao*s+1));
Gz_A = c2d(Gs_A,0.06)
rlocus(Gz_A)
step(Gs_A,Gz_a)

%% Asignacion de polos

%% Metodo LGR
rlocus(Gz_A)
KP_LGR_A = 0.339