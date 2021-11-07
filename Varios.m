clc 
clear
close all
%% Observador de lunenberg
A = [1.5114 -0.54881; 1 0]';
B = [0.03602; 0.0295];
C = [1 0];
sys = ss(A,B,C,0,10);

L = [0.3724; -0.9496];
xk = [0;0];
uk = 1;
ykk = [];
xk_est = [0;0];
ek = 0;
xkk1_est = []
xkk2_est = []
x1kest = 0;
x2kest = 0;
ekk = [];
for k = 0:59
    xk_1 = A*xk +B*uk;
    yk = C*xk;
    xk = xk_1;
    ykk = [ykk; yk];
    
    % obsservador
    
    xk_1est = A*xk_est + B*uk +L*ek;
    yk_est = C*xk_est;
    ek = yk-yk_est;
    xk_est = xk_1est;
    x1k_est = xk_est(1);
    x2k_est = xk_est(2);
    xkk1_est = [xkk1_est; x1k_est];
    xkk2_est = [xkk2_est; x2k_est];
    ekk = [ekk; ek];
    
end

t = 0:1:59;
plot(t*10,xkk1_est, t*10,ykk,t*10,xkk2_est)
figure(2)
plot(t*10,ekk);

%% Calculo de los tiempos de muestreo

n=input('ENTRE EL NUMERADOR DEL SISTEMA=');
d=input('ENTRE EL DENOMINADOR DEL SISTEMA=');
[nw,dw]=cloop(n,d,-1);  %Calcula FT en lazo cerrado 
[mag,fase,w]=bode(nw,dw); %Calcula Magnitud, y fase 
mag1=mag(1,1);        % Magnitud a baja frecuencia      
mag2=0.707*mag1;    %Calcula el valor de la magnitud para wc 
wc=interp1(mag,w,mag2,'spline'); %Interpolacion para cálculo exacto 
wmin=8*wc;
wmax=12*wc;
Tmin=2*pi/wmax;
Tmax=2*pi/wmin;
fprintf(' RANGO PARA EL PERIODO : Wmin=%3.2f   Wmax=%3.2f',wmin, wmax) 
fprintf(' RANGO PARA EL PERIODO : Tmin=%3.2f   Tmax=%3.2f',Tmin, Tmax) 

%% Programa para obtener el LGR
clc 
n=[1];
d=[1 0.5 0];                   %planta continua
[nd,dd]=c2dm(n ,d,1,'zoh');    %discretizacion con t = 1
x=0 :0.1:2*pi;             
figure(1) 
plot(sin(x),cos(x),'.')        %dibuja el circulo unitario
hold
rlocus(nd ,dd)                 %Grafica del lugar geométrico de las raíces
axis([-3 1.5 -2 2])            %Escala para los ejes

%% Diagrama de bode
clc
n=input( 'Entre el numerador continuo n= '); 
d=input( 'Entre el denominador continuo d= '); 
theta=input ('Entre el retardo theta= '); 
T=input( 'Entre el periodo de muestreo T= '); 
G=tf (n,d, 'IODelay',theta )
GD=c2d(G ,T) % Discretiza la funcion
get(GD);     % Muuestra propiedades de la funcion
w=0.01:0.05:3; % Rango de frecuencia deseado (opcional) 
[mag,fase,w]=bode (GD,w); %Calcula la magnitud y el ángulo 
imargin(mag,fase,w) % Hace la gráfica"
grid
[Kmax,PHIPM,wpi,wc]= imargin(mag,fase,w) % Calcula Kmax"

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

%% Ziegler y Nichols limit gain
% Basic parameters
bode(Gs/(Gs+1))
%Gm es la ganancia limite, Pm es phase margin, Wcg Phase crossover
%frequency, Wcp Gain crossover frequency
[Gm,Pm,Wcg,Wcp] = margin(Gs);

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
