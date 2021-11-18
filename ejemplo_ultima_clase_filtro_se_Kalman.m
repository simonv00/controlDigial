clc
clear
close all
%%
DatosLevitador=load('Datos_Levitador.txt');
t=DatosLevitador(:,1);
ut=DatosLevitador(:,2);
ut = [0;ut];
ybt=DatosLevitador(:,3);
vbt=DatosLevitador(:,4);
%%
%Sistema No Lineal
syms yb Vb Va V w1 w2 w3 v
f1=Vb+w1;
f2=0.1649*(Va-Vb)^2-9.81+w2;
f3=(93.9242*V+3.0098)/yb-Va/2+w3;
hx=yb+v;
f=[f1;f2;f3];
x=[yb;Vb;Va];
w=[w1;w2;w3];
A=vpa(jacobian(f,x),5);
B=vpa(jacobian(f,V),5);
C=jacobian(hx,x);
L=vpa(jacobian(f,w),5);
M=vpa(jacobian(hx,v),5);
%%
%Sistema lineal en tiempo continuo
yb=25.135; Vb=0; Va=7.713; V=1;% Punto de Equilibrio
Ac=eval(A);
Bc=eval(B);
Cc=eval(C);
sysc=ss(Ac,Bc,Cc,0);
lambdac=eig(Ac);
Taoc=1/abs(real(min(lambdac)));
Tm=Taoc*0.02147;
%%
%Sistema lineal en tiempo discreto
sysd=c2d(sysc,Tm,'zoh');
[Ad,Bd,Cd,Dd]=ssdata(sysd);
lambdad=eig(Ad);
Taod=-0.3/(log(max(abs(lambdad))));

%%
%Filtro se Kalman

Q = eye(3)*0.1;
R = 10;
xhat_pos = [3;-1;-5];
P_pos = Q;
ybhat = [];
vbhat = [];
vahat = [];

for i=1:1603
    % Etapa de predicción
    xhat_pri = Ad*xhat_pos + Bd*ut(i);
    P_pri = Ad*P_pos*Ad' + Q;
    
    % Etapa de actualización y corrección
    Kk = P_pri*Cd'*(Cd*P_pri*Cd'+R)^-1;
    xhat_pos = xhat_pri + Kk*(ybt(i) - Cd*xhat_pri)
    P_pos = (eye(3)-Kk*Cd)*P_pri;
    
    xhat1 = xhat_pos(1);
    xhat2 = xhat_pos(2);
    xhat3 = xhat_pos(3);
    
    ybhat = [ybhat xhat1];
    vbhat = [vbhat xhat2];
    vahat = [vahat xhat3]; 
end
k = 0:1:1602;
k =k*0.1;
plot(k,ybhat,'m',k,vbhat,'c',k,vahat,'k',t,ybt,'b',t,vbt,'g','linewidth',1.5)
xlim([0 160])

