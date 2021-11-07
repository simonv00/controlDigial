function LGR(n,d,t)
% n -> numerador del sistema (vector)
% d -> numerador del sistema (vector)
[nd,dd]=c2dm(n ,d,t,'zoh');    %discretizacion con t = 1
x=0 :0.1:2*pi;             
figure(1) 
plot(sin(x),cos(x),'.')        %dibuja el circulo unitario
hold
rlocus(nd ,dd)                 %Grafica del lugar geométrico de las raíces
axis([-3 1.5 -2 2])            %Escala para los ejes

end

