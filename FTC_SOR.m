function Gso = FTC_SOR(K, theta, zita, Wn)
% k = (delta salida/ delta entrada) -> La ganancia del sistema
% Wn -> frecuencia natural de la planta
% theta -> retraso del sistema (10% del tao aprox)
% zita -> coeficiente de amortiguamiento (>1 sobre amortiguado, <1
% subamortiguado.
s = tf('s');
Gso = (K*Wn^2*exp(-theta*s))/(s^2+2*Wn*zita*s+Wn^2);
end

