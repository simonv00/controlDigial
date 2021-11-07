function Gso = FTC_SOR_Completo(K, theta, t1, t2, t3)
% k = (delta salida/ delta entrada) -> La ganancia del sistema
% t1, t2, t3 -> Tiempo requerido para que la respuesta alcance el
% (15%,45%,75%) del cambio total
% Wn -> frecuencia natural de la planta
% theta -> retraso del sistema (10% del tao aprox)
% zita -> coeficiente de amortiguamiento (>1 sobre amortiguado, <1
% subamortiguado.
x = (t2-t1)/(t3-t1);
zita = (0.0805-5.547*(0.475-x)^2)/(x-0.356);
F2 = 0;
if(zita >= 1)
    F2 = 2.6*zita-0.6;
elseif (zita < 1)
    F2 = 0.708*(0.2811)^zita;
end
Wn = F2/(t2-t2);

Gso = FTC_SOR(K, theta, zita, Wn);
end

