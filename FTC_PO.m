function Gpo = FTC_PO(K,theta,tao)
% k = (delta salida/ delta entrada) -> La ganancia del sistema
% theta -> retraso del sistema (10% del tao aprox)
% Tao -> tiempo que le toma al sistema a llegar al 63.2% de su valor final
s = tf('s');
Gpo = K*exp(-theta*s)/(tao*s+1);
end

