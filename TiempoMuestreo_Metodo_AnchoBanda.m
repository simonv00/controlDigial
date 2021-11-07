function [Tmin,Tmax] = TiempoMuestreo_Metodo_AnchoBanda(n,d)
% n -> numerador del sistema (vector)
% d -> numerador del sistema (vector)

[nw,dw]=cloop(n,d,-1);  %Calcula FT en lazo cerrado 
[mag,fase,w]=bode(nw,dw); %Calcula Magnitud, y fase 
mag1=mag(1,1);        % Magnitud a baja frecuencia      
mag2=0.707*mag1;    %Calcula el valor de la magnitud para wc 
wc=interp1(mag,w,mag2,'spline'); %Interpolacion para cálculo exacto 
wmin=8*wc;
wmax=12*wc;
Tmin=2*pi/wmax;
Tmax=2*pi/wmin;
%fprintf(' RANGO PARA EL PERIODO : Wmin=%3.2f   Wmax=%3.2f',wmin, wmax) 
fprintf(' RANGO PARA EL PERIODO : Tmin=%3.2f   Tmax=%3.2f',Tmin, Tmax) 

end

