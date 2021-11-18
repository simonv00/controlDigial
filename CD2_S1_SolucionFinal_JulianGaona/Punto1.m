clc
clear all
close all

s=tf('s')
G=3.5/(6.25*s+1);
Gc=feedback(G,1)
step(Gc)
Gz=c2d(G,0.2715)
wn=0.265;
Ts=0.2715;
zita=0.69;
Pz=exp(-zita*wn*Ts);
Rp=Pz*cos(wn*Ts*sqrt(1-zita^2));
Ip=Pz*sin(wn*Ts*sqrt(1-zita^2));
c=Rp^2+Ip^2;
P1=-2*exp(-zita*wn*Ts)*cos(wn*Ts*sqrt(1-zita^2))
P2=exp(-2*zita*wn*Ts)


