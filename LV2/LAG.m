% LAG
clear all; close all; clc;
PMspec = 40 + 10;
s = tf('s');
w = 0.01:0.001: 110;
Gp = 60/((s+1)*(s+2)*(s+3));
G = Gp * 2 / s;
Kc = 190;
G = Kc*0.1/((0.2*s+1)*(0.02*s+1))
[mag,phase] = bode(G, w); % ne vrati u dB (samo plota u dB), phase je u stepenima
mag = squeeze(mag);
phase = squeeze(phase);


 % 1.korak => pronaci wx gdje je angle(G) = -180 + PMspec
 wx = interp1( phase, w, -180 + PMspec)
 
 % 2.korak => pronaci |G(jwx)| 
 wx = 10;
 modulGwx = interp1(w, 20*log10(mag), wx)
 
 % 3.korak => izracunati alpha
 
 alpha = 10 ^ ( modulGwx / 20)
 
 % 4. izracunati zc i pc, zc < wx/10
 
 zc = wx / 10
 pc = zc / alpha
 
 Gc = (s/zc + 1) / ( s/pc + 1);
Gtotal = Gc * G;
bode(Gtotal, w)
hold on; grid on;
bode(G, w)
legend('Gtotal', 'G')


