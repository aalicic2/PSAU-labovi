% koraci iz predavanja
clear all; close all; clc;
s = tf('s');
w = 0.01:0.001: 110;
PMspec = 50; % zeljena fazna margina
Gp = 2/((s+1)*(s+2)*(s+3)); % proces

% 1.korak => rucno odredimo broj integratora i pojacanje kc na osnovu e(t) i dobijemo
% novi sistem koji posmatramo G
kc = 2.5;  
G = Gp * kc / s;

% 2.korak => očitamo fazu nekompenziranog sistema i wx
[mag,phase] = bode(G, w); % ne vrati u dB (samo plota u dB), phase je u stepenima
mag = squeeze(mag);
phase = squeeze(phase);
wx = interp1(20*log10(mag), w, 0) % tamo gdje je 20log10(mag) = 0
phaseGwx = interp1(w, phase, wx); % tamo gdje je w = wx
PMuncom = 180 + phaseGwx % faza nekompenziranog sistema

% 3. korak  => odredimo FImax i alpha
FImax = PMspec - PMuncom + 10
radijan = pi / 180;
alpha = (1 - sin(FImax * radijan)) / (1 + sin(FImax*radijan))

% 4. korak => trazimo wMax
wMax = interp1(20*log10(mag), w, -10*log10(1/alpha))

% 5.korak => odraditi pc i zc jer je wMax = sqrt(pc*zc) i alpha = zc/pc

zc = wMax * sqrt(alpha)
pc = zc / alpha


Gc = (s/zc + 1) / ( s/pc + 1);
Gtotal = Gc * G;
bode(Gtotal, w)




