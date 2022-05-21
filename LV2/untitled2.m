% kada je poznata wc
clear all; close all; clc;
s = tf('s');
w = 0.01:0.001: 110;
PMspec = 40; % zeljena fazna margina
wc = 10; % zeljena cross-over freknvecija
G = 1/(s*(s+1)); % proces

[mag,phase] = bode(G, w); % ne vrati u dB (samo plota u dB), phase je u stepenima
mag = squeeze(mag);
phase = squeeze(phase);
phaseGwc = interp1(w, phase, wc); % tamo gdje je w = wc

FImax = PMspec - (180 + phaseGwc) + 2
radijan = pi / 180;
alpha = (1 - sin(FImax * radijan)) / (1 + sin(FImax*radijan))

zc = wc * sqrt(alpha)
pc = zc / alpha

Gc = (s/zc + 1) / ( s/pc + 1);
Gtotal = Gc * G;

[mag,phase] = bode(Gtotal, w); % ne vrati u dB (samo plota u dB), phase je u stepenima
mag = squeeze(mag);

kc = 10 ^ (-interp1(w, 20*log10(mag), 10) / 20)

bode(kc*Gtotal, w)
hold on;
grid on;
bode(G,w)
legend('Gtotal','G')