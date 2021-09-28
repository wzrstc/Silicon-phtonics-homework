clc
clear all
close all

lam0 = 1.55e-6;             %wavelength
k0 = 2*pi/lam0;             %k vector
nf = 3.5;                   %core 
ns = 1.5;                   %lower clad
nc = 1;                     %upper clad
D = 3e-6;                   
H = 5e-6;
W = 3.5e-6;
m=0;                        %foudamental mode
initial1 = 3.4966*k0;       %initial value

%%
%solve foundamental TM mode
funTM_clad = @(beta)sqrt(k0^2*nf^2-beta^2)*D-m*pi-atan((nf^2/ns^2)*(sqrt((beta^2-k0^2*ns^2)/(k0^2*nf^2-beta^2))))...
-atan((nf^2/nc^2)*(sqrt((beta^2-k0^2*nc^2)/(k0^2*nf^2-beta^2))));
out = fzero(funTM_clad,initial1);

neff_TM_clad = out/k0;

%the only difference between funTM_core and funTM_clad is D and H
funTM_core = @(beta)sqrt(k0^2*nf^2-beta^2)*H-m*pi-atan((nf^2/ns^2)*(sqrt((beta^2-k0^2*ns^2)/(k0^2*nf^2-beta^2))))...
-atan((nf^2/nc^2)*(sqrt((beta^2-k0^2*nc^2)/(k0^2*nf^2-beta^2))));
out = fzero(funTM_core,initial1);

neff_TM_core = out/k0;


%%
%solve foundamental TM mode
initial2 = [3.4906*k0 3.49659397*k0];
funTE = @(beta)sqrt(k0^2*neff_TM_core^2-beta^2)*W-m*pi-atan(sqrt((beta^2-k0^2*neff_TM_clad^2)/(k0^2*neff_TM_core^2-beta^2)))...
-atan(sqrt((beta^2-k0^2*neff_TM_clad^2)/(k0^2*neff_TM_core^2-beta^2)));
out = fzero(funTE,initial2);

neff_TE = out/k0;

neff_TM0 = neff_TE;