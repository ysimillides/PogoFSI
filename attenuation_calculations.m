clear all
close all
%%
lochist = "D:\SIMULATIONS\crunch_scripts\pogofluidw.pogo-hist" %location for Hist file. 
hist1 = loadPogoHist(lochist)
vals= hist1.sets.midpointy.histTraces;
%vals = vals(1000:22001,1000:3000);
[FreqMHz2,atten2] = importfile3('a0atten.txt',4,685);

d = size(vals);

%% visualizes the time dependent spatial positions
close all
figure()
x= 5e-14
ylim([-x,x])
grid minor
for i=1:10:8000
      grid minor
      plot(vals(i,:))
      grid minor
      ylim([-x,x])
      %pause(0.2)
      w = waitforbuttonpress;
      i
      grid on 
end
msg = "done"

%%
close all

valsnew = vals(3881:3991,2075:2080); %this worked for aaa


mm = 1
nn = 2
dxx =hist1.sets.midpointx.nodePos(1,nn)-hist1.sets.midpointx.nodePos(1,mm)
pos4 = valsnew(:,mm);
pos5 = valsnew(:,nn);
figure()
plot(pos4)
title('First Signal, Rectangular Windowed')

figure()
plot(pos5)
title('Second Signal, Rectangular Windowed')

sizenew = size(pos4)
%%
Y4 = fft(pos4,2^(nextpow2(length(pos4))+4));
L4 = size(Y4)
L4 = L4(1)
P24 = abs(Y4/L4);
P14 = P24(1:L4/2+1);
P14(2:end-1) = 2*P14(2:end-1);
f = Fs*(0:(L4/2))/L4;
plot(f,P14) 
title('Single-Sided Amplitude Spectrum of A0 at first Measurement')
xlabel('f (Hz)')
ylabel('Amplitude of FFT')
xlim([0,0.3e7])
%%
figure()

Y5 = fft(pos5,2^(nextpow2(length(pos5))+4));
L5 = size(Y5)
L5 = L5(1)
P25 = abs(Y5/L5);
P15 = P25(1:L5/2+1);
P15(2:end-1) = 2*P15(2:end-1);
f = Fs*(0:(L5/2))/L5;
plot(f,P15) 
title('Single-Sided Amplitude Spectrum of A0 at second Measurement')
xlabel('f (Hz)')
ylabel('Amplitude of FFT')
xlim([0,0.3e7])
%%
figure()
plot(f,P14)
hold on
plot(f,P15)
title('Amplitude Spectrum of X(t), both signals')
xlabel('f (Hz)')
xlim([0,0.3e7])
%ylim([0,1.5e-13])
legend(["First Signal","Second Signal"])
grid on
ylabel('Amplitude')

%%
figure()
calc = abs(P14)./abs(P15);
attenuation = log(calc);
attenuation = attenuation / (dxx);
attenuation = attenuation * 20*log10(exp(1));
plot(f,attenuation)
title('Attenuation of signal  in two consecutive points')
xlabel('f (Hz)')
grid on
ylabel('Attenuation (dB/m)')
hold on

plot([FreqMHz2*1e6],[atten2])
legend(["Calculated Attenuation","Disperse Attenuation"])

xlim([0.04e7,0.16e7])
%ylim([450,550])
