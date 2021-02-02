clear all
close all
clc
%% basic visualization (we look at the whole spatial locations at a timestep each time
lochist = "D:\watersteelwater10mm.pogo-hist" %location for Hist file. 
hist1 = loadPogoHist(lochist)
vals= hist1.sets.midpointy.histTraces;



%%

figure()
x= 5e-14
ylim([-x,x])
grid minor
      grid on 
      title('Wave Propagation')
      xlabel('Node numbering')
      grid on
      ylabel('Amplitude on the second DOF')
for i=1:50:20000
      grid minor
            grid on 
      title('Wave Propagation')
      xlabel('Node numbering')
      grid on
      ylabel('Amplitude on the second DOF')
      plot(vals(i,:))
      grid minor
      ylim([-x,x])
      %pause(0.2)
      w = waitforbuttonpress;
      i
      grid on 
      title('Wave Propagation in solid')
      xlabel('Node numbering')
      grid on
      ylabel('Amplitude on the second DOF')
end
