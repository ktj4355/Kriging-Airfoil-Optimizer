%% Airfoil Optimizer
% Regression and Kriging Model
% kim tae jong  |   ktj4355@gmail.com   |   010 4355 1390
% Sejong University |  Propulsion Aerodynamic Lab.

close all; clear; clc
disp("LHS Sampling")
disp("  Please input Number of Sample Data      (1/5)")
sampleNum=input(">>>  ")
clc;
disp("LHS Sampling")
disp("  Input Initial AF Camber                 (2/5)")
Target_camber=input(">>>  ")
clc;
disp("LHS Sampling")
disp("  Input Initial AF Camber position        (3/5)")
Target_camberPos=input(">>>  ")
clc;
disp("LHS Sampling")
disp("  Input Initial AF Thickness              (4/5)")
Target_Thickness=input(">>>  ")
clc;
disp("LHS Sampling")
disp("  Input Initial AF Thickness position     (5/5)")
Target_ThickPos=input(">>>  ")
clc;

CamberRange = [0 10] ;      %x1
CamberPosRange = [Target_camberPos-15 Target_camberPos+15];   %x2
if Target_camberPos-15<10;CamberPosRange(1)=10;end
if Target_camberPos+15>90;CamberPosRange(2)=90;end



ThickRange = [7.5 Target_Thickness*1.5];        %x3
ThickPosRange = [Target_ThickPos-15 Target_ThickPos+15];    %x4
LECircleRange=[0.8 1.2];      %x50
if Target_ThickPos+15>90;ThickPosRange(2)=90;end


NomalF=@(x,Lb,Ub) (x-Lb)./(Ub-Lb);
unNomalF=@(k,Lb,Ub) Lb+k.*(Ub-Lb);

NomSample=lhsdesign(sampleNum,5);

Camber_Sampel_Nomal=NomSample(:,1);
CamberPos_Sampel_Nomal=NomSample(:,2);
Thk_Sampel_Nomal=NomSample(:,3);
ThkPos_Sampel_Nomal=NomSample(:,4);
LEc_Sampel_Nomal=NomSample(:,5);

Camber_Sampel=unNomalF(Camber_Sampel_Nomal,CamberRange(1),CamberRange(2));
CamberPos_Sampel=unNomalF(CamberPos_Sampel_Nomal,CamberPosRange(1),CamberPosRange(2));
Thk_Sampel=unNomalF(Thk_Sampel_Nomal,ThickRange(1),ThickRange(2));
ThkPos_Sampel=unNomalF(ThkPos_Sampel_Nomal,ThickPosRange(1),ThickPosRange(2));
LEc_Sampel=unNomalF(LEc_Sampel_Nomal,LECircleRange(1),LECircleRange(2));

SampleData=[Camber_Sampel CamberPos_Sampel Thk_Sampel ThkPos_Sampel LEc_Sampel];
%% plot
figure(1)
xtemp=1:size(SampleData,1);

ytemp=SampleData(:,1)';
plot(xtemp,ytemp,'kx');
hold on
grid on
yline(mean(ytemp),'r');
title("Sample of Max Camber (%)");
xlabel("Sample #"); ylabel("Max Camber (%)");
legend("Sample", "MeanValue")
grid minor

figure(2)
ytemp=SampleData(:,2)';
plot(xtemp,ytemp,'kx');
hold on
grid on

yline(mean(ytemp),'r');
title("Sample of Max Camber Position (%)");
xlabel("Sample #"); ylabel("Max Camber Position (%)");
legend("Sample", "MeanValue")
grid minor

figure(3)
ytemp=SampleData(:,3)';
plot(xtemp,ytemp,'kx');
hold on
grid on

yline(mean(ytemp),'r');
title("Sample of Max Thicknees (%)");
xlabel("Sample #"); ylabel("Max Thicknees (%)");
legend("Sample", "MeanValue")
grid minor

figure(4)
ytemp=SampleData(:,4)';
plot(xtemp,ytemp,'kx');
hold on
grid on

yline(mean(ytemp),'r');
title("Sample of Max Thicknees Position (%)");
xlabel("Sample #");
ylabel("Max Thicknees Position (%)");
legend("Sample", "MeanValue")
grid minor

figure(5)
grid on
ytemp=SampleData(:,5)';
plot(xtemp,ytemp,'kx');
hold on
grid on
grid minor
yline(mean(ytemp),'r');
title("Sample of Leading Edge Circle Ratio");
xlabel("Sample #"); ylabel("Leading Edge Circle Ratio");
legend("Sample", "MeanValue");
save("sampledata.mat",'SampleData','-mat');
disp("   ")

disp("=====================================================================   ")

disp("LHS Sampling Result")
disp("    Sample Count  : "+ sampleNum)
disp("    Initial AF")
disp("        Camber      : "+ Target_camber(1)+" (%)")
disp("        Camber Pos. : "+ Target_camberPos(1)+" (%)")
disp("        Thickness   : "+ Target_Thickness(1)+" (%)")
disp("        Thick. Pos. : "+ Target_ThickPos(1)+" (%)")
disp("    Sample Range")
disp("        Camber      : "+ CamberRange(1)+" ~ "+ CamberRange(2) + " (%)")
disp("        Camber Pos. : "+ CamberPosRange(1)+" ~ "+ CamberPosRange(2) + " (%)")
disp("        Thickness   : "+ ThickRange(1)+" ~ "+ ThickRange(2) + " (%)")
disp("        Thick. Pos. : "+ ThickPosRange(1)+" ~ "+ ThickPosRange(2) + " (%)")
disp("        LE Ratio      : "+ LECircleRange(1)+" ~ "+ LECircleRange(2) + " (%)")
disp("   ")

input("Finish!, Press Enter to Mainmenu")
