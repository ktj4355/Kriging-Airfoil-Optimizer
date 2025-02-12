%% Airfoil Optimizer
% Regression and Kriging Model
% kim tae jong  |   ktj4355@gmail.com   |   010 4355 1390
% Sejong University |  Propulsion Aerodynamic Lab.
clc; clear; close all;
extcond=0

while(extcond==0)
clc;
    disp("Airfoil Optimizer (by Taejong Kim)")
    disp("Based on Regression and Kriging Optimizer")
    disp("Read Manual Before Using Tool")
    disp(" ")

    disp("Please Seclect Mode")
    disp("  1) Sampling by Latin Hypercube Sampling ")
    disp("  2) Airfoil Design by Sample Point")
    disp("  3) Calculate Data about XFLR5 OUTPUT")
    disp("  4) Optimizing and Export Airfoil ")
    disp("  5) Exit")

    inpNum=input(">>>   ","s");
    inpNum=int16(str2double(inpNum));
    if ~isempty(inpNum) & isinteger(inpNum)
        switch inpNum
            case 1
                run("LHSample.m");
                extcond=0;

            case 2
                run("DesiginAirfoil.m");
                extcond=0;
            case 3
                 run("findObjectValue.m");
                extcond=0;
            case 4
                run("Optimizer.m");
                extcond=0;
            case 5
                extcond=1;
            otherwise
                continue;

        end
    end

end

