%% Airfoil Optimizer
% Regression and Kriging Model
% kim tae jong  |   ktj4355@gmail.com   |   010 4355 1390
% Sejong University |  Propulsion Aerodynamic Lab.

clc; clear all ;close
global minCL_Const;
global maxCM_Const
options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "EnableFeasibilityMode",true,...
    "UseParallel",false);


disp("Optimizer Setting")
disp("  [1] Please input Minimum Constraint about CL at 0 AOA")
minCL_Const=input(">>> Minimum CL :  ");
disp("  [2] Please input Maximum Constraint about Cm (Inverse)")
maxCM_Const=input(">>> Maximum Cm :  ");
CVsw=0;
while(1)

    
disp("  [3] Run Test of 1 Point Cross Validation? [y/n]")
CVsw_string=string(input(">>> [y/n]      :  ",'s'));
if CVsw_string=="Y" | CVsw_string=="y" 
    CVsw=1;
    break;
elseif CVsw_string=="N" | CVsw_string=="n"
    CVsw=0;
    break;
else
    disp("Worng Input, Please input only y, Y, n, N ")
end
end
data=[];
% setting Sample
solverData=importdata("XFLR5CalculatedData.mat");
sampleData=importdata("sampledata.mat");
for ii= 1:size(solverData)
    sampelIDX=solverData(ii,1);
    data=[data;sampleData(sampelIDX,:),solverData(ii,2:end)];
end
%% Nomalize Sample Data
nomFunc=@(X, min, max) (X-min)./(max-min);
unNomFunc=@(K, min, max) K*(max-min)+min;
x1min=min(data(:,1));
x1max=max(data(:,1));
x2min=min(data(:,2));
x2max=max(data(:,2));
x3min=min(data(:,3));
x3max=max(data(:,3));
x4min=min(data(:,4));
x4max=max(data(:,4));
x5min=min(data(:,5));
x5max=max(data(:,5));

dalta=(max(data)-min(data));
nom_data=data;
nom_data(:,1)=nomFunc(nom_data(:,1),x1min,x1max);
nom_data(:,2)=nomFunc(nom_data(:,2),x2min,x2max);
nom_data(:,3)=nomFunc(nom_data(:,3),x3min,x3max);
nom_data(:,4)=nomFunc(nom_data(:,4),x4min,x4max);
nom_data(:,5)=nomFunc(nom_data(:,5),x5min,x5max);
caldata=[];


%% Kriging Method with DACE toolbox
S=nom_data(:,1:5);
LD_f=nom_data(:,6);
CL_0_f=nom_data(:,9);
inv_LD_f=-1*nom_data(:,6);
CmMax=nom_data(:,10);

theta0=0.5*ones(1,5);
lb_theta=0.001*ones(1,5);
ub_theta=1*ones(1,5);

[dmodel_LD, perf_LD] = dacefit(S,LD_f,@regpoly2,@corrgauss,theta0, lb_theta, ub_theta);
[dmodel_CL, perf_CL] = dacefit(S,CL_0_f,@regpoly2,@corrgauss,theta0, lb_theta, ub_theta);
[inv_dmodel_LD, perf_INVLD] = dacefit(S,inv_LD_f,@regpoly2,@corrgauss,theta0, lb_theta, ub_theta);
[dmodel_CmMax, perf_CmMax] = dacefit(S,CmMax,@regpoly2,@corrgauss,theta0, lb_theta, ub_theta);



global model_krig_LD;
global model_krig_CL;
global model_krig_CmMax;

global model_krig_invLD;

model_krig_LD=@(X)predictor(X,dmodel_LD);
model_krig_CL=@(X)predictor(X,dmodel_CL);
model_krig_invLD=@(X)predictor(X,inv_dmodel_LD);
model_krig_CmMax=@(X)predictor(X,dmodel_CmMax);

%% Cross Validation Kriging model
LD_CV_err=0;
CL_CV_err=0;
CM_CV_err=0;

if CVsw==1
    LD_err=[];
    CL_err=[];
    CM_err=[];
    for ii = 1:size(sampleData,1)
        clc;
        disp("1 Point Cross Validation Test Run...")
        disp("Test Point : "+ii+"   Progress... "+100*ii/(size(sampleData,1))+"%");
        sample_CV=S;
        LD_f_CV=LD_f;
        CL_0_f_CV=CL_0_f;
        CmMax_CV=CmMax;

        sample_CV(ii,:)=[];
        LD_f_CV(ii,:)=[];
        CL_0_f_CV(ii,:)=[];
        CmMax_CV(ii,:)=[];

        [CV_dmodel_LD, perf_LD] = dacefit(sample_CV,LD_f_CV,@regpoly2,@corrgauss,theta0, lb_theta, ub_theta);
        [CV_dmodel_CL, perf_CL] = dacefit(sample_CV,CL_0_f_CV,@regpoly2,@corrgauss,theta0, lb_theta, ub_theta);
        [CV_dmodel_CmMax, perf_CmMax] = dacefit(sample_CV,CmMax_CV,@regpoly2,@corrgauss,theta0, lb_theta, ub_theta);
        CV_model_krig_LD=@(X)predictor(X,dmodel_LD);
        CV_model_krig_CL=@(X)predictor(X,dmodel_CL);
        CV_model_krig_CmMax=@(X)predictor(X,dmodel_CmMax);
        CV_X=S(ii,:);
        LD_Real=LD_f(ii);
        LD_CV=CV_model_krig_LD(CV_X);
        LD_err=[LD_err;(LD_Real-LD_CV).^2];

        CL_Real=CL_0_f(ii);
        CL_CV=CV_model_krig_CL(CV_X);
        CL_err=[CL_err;(CL_Real-CL_CV).^2];

        CM_Real=CmMax(ii);
        CM_CV=CV_model_krig_CmMax(CV_X);
        CM_err=[CM_err;(CM_Real-CM_CV).^2];

    end
    LD_Real_minmax=[min(LD_f) max(LD_f)];
    CL_Real_minmax=[min(CL_0_f) max(CL_0_f)];
    CM_Real_minmax=[min(CmMax) max(CmMax)];

    LD_CV_err=(sum(LD_err)/(length(LD_err)^2))/(LD_Real_minmax(2)-LD_Real_minmax(1));
    CL_CV_err=(sum(CL_err)/(length(CL_err)^2))/(CL_Real_minmax(2)-CL_Real_minmax(1));
    CM_CV_err=(sum(CM_err)/(length(CM_err)^2))/(CM_Real_minmax(2)-CM_Real_minmax(1));

    disp("1 Point Crossvalidation Error")
    disp("     LD model   : "+ LD_CV_err)
    disp("     CL model   : "+ CL_CV_err)
    disp("     Cm model   : "+ CM_CV_err)
end

%% Test of kriging model
calkrigdata=[];
for kk=1:size(nom_data)
    x1=nom_data(kk,1);
    x2=nom_data(kk,2);
    x3=nom_data(kk,3);
    x4=nom_data(kk,4);
    x5=nom_data(kk,5);
    X=[x1 x2 x3 x4 x5];
    cal_krig_LD=model_krig_LD(X);
    cal_krig_CL=model_krig_CL(X);
    cal_krig_CmMax=model_krig_CmMax(X);
    calkrigdata=[calkrigdata;cal_krig_LD, cal_krig_CL,cal_krig_CmMax];
end
sampleSize=1:size(calkrigdata(:,1),1);
realData=nom_data(:,[6,9,10]);
figure(2)
sgtitle("Comparison between Real and kriging Model")
hold on
subplot(3,1,1)
hold on
plot(sampleSize, realData(:,1),"kx")
plot(sampleSize, calkrigdata(:,1),"ko")
xlabel("# of Sample")
ylabel("L/D")
legend("DOE Data", "Kriging Data")
subplot(3,1,2)
hold on
plot(sampleSize, realData(:,2),"kx")
plot(sampleSize, calkrigdata(:,2),"ko")
xlabel("# of Sample")
ylabel("Cl at AOA 0 (deg)")
legend("DOE Data", "Kriging Data")
subplot(3,1,3)
hold on
plot(sampleSize, realData(:,3),"kx")
plot(sampleSize, calkrigdata(:,3),"ko")
xlabel("# of Sample")
ylabel("inverse Max CM")
legend("DOE Data", "Kriging Data")
%% Optimizing with fmincon
x0=[0.5 0.5 0.5 0.5 0.5];

[opt_X_LD_nom_krig, opt_LD_krig]=fmincon(model_krig_invLD,x0,[],[],[],[],[0 0 0 0 0],[1 1 1 1 1],'constraintB',options);
%[opt_X_CL_nom_krig, opt_CL_krig]=fmincon(model_krig_CL,x0,[],[],[],[],[0 0 0 0 0],[1 1 1 1 1],'constraintB');
opt_LD_krig=-opt_LD_krig;
CL_at_optLD_krig=model_krig_CL(opt_X_LD_nom_krig);
CmMax_at_optLD_krig=model_krig_CmMax(opt_X_LD_nom_krig);

disp("LD MAX optimization (kriging model)")
disp("    Maximum L/D   : "+ opt_LD_krig)
disp("    CL at point   : "+ CL_at_optLD_krig)
disp("    CmMax         : "+ CmMax_at_optLD_krig)
disp("    x1(Nomalized) : "+ opt_X_LD_nom_krig(1))
disp("    x2(Nomalized) : "+ opt_X_LD_nom_krig(2))
disp("    x3(Nomalized) : "+ opt_X_LD_nom_krig(3))
disp("    x4(Nomalized) : "+ opt_X_LD_nom_krig(4))
disp("    x5(Nomalized) : "+ opt_X_LD_nom_krig(5))

disp("    x1(Real)      : "+ unNomFunc(opt_X_LD_nom_krig(1),x1min,x1max))
disp("    x2(Real)      : "+ unNomFunc(opt_X_LD_nom_krig(2),x2min,x2max))
disp("    x3(Real)      : "+ unNomFunc(opt_X_LD_nom_krig(3),x3min,x3max))
disp("    x4(Real)      : "+ unNomFunc(opt_X_LD_nom_krig(4),x4min,x4max))
disp("    x5(Real)      : "+ unNomFunc(opt_X_LD_nom_krig(5),x5min,x5max))

disp(" ")



X1R_krg= unNomFunc(opt_X_LD_nom_krig(1),x1min,x1max);
X2R_krg= unNomFunc(opt_X_LD_nom_krig(2),x2min,x2max);
X3R_krg= unNomFunc(opt_X_LD_nom_krig(3),x3min,x3max);
X4R_krg= unNomFunc(opt_X_LD_nom_krig(4),x4min,x4max);
X5R_krg= unNomFunc(opt_X_LD_nom_krig(5),x5min,x5max);
opt_out_krg= [X1R_krg X2R_krg X3R_krg X4R_krg X5R_krg];

%% Save Airfoil and Optimizing Data
status = rmdir('OUTPUT DATA','s');
mkdir 'OUTPUT DATA'\;
delete("OUTPUT DATA\*.*")
delete("OUTPUT DATA\*.*")
name="output_krig_airfoil";
outAirfoil=opt_out_krg;
%insert base Coord data
disp("input Base Airfoil")
[fileFullname,Wd] = uigetfile({'*.dat';'*.*'});
filename=split(fileFullname,".");
for ii=1
    inpFileName= "output.inp";
    fid = fopen(inpFileName,'w');
    modifingDATA=outAirfoil(ii,:);

    if (fid<=0)
        error([mfilename ':io'],'Unable to create xfoil.inp file');
        continue;
    end
    fprintf(fid,'load %s\n',fileFullname);
    fprintf(fid,'\nppar\n');
    fprintf(fid,'N\n200\n');
    fprintf(fid,'\n\ngdes\n');
    %camber and Thickness
    fprintf(fid,'tset %f %f\n',modifingDATA(3)/100.0,modifingDATA(1)/100.0);
    %camber and Thickreness xCordinate
    fprintf(fid,'high %f %f\n',modifingDATA(4)/100.0,modifingDATA(2)/100.0);
    % Leading Edge
    fprintf(fid,'lera %f\n\n',modifingDATA(5));
    fprintf(fid,'eXec\ngset\n\n');

    fprintf(fid,'\nNAME %s\n',name(ii));
    fprintf(fid,'\nSAVE %s\n',"OUTPUT DATA\"+name(ii)+".dat");
    fprintf(fid,'quit\n');

    fclose(fid);
    % input Xfoil.exe
    wd = fileparts(which(mfilename)); % working directory, where xfoil.exe needs to be
    cmd = sprintf('cd %s && xfoil.exe < %s > xfoil.out',wd,inpFileName);
    [status,result] = system(cmd);
    if (status~=0),
        disp(result);
        error([mfilename ':system'],'Xfoil execution failed! %s',cmd);
    end;
end

%% output FIle data
ExpFID = fopen('Export Data.txt','w');


fprintf(ExpFID,"Kriging Model Report \n====================================================\n");
fprintf(ExpFID,"\tModel :\tKriging (DACE) Model\n");
fprintf(ExpFID,"\n\t1 point Cross Validation Error (Nomalized)\n");
fprintf(ExpFID,"\t\tL/D Model :\t%g\n",LD_CV_err);
fprintf(ExpFID,"\t\tCL Model  :\t%g\n",CL_CV_err);
fprintf(ExpFID,"\t\tCm Model  :\t%g\n\n",CM_CV_err);
fprintf(ExpFID,"\tDesign Point Value\n");
fprintf(ExpFID,"\t\tMaximum L/D  :\t%f\n",opt_LD_krig);
fprintf(ExpFID,"\t\tCL at 0 AOA  :\t%f\n",CL_at_optLD_krig);
fprintf(ExpFID,"\t\tMaximum Cm   :\t%f\n\n",CmMax_at_optLD_krig);
fprintf(ExpFID,"\tOptimized Design Value\n");
fprintf(ExpFID,"\t\tMax Camber \t at Position \t Max Thickness \t at Position \t LE Circle Ratio\n");
fprintf(ExpFID,"\t----------------------------------------------------------------------------------\n");
fprintf(ExpFID,"\t\t%f \t %f \t\t %f \t\t %f \t\t %f\n",opt_out_krg);



fprintf(ExpFID,"\n\n\nOptimizing RAW Data \n====================================================\n");

fprintf(ExpFID,"\n\nLD MAX optimization (Kriging Model)");
fprintf(ExpFID,"\n\tMaximum L/D   : "+ opt_LD_krig);
fprintf(ExpFID,"\n\tCL at point   : "+ CL_at_optLD_krig);
fprintf(ExpFID,"\n\tCmMax         : "+ CmMax_at_optLD_krig);
fprintf(ExpFID,"\n\tx1(Nomalized) : "+ opt_X_LD_nom_krig(1));
fprintf(ExpFID,"\n\tx2(Nomalized) : "+ opt_X_LD_nom_krig(2));
fprintf(ExpFID,"\n\tx3(Nomalized) : "+ opt_X_LD_nom_krig(3));
fprintf(ExpFID,"\n\tx4(Nomalized) : "+ opt_X_LD_nom_krig(4));
fprintf(ExpFID,"\n\tx5(Nomalized) : "+ opt_X_LD_nom_krig(5));

fprintf(ExpFID,"\n\tx1(Real)      : "+ unNomFunc(opt_X_LD_nom_krig(1),x1min,x1max));
fprintf(ExpFID,"\n\tx2(Real)      : "+ unNomFunc(opt_X_LD_nom_krig(2),x2min,x2max));
fprintf(ExpFID,"\n\tx3(Real)      : "+ unNomFunc(opt_X_LD_nom_krig(3),x3min,x3max));
fprintf(ExpFID,"\n\tx4(Real)      : "+ unNomFunc(opt_X_LD_nom_krig(4),x4min,x4max));
fprintf(ExpFID,"\n\tx5(Real)      : "+ unNomFunc(opt_X_LD_nom_krig(5),x5min,x5max));
fprintf(ExpFID,"\n====================================================\n\n");

fclose(ExpFID);
winopen("Export Data.txt")

input("Finish!, Press Enter to Mainmenu");

