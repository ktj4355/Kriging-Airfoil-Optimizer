%% Airfoil Optimizer
% Regression and Kriging Model
% kim tae jong  |   ktj4355@gmail.com   |   010 4355 1390
% Sejong University |  Propulsion Aerodynamic Lab.

clc; clear; close all;
caseoutput=[];

A=dir("XFLR5 DATA");
B=struct2cell(A)';
C=B(:,1);
C(1:2,:)=[];

for Ci=1:size(C);
Dtmp=split(C(Ci,1),"_");
Ddtmp=split(Dtmp(1),"case");
C{Ci,2}=str2num(Ddtmp{2});
end
C=sortrows(C,2) 

for caseII=1:size(C,1)
  
    fid=fopen("XFLR5 DATA\"+C{caseII,1});
    if fid<0;break;end
    clc;
    disp("XFLR5 Data Sort and Calcuator")
    disp("If you have Error, Please Check that you have been Save XFLR5 data on ""XFLR5 DATA"" Folder")
    disp("Run..... "+caseII+"Case ")

    Rawdata = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);

    Rawdata=Rawdata{1};
    Rawdata(1:11)=[];
    Rawdata=string(Rawdata);
    data=[];
    %data ;   alpha | CL | CD | CDp | Cm |Top Xtr | Bot Xtr | Cpmin | Chinge | XCp

    for ii=1:size(Rawdata,1)

        datastring=Rawdata{ii};
        dataLine=str2double(split(datastring));
        if isnan(dataLine)
            break;
        end
        data=[data;dataLine'];

    end
    aoaRange=[0 4];
    if or(isempty(data),size(data,1)<5);continue;end;
    intp_Space=transpose(-2:0.1:10);
    data_Cl=interp1(data(:,1),data(:,2),intp_Space);
    data_Cd=interp1(data(:,1),data(:,3),intp_Space);
    data_LD=data_Cl./data_Cd;
    data_Cm=interp1(data(:,1),data(:,5),intp_Space);
    MaxABSCm=max(abs(data_Cm));
    MODDATA=[intp_Space,data_Cl,data_Cd,data_LD];
    ind_aoa_fst=find(intp_Space==aoaRange(1));
        ind_aoa_end=find(intp_Space==aoaRange(2));
    ldmaxdata=data_LD(ind_aoa_fst:ind_aoa_end);
    [maxLD maxIND]=max(ldmaxdata);
    
    cl0=[0 data_Cl(find(intp_Space==0))];
    if isnan(cl0(2)), continue;end
    caseoutput=[caseoutput;caseII,maxLD, intp_Space(maxIND+ind_aoa_fst-1),cl0,MaxABSCm];
end
%% test
meanOutput=mean(caseoutput);
stdOutput=mad(caseoutput);
range=meanOutput(2)+2.5.*(stdOutput(2)).*[-1,+1]

filtered_Output=[];
for jj=1:size(caseoutput,1)
    
    if (caseoutput(jj,2)<range(1))|(caseoutput(jj,2)>range(2))
        continue;
    end
    filtered_Output=[filtered_Output;caseoutput(jj,:)];
end
close
figure(1)
hold on
plot(1:size(caseoutput,1),caseoutput(:,2),"ko")
plot(1:size(caseoutput,1),meanOutput(:,2),"kx")
plot(1:size(caseoutput,1),range(1),"ro")
plot(1:size(caseoutput,1),range(2),"ro")


%% 

save("XFLR5CalculatedData.mat",'filtered_Output','-mat');
input("Finish!, Press Enter to Mainmenu")

