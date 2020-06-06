%% Siemens:Internship %%
% Author: Jason Le
% Date created: June 5, 2020
% Date modified: June 6, 2020
% Purpose: To analyze NACA63(3)-618 airfoil at different Reynolds numbers
% between -10 to 20 AoA
% Inputs: Certain Reynolds numbers and data from Xfoil.exe
% Outputs: Plots of Aerodynamic coefficients,CP Distributions, and Boundary
% Layer properties at the trailing edge
% Assumptions: No friction except in braking section, cart is point mass

%% House Keeping
clear
close all
clc

%% Pre-Allocated Variables
RE = ['3000000 ';'10000000';'15000000']; %Reynolds numbers
airfoildata = 'Airfoildata.txt'; %Name of file that matlab reads from xfoil
AoACP = ['0 ';'4 ';'8 ';'12']; %AoA for CP Distributions
AoABL = ['-10';'-9 ';'-8 ';'-7 ';'-6 ';'-5 ';'-4 ';'-3 ';'-2 ';'-1 ';'0  '...
         ;'1  ';'2  ';'3  ';'4  ';'5  ';'6  ';'7  ';'8  ';'9  ';'10 '...
         ;'11 ';'12 ';'13 ';'14 ';'15 ';'16 ';'17 ';'18 ';'19 ';'20 ']; %string for parsing xfoil data
AoABLnum = -10:20; %double data type AoA
%Pre-Allocated Variables for Boundary Layer Properties
DStar = zeros(size(AoABL,1),size(RE,1));
Theta = zeros(size(AoABL,1),size(RE,1));
H = zeros(size(AoABL,1),size(RE,1));

%% Edge Case 
% Delete files if they exist
if (exist(airfoildata,'file'))
    delete(airfoildata);
end

%% Getting Aerodynamic Coefficients
%Create notepad to read into xfoil for Reynolds number 3e6 
fid = fopen('xfoil_input.txt','w');
fprintf(fid,'load NACA633618.txt \n');
fprintf(fid,'NACA633618 \n');
fprintf(fid,'pane \n');
fprintf(fid,'oper \n');
fprintf(fid,'iter 250 \n');
fprintf(fid,['visc ' RE(1,:) '\n']);
fprintf(fid,'seqp \n');
fprintf(fid,'pacc \n');
fprintf(fid,[airfoildata '\n']);
fprintf(fid,'\n');
fprintf(fid,'aseq -10 20 1 \n');
fclose(fid);
%Load notepad commands into xfoil
cmd = 'xfoil.exe < xfoil_input.txt';
[~,~] = system(cmd);
%Parsing Data from xfoil
fidData = fopen(airfoildata);
Dataraw = textscan(fidData,'%f %f %f %f %f %f %f','HeaderLines',12,'CollectOutput',1,'Delimiter','');
fclose(fidData);
delete(airfoildata);
AoA1 = Dataraw{1}(:,1);
CL1 = Dataraw{1}(:,2);
CD1 = Dataraw{1}(:,3);
CM1 = Dataraw{1}(:,5);
CLCD1 = CL1./CD1;

%Create notepad to read into xfoil for Reynolds number 10e6 
fid = fopen('xfoil_input.txt','w');
fprintf(fid,'load NACA633618.txt \n');
fprintf(fid,'NACA633618 \n');
fprintf(fid,'pane \n');
fprintf(fid,'oper \n');
fprintf(fid,'iter 250 \n');
fprintf(fid,['visc ' RE(2,:) '\n']);
fprintf(fid,'seqp \n');
fprintf(fid,'pacc \n');
fprintf(fid,[airfoildata '\n']);
fprintf(fid,'\n');
fprintf(fid,'aseq -10 20 1 \n');
fclose(fid);
%Load notepad commands into xfoil
cmd = 'xfoil.exe < xfoil_input.txt';
[~,~] = system(cmd);
%Parsing Data from xfoil
fidData = fopen(airfoildata);
Dataraw = textscan(fidData,'%f %f %f %f %f %f %f','HeaderLines',12,'CollectOutput',1,'Delimiter','');
fclose(fidData);
delete(airfoildata);
AoA2 = Dataraw{1}(:,1);
CL2 = Dataraw{1}(:,2);
CD2 = Dataraw{1}(:,3);
CM2 = Dataraw{1}(:,5);
CLCD2 = CL2./CD2;

%Create notepad to read into xfoil for Reynolds number 15e6 
fid = fopen('xfoil_input.txt','w');
fprintf(fid,'load NACA633618.txt \n');
fprintf(fid,'NACA633618 \n');
fprintf(fid,'pane \n');
fprintf(fid,'oper \n');
fprintf(fid,'iter 250 \n');
fprintf(fid,['visc ' RE(3,:) '\n']);
fprintf(fid,'seqp \n');
fprintf(fid,'pacc \n');
fprintf(fid,[airfoildata '\n']);
fprintf(fid,'\n');
fprintf(fid,'aseq -10 20 1 \n');
fclose(fid);
%Load notepad commands into xfoil
cmd = 'xfoil.exe < xfoil_input.txt';
[~,~] = system(cmd);
%Parsing Data from xfoil
fidData = fopen(airfoildata);
Dataraw = textscan(fidData,'%f %f %f %f %f %f %f','HeaderLines',12,'CollectOutput',1,'Delimiter','');
fclose(fidData);
delete(airfoildata);
AoA3 = Dataraw{1}(:,1);
CL3 = Dataraw{1}(:,2);
CD3 = Dataraw{1}(:,3);
CM3 = Dataraw{1}(:,5);
CLCD3 = CL3./CD3;

%% Getting CP Distribution
for i = 1:size(AoACP,1)
%This for loop goes from 0, 4, 8, and 12 AoA for Reynolds number 3e6
    %Create notepad to read into xfoil for Reynolds number 3e6
    fid = fopen('xfoil_input.txt','w');
    fprintf(fid,'load NACA633618.txt \n');
    fprintf(fid,'NACA633618 \n');
    fprintf(fid,'pane \n');
    fprintf(fid,'oper \n');
    fprintf(fid,'iter 250 \n');
    fprintf(fid,['visc ' RE(1,:) '\n']);
    fprintf(fid,['alfa ' AoACP(i,:) '\n']);
    fprintf(fid,['CPWR ' airfoildata '\n']);
    fclose(fid);
    %Load notepad commands into xfoil
    cmd = 'xfoil.exe < xfoil_input.txt';
    [~,~] = system(cmd);
    %Parsing Data from xfoil
    fidData = fopen(airfoildata);
    Dataraw = textscan(fidData,'%f %f %f','HeaderLines',3,'CollectOutput',1,'Delimiter','');
    fclose(fidData);
    delete(airfoildata);
    if i == 1
    %if statement is for the first iteration pre-allocation
        X1 = zeros(numel(Dataraw{1}(:,1)),size(AoACP,1));
        CP1 = zeros(numel(Dataraw{1}(:,1)),size(AoACP,1));
    end
    X1(:,i) = Dataraw{1}(:,1);
    CP1(:,i) = Dataraw{1}(:,3);
end

for k = 1:size(AoACP,1)
%This for loop goes from 0, 4, 8, and 12 AoA for Reynolds number 10e6
    %Create notepad to read into xfoil for Reynolds number 10e6
    fid = fopen('xfoil_input.txt','w');
    fprintf(fid,'load NACA633618.txt \n');
    fprintf(fid,'NACA633618 \n');
    fprintf(fid,'pane \n');
    fprintf(fid,'oper \n');
    fprintf(fid,'iter 250 \n');
    fprintf(fid,['visc ' RE(2,:) '\n']);
    fprintf(fid,['alfa ' AoACP(k,:) '\n']);
    fprintf(fid,['CPWR ' airfoildata '\n']);
    fclose(fid);
    %Load notepad commands into xfoil
    cmd = 'xfoil.exe < xfoil_input.txt';
    [~,~] = system(cmd);
    %Parsing Data from xfoil
    fidData = fopen(airfoildata);
    Dataraw = textscan(fidData,'%f %f %f','HeaderLines',3,'CollectOutput',1,'Delimiter','');
    fclose(fidData);
    delete(airfoildata);
    if k == 1
    %if statement is for the first iteration pre-allocation
        X2 = zeros(numel(Dataraw{1}(:,1)),size(AoACP,1));
        CP2 = zeros(numel(Dataraw{1}(:,1)),size(AoACP,1));
    end
    X2(:,k) = Dataraw{1}(:,1);
    CP2(:,k) = Dataraw{1}(:,3);
end

for l = 1:size(AoACP,1)
%This for loop goes from 0, 4, 8, and 12 AoA for Reynolds number 15e6
    %Create notepad to read into xfoil for Reynolds number 15e6
    fid = fopen('xfoil_input.txt','w');
    fprintf(fid,'load NACA633618.txt \n');
    fprintf(fid,'NACA633618 \n');
    fprintf(fid,'pane \n');
    fprintf(fid,'oper \n');
    fprintf(fid,'iter 250 \n');
    fprintf(fid,['visc ' RE(3,:) '\n']);
    fprintf(fid,['alfa ' AoACP(l,:) '\n']);
    fprintf(fid,['CPWR ' airfoildata '\n']);
    fclose(fid);
    %Load notepad commands into xfoil
    cmd = 'xfoil.exe < xfoil_input.txt';
    [~,~] = system(cmd);
    fidData = fopen(airfoildata);
    %Parsing Data from xfoil
    Dataraw = textscan(fidData,'%f %f %f','HeaderLines',3,'CollectOutput',1,'Delimiter','');
    fclose(fidData);
    delete(airfoildata);
    if l == 1
    %if statement is for the first iteration pre-allocation
        X3 = zeros(numel(Dataraw{1}(:,1)),size(AoACP,1));
        CP3 = zeros(numel(Dataraw{1}(:,1)),size(AoACP,1));
    end
    X3(:,l) = Dataraw{1}(:,1);
    CP3(:,l) = Dataraw{1}(:,3);
end

%% Getting Boundary Layer Properties
for b = 1:size(AoABL,1)
%This for loop goes from -10 to 20 AoA for Reynolds number 3e6
    %Create notepad to read into xfoil for Reynolds number 3e6
    fid = fopen('xfoil_input.txt','w');
    fprintf(fid,'load NACA633618.txt \n');
    fprintf(fid,'NACA633618 \n');
    fprintf(fid,'pane \n');
    fprintf(fid,'oper \n');
    fprintf(fid,'iter 250 \n');
    fprintf(fid,['visc ' RE(1,:) '\n']);
    fprintf(fid,['alfa ' AoABL(b,:) '\n']);
    fprintf(fid,['DUMP ' airfoildata '\n']);
    fclose(fid);
    %Load notepad commands into xfoil
    cmd = 'xfoil.exe < xfoil_input.txt';
    [~,~] = system(cmd);
    %Parsing Data from xfoil
    fidData = fopen(airfoildata);
    Dataraw = textscan(fidData,'%f %f %f %f %f %f %f %f','HeaderLines',1,'CollectOutput',1,'Delimiter','');
    fclose(fidData);
    delete(airfoildata);
    XBL = Dataraw{1}(:,2);
    idx = find(XBL ==  1);
    DStar(b,1) = Dataraw{1}(idx(2),5);
    Theta(b,1) = Dataraw{1}(idx(2),6);
    H(b,1) = Dataraw{1}(idx(2),8);
end

for n = 1:size(AoABL,1)
%This for loop goes from -10 to 20 AoA for Reynolds number 10e6
    %Create notepad to read into xfoil for Reynolds number 10e6
    fid = fopen('xfoil_input.txt','w');
    fprintf(fid,'load NACA633618.txt \n');
    fprintf(fid,'NACA633618 \n');
    fprintf(fid,'pane \n');
    fprintf(fid,'oper \n');
    fprintf(fid,'iter 250 \n');
    fprintf(fid,['visc ' RE(2,:) '\n']);
    fprintf(fid,['alfa ' AoABL(n,:) '\n']);
    fprintf(fid,['DUMP ' airfoildata '\n']);
    fclose(fid);
    %Load notepad commands into xfoil
    cmd = 'xfoil.exe < xfoil_input.txt';
    [~,~] = system(cmd);
    %Parsing Data from xfoil
    fidData = fopen(airfoildata);
    Dataraw = textscan(fidData,'%f %f %f %f %f %f %f %f','HeaderLines',1,'CollectOutput',1,'Delimiter','');
    fclose(fidData);
    delete(airfoildata);
    XBL = Dataraw{1}(:,2);
    idx = find(XBL ==  1);
    DStar(n,2) = Dataraw{1}(idx(2),5);
    Theta(n,2) = Dataraw{1}(idx(2),6);
    H(n,2) = Dataraw{1}(idx(2),8);
end

for m = 1:size(AoABL,1)
%This for loop goes from -10 to 20 AoA for Reynolds number 15e6
    %Create notepad to read into xfoil for Reynolds number 15e6
    fid = fopen('xfoil_input.txt','w');
    fprintf(fid,'load NACA633618.txt \n');
    fprintf(fid,'NACA633618 \n');
    fprintf(fid,'pane \n');
    fprintf(fid,'oper \n');
    fprintf(fid,'iter 250 \n');
    fprintf(fid,['visc ' RE(3,:) '\n']);
    fprintf(fid,['alfa ' AoABL(m,:) '\n']);
    fprintf(fid,['DUMP ' airfoildata '\n']);
    fclose(fid);
    %Load notepad commands into xfoil
    cmd = 'xfoil.exe < xfoil_input.txt';
    [~,~] = system(cmd);
    %Parsing Data from xfoil
    fidData = fopen(airfoildata);
    Dataraw = textscan(fidData,'%f %f %f %f %f %f %f %f','HeaderLines',1,'CollectOutput',1,'Delimiter','');
    fclose(fidData);
    delete(airfoildata);
    XBL = Dataraw{1}(:,2);
    idx = find(XBL ==  1);
    DStar(m,3) = Dataraw{1}(idx(2),5);
    Theta(m,3) = Dataraw{1}(idx(2),6);
    H(m,3) = Dataraw{1}(idx(2),8);
end

%% Plots
%Plot for CL vs AoA
figure()
hold on
plot(AoA1,CL1,'g')
plot(AoA2,CL2,'b')
plot(AoA3,CL3,'r')
title('CL vs AoA','fontsize',15,'interpreter','latex')
xlabel('Angle of Attack [deg]','fontsize',15,'interpreter','latex')
ylabel('CL','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000','Location','northwest')
hold off

%Plot for CD vs AoA
figure()
hold on
plot(AoA1,CD1,'g')
plot(AoA2,CD2,'b')
plot(AoA3,CD3,'r')
title('CD vs AoA','fontsize',15,'interpreter','latex')
xlabel('Angle of Attack [deg]','fontsize',15,'interpreter','latex')
ylabel('CD','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000','Location','northwest')
hold off

%Plot for CM vs AoA
figure()
hold on
plot(AoA1,CM1,'g')
plot(AoA2,CM2,'b')
plot(AoA3,CM3,'r')
title('CM vs AoA','fontsize',15,'interpreter','latex')
xlabel('Angle of Attack [deg]','fontsize',15,'interpreter','latex')
ylabel('CM','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000','Location','northwest')
hold off

%Plot for CL/CD vs AoA
figure()
hold on
plot(AoA1,CLCD1,'g')
plot(AoA2,CLCD2,'b')
plot(AoA3,CLCD3,'r')
title('CL/CD vs AoA','fontsize',15,'interpreter','latex')
xlabel('Angle of Attack [deg]','fontsize',15,'interpreter','latex')
ylabel('CL/CD','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000','Location','northwest')
hold off

%Plot for CP Distribution of 0 deg AoA
figure()
hold on
plot(X1(:,1),CP1(:,1),'g')
plot(X2(:,1),CP2(:,1),'b')
plot(X3(:,1),CP3(:,1),'r')
title('CP Distribution at 0 deg','fontsize',15,'interpreter','latex')
xlabel('Percent Cord x/c','fontsize',15,'interpreter','latex')
ylabel('CP','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000')
set(gca,'Ydir','reverse')

%Plot for CP Distribution of 4 deg AoA
figure()
hold on
plot(X1(:,2),CP1(:,2),'g')
plot(X2(:,2),CP2(:,2),'b')
plot(X3(:,2),CP3(:,2),'r')
title('CP Distribution at 4 deg','fontsize',15,'interpreter','latex')
xlabel('Percent Cord x/c','fontsize',15,'interpreter','latex')
ylabel('CP','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000')
set(gca,'Ydir','reverse')

%Plot for CP Distribution of 8 deg AoA
figure()
hold on
plot(X1(:,3),CP1(:,3),'g')
plot(X2(:,3),CP2(:,3),'b')
plot(X3(:,3),CP3(:,3),'r')
title('CP Distrubition at 8 deg','fontsize',15,'interpreter','latex')
xlabel('Percent Cord x/c','fontsize',15,'interpreter','latex')
ylabel('CP','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000')
set(gca,'Ydir','reverse')

%Plot for CP Distribution of 12 deg AoA
figure()
hold on
plot(X1(:,4),CP1(:,4),'g')
plot(X2(:,4),CP2(:,4),'b')
plot(X3(:,4),CP3(:,4),'r')
title('CP Distribution at 12 deg','fontsize',15,'interpreter','latex')
xlabel('Percent Cord x/c','fontsize',15,'interpreter','latex')
ylabel('CP','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000')
set(gca,'Ydir','reverse')

%Plot for Displacement Thickness vs AoA
figure()
hold on
plot(AoABLnum,DStar(:,1),'g')
plot(AoABLnum,DStar(:,2),'b')
plot(AoABLnum,DStar(:,3),'r')
title('Displacement Thickness vs AoA','fontsize',15,'interpreter','latex')
xlabel('Angle of Attack [deg]','fontsize',15,'interpreter','latex')
ylabel('Displacement Thickness [delta/c]','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000')
hold off

%Plot for Momentum Thickness vs AoA
figure()
hold on
plot(AoABLnum,Theta(:,1),'g')
plot(AoABLnum,Theta(:,2),'b')
plot(AoABLnum,Theta(:,3),'r')
title('Momentum Thickness vs AoA [theta/c]','fontsize',15,'interpreter','latex')
xlabel('Angle of Attack [deg]','fontsize',15,'interpreter','latex')
ylabel('Momentum Thickness','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000')
hold off

%Plot for Shape Factor vs AoA
figure()
hold on
plot(AoABLnum,H(:,1),'g')
plot(AoABLnum,H(:,2),'b')
plot(AoABLnum,H(:,3),'r')
title('Shape Factor vs AoA','fontsize',15,'interpreter','latex')
xlabel('Angle of Attack [deg]','fontsize',15,'interpreter','latex')
ylabel('Shape Factor','fontsize',15,'interpreter','latex')
legend('RE = 3000000','RE = 10000000','RE = 15000000')
hold off