clear all; clc; close all;

%% Import data
Data = load('output.dat');

X = Data(:,1);  Y = Data(:,2);      % load x and y data
NPI = sum(X == X(1));               % Get original matrix size from x and y
NPJ = sum(Y == Y(1));
X = reshape(X,[NPI, NPJ]);          % Reshape x
Y = reshape(Y,[NPI, NPJ]);          % Reshape y

u = reshape(Data(:,3),[NPI, NPJ]);      v = reshape(Data(:,4),[NPI, NPJ]);
p = reshape(Data(:,5),[NPI, NPJ]);      T = reshape(Data(:,6),[NPI, NPJ]);
rho = reshape(Data(:,7),[NPI, NPJ]);    mu = reshape(Data(:,8),[NPI, NPJ]);

% v = Data(:,4);          p = Data(:,5);          T = Data(:,6);
% rho = Data(:,7);        mu = Data(:,8);         Gamma = Data(:,9);
% k = Data(:,10);         eps = Data(:,11);       uplus = Data(:,12); 
% yplus = Data(:,13);     yplus_u = Data(:,14);   yplus_v = Data(:,15);
% uplus_u = Data(:,16);   uplus_v = Data(:,17);

%% Load data from constraints file



%% Data analysis

surf(X, Y, T)
colorbar
%view(0,90)

ReadLine('constraints.dat',1)

function out = ReadLine(filename, linenum);
    fileID = fopen(filename,'r');
    C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    fseek(fileID,0,'bof');
    out = strsplit(string(C{1}))
    fclose(fileID)
end
