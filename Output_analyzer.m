clear all; clc; close all;

%% Import data
Data = load('output.dat');

X = Data(:,1);  Y = Data(:,2);      % load x and y data
NPJ = sum(X == X(1));               % Get original matrix size from x and y
NPI = sum(Y == Y(1));
X = reshape(X,[NPJ, NPI]);          % Reshape x
Y = reshape(Y,[NPJ, NPI]);          % Reshape y

u = reshape(Data(:,3),[NPJ, NPI]);      v = reshape(Data(:,4),[NPJ, NPI]);
p = reshape(Data(:,5),[NPJ, NPI]);      T = reshape(Data(:,6),[NPJ, NPI]);
rho = reshape(Data(:,7),[NPJ, NPI]);    mu = reshape(Data(:,8),[NPJ, NPI]);

% v = Data(:,4);          p = Data(:,5);          T = Data(:,6);
% rho = Data(:,7);        mu = Data(:,8);         Gamma = Data(:,9);
% k = Data(:,10);         eps = Data(:,11);       uplus = Data(:,12); 
% yplus = Data(:,13);     yplus_u = Data(:,14);   yplus_v = Data(:,15);
% uplus_u = Data(:,16);   uplus_v = Data(:,17);

%% Load data from constraints file



%% Data analysis


figure(2)
quiver(X,Y,u,v)

Q = AddedHeat([0,0.3,0.5,0.6,0.7, 0.8, 0.9]);


function Q = AddedHeat(X_pos)

    % Initiate Q array
    Q = zeros(length(X_pos),1);
    
    % Load output data
    Data = load('output.dat');
    X = Data(:,1);  Y = Data(:,2);      % load x and y data
    NPJ = sum(X == X(1));               % Get original matrix size from x and y
    NPI = sum(Y == Y(1));
    X = reshape(X,[NPJ, NPI]);
    Y = reshape(Y,[NPJ, NPI]);
    u = reshape(Data(:,3),[NPJ, NPI]);
    T = reshape(Data(:,6),[NPJ, NPI]);
    
    % plot T
    figure(1)
    surf(X, Y, T)
    hold on
    colorbar
    view(0,90)
    
    % Get values from constraints file
    YMAX = ReadLine('constraints.dat',2);
    NPJX = ReadLine('constraints.dat',4);
    rho = ReadLine('constraints.dat',19);
    Cp = ReadLine('constraints.dat',21);
    TZERO = ReadLine('constraints.dat',16);
    % height of grid cell
    DY = YMAX / NPJX;
    
    % loop through x-positions given as input
    for i = 1 : length(X_pos)
        % find node close to given argument X_pos
        [~, I] = min(abs(X(1,:) - X_pos(i)));
        
        % loop over Y line
        for J = 1 : NPJ
            % mass flux: A * u * rho ( A = DY * 1 )
            % Q = mass flux * Cp * dT
            Q(i) = Q(i) + (DY * u(J,I) * rho) * Cp * (T(J,I) - TZERO);
        end
        
        % plot line where temperature is measured   
        line(X(:,I),Y(:,I),T(:,I),'Color','red')
        fprintf('Q_added at [x=%4.2f] = %f [W]\n',X_pos(i),Q(i))
    end
end

function out = ReadLine(filename, linenum)
    fileID = fopen(filename,'r');
    C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    fseek(fileID,0,'bof');
    out = strsplit(string(C{1}));
    out = double(out(2));
    fclose(fileID);
end
