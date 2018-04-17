clear all; clc; close all;

%% Import data
Data = load('output.dat');

X = Data(:,1);  Y = Data(:,2);      % load x and y data
NPJ = sum(X == X(1));               % Get original matrix size from x and y
NPI = sum(Y == Y(1));
X = reshape(X,[NPJ, NPI]);          % Reshape x
Y = reshape(Y,[NPJ, NPI]);          % Reshape y

%
%vq = griddata(x,y,v,xgrid,ygrid);
%uq = griddata(x,y,u,xgrid,ygrid);

YMAX = ReadLine('constraints.dat',2);
NPJX = ReadLine('constraints.dat',4);
DY = YMAX / NPJX;

XMAX = ReadLine('constraints.dat',1);
NPIX = ReadLine('constraints.dat',3);
DX = XMAX / NPIX;

[xgrid,ygrid] = meshgrid(0:DX:XMAX,0:DY:YMAX);


u = reshape(Data(:,3),[NPJ, NPI]);      v = reshape(Data(:,4),[NPJ, NPI]);
p = reshape(Data(:,5),[NPJ, NPI]);      T = reshape(Data(:,6),[NPJ, NPI]);
rho = reshape(Data(:,7),[NPJ, NPI]);    mu = reshape(Data(:,8),[NPJ, NPI]);

eps = reshape(Data(:,11),[NPJ, NPI]);


% plot T
figure(3)
surf(X, Y, eps)
hold on

colorbar
view(0,90)


%% Load data from constraints file



%% Data analysis


figure(2)
quiver(X,Y,u,v)

Q = AddedHeat([0,0.15,0.3,0.45]);


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
    %shading interp
    hold on
    set(gca,'ZScale','log')
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
        fprintf('Q_added at [x=%4.2f] = %f [W] Tavg = %6.2f [K]\n',X_pos(i),Q(i),mean(T(:,I)))
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
