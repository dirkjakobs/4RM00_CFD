clear all; clc; close all;

%% Import data
Data = load('output.dat');

X = Data(:,1);  Y = Data(:,2);      % load x and y data
NPJ = sum(X == X(1));               % Get original matrix size from x and y
NPI = sum(Y == Y(1));
X = reshape(X,[NPJ, NPI]);
Y = reshape(Y,[NPJ, NPI]);
u = reshape(Data(:,3),[NPJ, NPI]);

%plot with imagplot command: imagplot( X, Y, data, 'title');
FIG = imagplot(X,Y,u,'Velocity in u-direction [m/s]');

Q = AddedHeat([0,0.15,0.3,0.45]);

%% Functions

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
    TPLOT = imagplot(X,Y,T,'Temperature [K]');
    
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
        figure(TPLOT)
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

function FIG = imagplot(X,Y,data,titel)
    XMIN = min(X(:));
    XMAX = max(X(:));
    YMIN = min(Y(:));
    YMAX = max(Y(:));    
    FIG = figure();
    surf(X, Y, data)
    hold on

    shading interp
    %colormap(jet(256))
    colorbar
    view(0,90)
    axis equal
    grid off
    xlim([-0.05*XMAX 1.05*XMAX])
    ylim([-YMAX 2*YMAX])
    set(gca,'YTick',[0 : 0.015 : 0.045]);
    set(gca,'XTick',[0 : 0.15 : 0.45]);
    title(titel)

end
