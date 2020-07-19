clear;
clc;
global R X Y;
global cnodes;
lambda_c = 0.0000001;           % Central node density per m^2

Y = 1e4;                       % Length (vertical) of R^2 (km)
X = 1e4;                       % Width (horizontal) of R^2 (km)

show = 0;
Nnode = 10;                     % Number of member nodes per central node

R = 1e2;                        % Radius of each central node (m)

FIG_ON = 1;                     %Determines if initial figure is displayed
timeTick = 1;                   %Time taken per frame
frameRate = 30;                 %Frame rate in ms

Vvel = 100;                     %Vehicle velocity (m/s)
Vrad = 2*R;                       %Vehicle Effective Radius (m)
global VehW VehH
VehW = 25;
VehH = 50;

p = .1;                         %Node failure rate
maxt = 20;                      %Max runtime(seconds)

%% Initializations and Calculations
% Initial Calculations
Ncnode = ceil(lambda_c*X*Y);  %Calculates the number of central nodes needed
VarTitle = ['\lambda_{c} = ' num2str(lambda_c) ' / km^{2},   N_{nodes} = ' num2str(Ncnode) '/   R = ' num2str(R)];
% ^ Title for initial variables
connections = cell(2,Ncnode);
% Time calculations

%% cnode distribution
% central nodes
cnodeX = X*rand(1,Ncnode);          %Uniform Distribution 0-X
cnodeY = Y*rand(1,Ncnode);          %Uniform Distribution 0-Y
cnodes = complex(cnodeX, cnodeY);   %Sets array of coordinates for center of central nodes complex: Real = X,Imag = Y
i = 1;
%%Frame calculations
flags = zeros(1,length(cnodes));
start = ceil(length(cnodes)*rand);
flags(start) = 1;
for i = 1:length(flags)
    scat(i,:) = [0 0 1];
end
scat(start,:) = [0 1 0];
Frames{1,1} = scat;
plotx = zeros(1,2);
ploty = zeros(1,2);
plotc = zeros(1,3);
plotc = ['n' 'o' 'n'];
Frames{2,1} = {plotx,ploty,plotc};
w = 1;
t = 1;
go = 1;
tic
while(go)
    if ~length(find(flags==0))
        go = 0;
    end
    if (toc>=maxt)
        fprintf('runtime error\n');
        break;
    end
    Trans = find(flags==1);
    pickt = Trans(ceil(rand*length(Trans)));
    pickr = ceil(rand*length(cnodes));
    if flags(pickr)>=1
        continue;
    end
    t = t+1;
    fail = rand<p;
    flags(pickr) = fail+1;
    if fail
        linecolor = 'red';
    else
        linecolor = 'gre';
    end
    for i = 1:length(flags)
        switch(flags(i))
            case 0,
                scat(i,:)= [0 0 1]; %blue
            case 1,
                scat(i,:)= [0 1 0]; %green
            case 2,
                scat(i,:)= [1 0 0]; %red
        end
    end
    plotx(w,1) = real(cnodes(pickt));
    plotx(w,2) = real(cnodes(pickr));
    ploty(w,1) = imag(cnodes(pickt));
    ploty(w,2) = imag(cnodes(pickr));
    plotc(w,:) = linecolor;
    w = w+1;
    Frames{1,t} = scat;
    Frames{2,t} = {plotx,ploty,plotc};
end
%%
if FIG_ON
    figure;
    hold on;
    plotSetup(VarTitle);
    plotMap(Frames{1,1},Frames{2,1});
end
flag = 1;
flagskip = 1;
while flag
    flag = input('Enter a 1 to show the animation\nEnter in a 2 to zoom into specific areas\n');
    if flag == 1
        clf;
        plotSetup(VarTitle);
        fig = line(Frames{2,t}{1}',Frames{2,t}{2}','Color','non');
        nodes = scatter(real(cnodes)',imag(cnodes)',50,Frames{1,1},'Filled');
        for i = 2:t
            for e = 1:size(Frames{2,i}{3},1)
                fig(e).Color = char(Frames{2,i}{3}(e,:));
            end
            nodes.CData = Frames{1,i};
            pause(frameRate/1000);
        end
    elseif flag == 2
        [x] = input('Enter in window range:\nx values: [min,max] = ');
        if x(1) >= x(2) || length(x)~=2
            fprintf('invalid input');
            continue;
        end
        [y] = input('Enter in window range:\ny values: [min,max] = ');
        if y(1) >= y(2) || length(y)~=2
            fprintf('invalid input');
            continue;
        end
        clf;
        plotSetup(VarTitle);
        xlim(x);
        ylim(y);
        fig = line(Frames{2,t}{1}',Frames{2,t}{2}','Color','non');
        nodes = scatter(real(cnodes)',imag(cnodes)',50,Frames{1,1},'Filled');
        for i = 2:t
            for e = 1:size(Frames{2,i}{3},1)
                fig(e).Color = char(Frames{2,i}{3}(e,:));
            end
            nodes.CData = Frames{1,i};
            pause(frameRate/1000);
        end
    else
        break;
    end
end
% while flag
%     flag = input('Enter a 1 to show the animation\n');
%     if flag == 1
%         for i = 1:t
%             clf;
%             plotSetup(VarTitle);
%             plotMap(Frames{2,i},Frames{3,i});
%             pause(frameRate/1000);
%         end
%     else
%         break;
%     end
% end

%%
function plotSetup(plotTitle)
global X R Y;  %retrieves global variables
hold on
grid on
box on
xlim([-R X+R])
ylim([-R Y+R])
xlabel('x')
ylabel('y')
title(plotTitle)
end

% function plotMap(connections,scat)
% global cnodes;
% ax = gca;
% for i = 1:length(cnodes)
%     plotx = zeros(length(connections{1,i}),2);
%     ploty = zeros(length(connections{1,i}),2);
%     plotx(:,1) = real(cnodes(i));
%     plotx(:,2) = real(cnodes(connections{1,i}));
%     ploty(:,1) = imag(cnodes(i));
%     ploty(:,2) = imag(cnodes(connections{1,i}));
%     ax.ColorOrder = connections{2,i};
%     line(plotx',ploty');
% end
% scatter(real(cnodes)',imag(cnodes)',50,scat,'Filled');
% end
function plotMap(scat,lineval)
global cnodes col;
scatter(real(cnodes)',imag(cnodes)',50,scat,'Filled');
end

function out = distance(x1,y1,x2,y2)
out = sqrt((x1-x2)^2+(y1-y2)^2);
end