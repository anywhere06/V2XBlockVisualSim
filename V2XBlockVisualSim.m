clear;
clc;
global R X Y;
global cnodes;
lambda_c = 0.000005;           % Central node density per m^2

Y = 1e4;                       % Length (vertical) of R^2 (km)
X = 1e4;                       % Width (horizontal) of R^2 (km)

show = 0;
Nnode = 50;                     % Number of member nodes Radius
NnodePoisson = 10;
Pchildren = .3;

R = 1e2;                        % Radius of each central node (m)

conn = 10;              % # of connections for each central node
FIG_ON = 1;                     %Determines if initial figure is displayed
timeTick = 1;                   %Time taken per frame
frameRate = 300;                 %Frame rate in ms

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
failedmember = 0;
members = 0;
%% cnode distribution
% central nodes
cnodeX = X*rand(1,Ncnode);          %Uniform Distribution 0-X
cnodeY = Y*rand(1,Ncnode);          %Uniform Distribution 0-Y
cnodes = complex(cnodeX, cnodeY);   %Sets array of coordinates for center of central nodes complex: Real = X,Imag = Y
cChildren = find(Pchildren>rand(1,length(cnodes)));
for i = 1:length(cChildren)
    mem = poissrnd(NnodePoisson);
    members = members + mem;
    for q = 1:mem
        angle = 360*rand;
        radius = Nnode*rand;
        mnodeX(i,q) = radius*cosd(angle)+cnodeX(cChildren(i));
        mnodeY(i,q) = radius*sind(angle)+cnodeY(cChildren(i));
        mnodeC(i,q) = 'b';
    end
end
%%Frame calculations
flags = zeros(1,length(cnodes));
start = ceil(length(cnodes)*rand);
flags(start) = 1;
for i = 1:length(flags)
    scat(i,:) = [0 0 1];
end
scat(start,:) = [0 1 0];
Frames{1,1} = scat;
Frames{2,1} = {0};
Frames{3,1} = {mnodeX,mnodeY,mnodeC};
w = 1;
t = 1;
go = 1;
%%Connections
for i = 1:length(cnodes)
    dist = zeros(1,length(cnodes));
    for ii = 1:length(cnodes)
        if i == ii
            continue;
        end
        dist(ii) = distance(cnodeX(i),cnodeY(i),cnodeX(ii),cnodeY(ii));
    end
    distsort = sort(dist);
    for l = 2:(conn+1)
        choose(l-1) = find(dist==distsort(l));
    end
    connections{1,i} = choose;
end
%%
tic
k = 1;
Done = 0;
failcount = 0;
while(go)
    if ~length(find(flags==0))
        go = 0;
    end
    if (toc>=maxt)
        fprintf('runtime error\n');
        break;
    end
    t = t+1;
    Trans = find(flags==1);
    r = 1;
    linegraph = cell(1,1);
    for i = 1:length(Trans)
        if find(Done==Trans(i))
            continue;
        end
        for ii = 1:length(connections{1,Trans(i)})
            val = connections{1,Trans(i)}(ii);
            if (~isempty(find(Done==val))||(flags(val)))
                continue;
            end
            Done(k) = Trans(i);
            k = k+1;
            fail = rand<p;
            flags(val) = fail+1;
            if find(cChildren==val)
                nodechosen = find(cChildren==val);
                if fail
                    mnodeC(nodechosen,1:length(find(mnodeX(nodechosen,:)))) = 'r';
                    failedmember = failedmember+ length(find(mnodeX(nodechosen,:)));
                else
                    for q = 1:length(find(mnodeX(nodechosen,:)))
                        newfail = rand<p;
                        if newfail
                            mnodeC(nodechosen,q) = 'r';
                            failedmember = failedmember+1;
                        else
                            mnodeC(nodechosen,q) = 'g';
                        end
                    end
                end
            end
        if fail
            linecolor = 'red';
            failcount = failcount +1;
        else
            linecolor = 'gre';
        end
        linegraph{r,1} = val;
        linegraph{r,2} = Trans(i);
        linegraph{r,3} = linecolor;
        r = r+1;
        end
    end
    
    for e = 1:length(flags)
        switch(flags(e))
            case 0,
                scat(e,:)= [0 0 1]; %blue
            case 1,
                scat(e,:)= [0 1 0]; %green
            case 2,
                scat(e,:)= [1 0 0]; %red
        end
    end
    %     plotx(w,1) = real(cnodes(pickt));
    %     plotx(w,2) = real(cnodes(pickr));
    %     ploty(w,1) = imag(cnodes(pickt));
    %     ploty(w,2) = imag(cnodes(pickr));
    %     plotc(w,:) = linecolor;
    w = w+1;
    Frames{1,t} = scat;
    Frames{2,t} = linegraph;
    Frames{3,t} = {mnodeX,mnodeY,mnodeC};
end
%%
if FIG_ON
    figure;
    hold on;
    plotSetup(VarTitle);
    plotMap(Frames{1,1},Frames{2,1});
    Mnodes = scatter(mnodeX(find(mnodeX)),mnodeY(find(mnodeY)),10,'b','Filled');
end
flag = 1;
flagskip = 1;
while flag
    fprintf('Enter a 0 to terminate\n');
    flag = input('Enter a 1 to show the animation\nEnter in a 2 to zoom into specific areas\n');
    if flag == 1
        clf;
        plotSetup(VarTitle);
        %        fig = line(Frames{2,t}{1}',Frames{2,t}{2}','Color','non');
        nodes = scatter(real(cnodes)',imag(cnodes)',50,Frames{1,1},'Filled');
        Mnodes = scatter(mnodeX(find(mnodeX))',mnodeY(find(mnodeY))',10,'b','Filled');
%         nodes = scatter(Frames{3,:}(1:)
        for i = 2:t
            %             for e = 1:size(Frames{2,i}{3},1)
            %                 fig(e).Color = char(Frames{2,i}{3}(e,:));
            %             end
            nodes.CData = Frames{1,i};
            Data = Frames{3,i}{1,3}(find(Frames{3,i}{1,3}));
            for m = 1:length(Data)
                if Data(m) == 'r'
                    DataVector(m,:) = [1 0 0];
                elseif Data(m) == 'g'
                    DataVector(m,:) = [0 1 0];
                elseif Data(m) == 'b'
                    DataVector(m,:) = [0 0 1];
                end
            end
            Mnodes.CData = DataVector;
            if i == t
            else
                lines = Frames{2,i};
                for e = 1:size(lines,1)
                    plot([cnodeX(lines{e,1}),cnodeX(lines{e,2})],[cnodeY(lines{e,1}),cnodeY(lines{e,2})],'Color',lines{e,3})
                end
            end
            pause(frameRate/1000);
        end
        
        % Statistics
        fprintf('Statistics:\n%d Central nodes: %d fail %d succeed %f % success rate\n',length(cnodes),failcount,length(cnodes)-failcount),100-100*(failcount/length(cnodes));
        fprintf('\n%d Member nodes: %d fail %d succeed %f %\n',length(find(mnodeX(:))),failedmember,length(find(mnodeX(:)))-failedmember),100-100*failedmember/length(find(mnodeX(:)));
        fprintf('\n%d time slots taken for full dissemination\n',t);
        
        % break
        break;
        
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
        %         fig = line(Frames{2,t}{1}',Frames{2,t}{2}','Color','non');
        nodes = scatter(real(cnodes)',imag(cnodes)',50,Frames{1,1},'Filled');
        %         for i = 2:t
        %             for e = 1:size(Frames{2,i}{3},1)
        %                 fig(e).Color = char(Frames{2,i}{3}(e,:));
        %             end
        %             nodes.CData = Frames{1,i};
        %             pause(frameRate/1000);
        %         end
    else
        break;
    end
end

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