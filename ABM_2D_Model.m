function [TC_ColumnAverages, EC_ColumnAverages, TC_ColumnAverage, EC_ColumnAverage, NumberSelfLoops, AverageSelfLoops, NumberTipTipAnastomoses, NumberTipSproutAnastomoses, NumberofTipCellsTime, NumberofBranchEvents, SproutLengths, AverageSproutLength, StDevSproutLength, AverageBranchLength, StDevBranchLength, NumberBacktrackingLoops, BranchLengths, SelfLoopsNoBacktracking, NumberofLargeSelfLoops, Networks, PerfusedNetworks, Perfused_ColumnAverages, Perfused_ColumnAverage, Average2DTipCellNetwork, Average2DStalkCellNetwork, SurvivingSprouts_LateralMovement, SurvivingSprouts_LengthofSprout, SurvivingSprouts_RatioofLateralMovement, SproutsAliveYesNo, ReachedTumorYesNo, TC_Solution, EC_Solution, SummaryStatistics, Parameters, Time] = ABM_2D_Model
% Computes multiple realizations of Pillay_2017_CA.m (See other function in
% folder)
%--------------------------------------------------------------------------
%% Name Solution Files
% filename = '12may2020_PillayCAModel_Fixed_P_p0_P_m1_k100_NonlinearXY_TAFField_NoBranchingOrAnastomosis';
filename = 'ABM_Data';
%% Parameter values, Pre-allocation of data: 
%------------ Parameters-------------------------------------
Parameters.N = 200 + 1; % Number of Lattice points
Parameters.h = 1/(Parameters.N-1); % Spatial step size of lattice
Parameters.Q = 2*Parameters.N+50; % Number of cells to track for each branch/sprout
Parameters.P_m = 1; % Probability of Movement
Parameters.P_p = 0; % Probability of branching (note: proportional to TAF concentration also)
Parameters.dt = 1/160; % Time Step
% Parameters.dt = 1/320;
% Parameters.dt = 0.01;
Parameters.k = 100; % ABM parameter related to sensitivity of tip cells to TAF
Parameters.Tf = 2; % Final Time
% Parameters.Tf = 320*Parameters.dt;

Parameters.a_n = false; % Does Tip-to-Tip Anastomosis Occur? (Binary)
Parameters.a_e = false; % Does Tip-to-Sprout Anastomosis Occur? (Binary)
Parameters.EC_Branching_YesNo = false; % Will you add on a stalk cell after branching or not? (Binary)

Parameters.M = 1000; % Number of Realizations

% Add on number of realizations to filename:
filename = [filename,'_',num2str(Parameters.M),'Realizations'];
%------------ TAF Concentration-------------------------------------
% c(x,y) = x
c = repmat(linspace(0, 1, Parameters.N), Parameters.N, 1); % TAF Concentration Matrix

% c(x,y) = exp((-(x-1).^2)-(y-1).^2/2.5)
% c = exp((-(linspace(0,1,Parameters.N)-1).^2)./1+(-(linspace(0,1,Parameters.N)'-1).^2.)./2.5);

% c(x,y) = xy
% c = linspace(0,1,201).*linspace(0,1,201)';

% c(x,y) = 0.5*(x + y)
% c = 0.5*(linspace(0,1,Parameters.N)' + linspace(0,1,Parameters.N));

% c(x,y) = 1/2*y + 1/3*x
% c = 0.5*(linspace(0,1,201)') + 1/3*(linspace(0,1,201));

% c(x,y) = 1/2*x^2
% c = repmat(0.5*(linspace(0, 1, Parameters.N)).^2, Parameters.N, 1); % TAF Concentration Matrix

% c(x,y) = 1-(x-0.5)^2-(y-0.5)^2
% c = 1 - (linspace(0,1,201)-0.5).^2 - (linspace(0,1,201)'-0.5).^2;

% Graph c(x,y) to make sure that you have the correct TAF concentration
% profile:
surf(linspace(0, 1, Parameters.N), linspace(0, 1, Parameters.N), c);
xlabel('x')

%------------ Pre-allocate solutions -------------------------------------
Time = (0:Parameters.dt:Parameters.Tf)';
Lt = length(Time);

TC_ColumnAverages = zeros(Parameters.N, Lt, Parameters.M);
EC_ColumnAverages = zeros(Parameters.N, Lt, Parameters.M);
Perfused_ColumnAverages = zeros(Parameters.N, Lt, Parameters.M);
NumberSelfLoops = zeros(Lt,Parameters.M, 'uint16');
NumberTipTipAnastomoses = NumberSelfLoops;
NumberTipSproutAnastomoses = NumberSelfLoops;
NumberofTipCellsTime = zeros(Lt, Parameters.M, 'uint16');
NumberofBranchEvents = NumberSelfLoops;
NumberBacktrackingLoops = NumberSelfLoops;

SproutLengths = zeros(Parameters.Q, Lt, Parameters.M, 'uint16');
BranchLengths = SproutLengths;
NumberofLargeSelfLoops = zeros(Parameters.M,1, 'uint16');

SurvivingSprouts_LateralMovement = [];
SurvivingSprouts_LengthofSprout = [];
SurvivingSprouts_RatioofLateralMovement = [];

SproutsAliveYesNo = zeros(Parameters.M,1,'logical');
ReachedTumorYesNo = zeros(Parameters.M,1,'logical');
Average2DTipCellNetwork = zeros(Parameters.N,Parameters.N,Lt);
Average2DStalkCellNetwork = zeros(Parameters.N,Parameters.N,Lt);

Networks = cell(Parameters.M);
PerfusedNetworks = Networks;
TC_Solution = Networks;
EC_Solution = Networks;

Num_realizations = Parameters.M;
%--------------------Solve the ABM for M realizations----------------------
parfor R = 1:Parameters.M
    fprintf(['Solving for Realization ',num2str(R),'\n'])
    [~, ~, TC_ColumnAverages(:,:,R), EC_ColumnAverages(:,:,R), Example_Networks, NumberSelfLoops(:,R), NumberTipTipAnastomoses(:,R), NumberTipSproutAnastomoses(:,R), NumberofTipCellsTime(:,R), NumberofBranchEvents(:,R), SproutLengths(:,:,R), NumberBacktrackingLoops(:,R), BranchLengths(:,:,R), NumberofLargeSelfLoops(R), Example_PerfusedNetworks, Perfused_ColumnAverages(:,:,R), ~, Example_TC_Solution, Example_EC_Solution, ~, ~, ~, SproutsAliveYesNo(R), ReachedTumorYesNo(R)] = Pillay_2017_CA(c,Parameters);
    Average2DTipCellNetwork = Average2DTipCellNetwork + double(Example_TC_Solution);
    Average2DStalkCellNetwork = Average2DStalkCellNetwork + double(Example_EC_Solution);
    if R == Num_realizations
        Networks{R} = Example_Networks;
        PerfusedNetworks{R} = Example_PerfusedNetworks;
        TC_Solution{R} = Example_TC_Solution;
        EC_Solution{R} = Example_EC_Solution;
    else
        Networks{R} = [];
        PerfusedNetworks{R} = [];
        TC_Solution{R} = [];
        EC_Solution{R} = [];
    end % if
end % for R

Networks = Networks{end};
PerfusedNetworks = PerfusedNetworks{end};
TC_Solution = TC_Solution{end};
EC_Solution = EC_Solution{end};

%% Take averages
% Gives column x Time matrices of size (N x Lt)
Average2DTipCellNetwork = 1/Parameters.M*Average2DTipCellNetwork;
Average2DStalkCellNetwork = 1/Parameters.M*Average2DStalkCellNetwork;
TC_ColumnAverage = mean(TC_ColumnAverages,3);
EC_ColumnAverage = mean(EC_ColumnAverages,3);
Perfused_ColumnAverage = mean(Perfused_ColumnAverages,3);
% These give Lt x 1 matrices, but the final entry is empty (because you're
% counting the number of stuff that happens WITHIN a time step, not AFTER
% it
TotalLoopsFormedNoSelfLoops = double(sum(NumberTipTipAnastomoses + NumberTipSproutAnastomoses, 1));
AverageSelfLoops = mean(double(NumberSelfLoops),2);

AverageSproutLength = mean(double(SproutLengths), 3);
StDevSproutLength = std(double(SproutLengths),0,3);
AverageBranchLength = mean(double(BranchLengths), 3);
StDevBranchLength = std(double(BranchLengths),0,3);

%% Plotting
figure;
imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', (Networks(:,:,end)>0) + 2*(PerfusedNetworks(:,:,end)>0) + logical(TC_Solution(:,:,end)));
axis tight
colormap([1 1 1; 218/255, 112/255, 214/255; 0 0 1; 1 0 0]);
colorbar;
title(['CA Model Network Created (Tip and Stalk Cells), Single Realization, t = ',num2str(Parameters.Tf),' s'])
xlabel('x');
ylabel('y');
figurepalette('show');
saveas(gcf, [filename,'_1_ExampleNetworkwithOverlay@Tf'], 'fig');
saveas(gcf, [filename,'_1_ExampleNetworkwithOverlay@Tf'], 'png');

figure;
imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', Networks(:,:,end));
axis tight
colormap(jet);
colorbar;
title(['CA Model Network Created (Tip and Stalk Cells), Single Realization, t = ',num2str(Parameters.Tf),' s'])
xlabel('x');
ylabel('y');
figurepalette('show')
saveas(gcf, [filename,'_2_ExampleNetwork@Tf'], 'fig');
saveas(gcf, [filename,'_2_ExampleNetwork@Tf'], 'png');

figure;
imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', PerfusedNetworks(:,:,end));
map = [1 1 1; 1 0 0];
axis tight
colormap(map);
% colorbar;
title(['CA Model Perfused Network Created, Single Realization, t = ',num2str(Parameters.Tf),' s'])
xlabel('x');
ylabel('y');
figurepalette('show')
saveas(gcf, [filename,'_3_ExamplePerfusedNetwork@Tf'], 'fig');
saveas(gcf, [filename,'_3_ExamplePerfusedNetwork@Tf'], 'png');

figure;
imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', Average2DTipCellNetwork(:,:,end));
colormap(jet);
axis tight
colorbar;
title(['Average 2D Tip Cell Network Created by CA Model over ',num2str(Parameters.M),' realizations, t = ',num2str(Parameters.Tf),' s']);
xlabel('x');
ylabel('y');
figurepalette('show')
saveas(gcf, [filename,'_4_Average2DTCSolution@Tf'], 'fig');
saveas(gcf, [filename,'_4_Average2DTCSolution@Tf'], 'png');

figure;
imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', Average2DStalkCellNetwork(:,:,end));
colormap(jet);
axis tight
colorbar;
title(['Average 2D Stalk Cell Network Created by CA Model over ',num2str(Parameters.M),' realizations, t = ',num2str(Parameters.Tf),' s']);
xlabel('x');
ylabel('y');
figurepalette('show')
saveas(gcf, [filename,'_5_Average2DECSolution@Tf'], 'fig');
saveas(gcf, [filename,'_5_Average2DECSolution@Tf'], 'png');

figure;
plot(0:Parameters.h:Parameters.h*(Parameters.N-1), TC_ColumnAverage(:,32:32:end));
title(['CA Model Tip Cell Column Average over ',num2str(Parameters.M),' Realizations in increments of 0.5 seconds']);
xlabel('x');
ylabel('N(x,t)');
figurepalette('show')
saveas(gcf, [filename,'_6_TC_ColumnAverage'], 'fig');
saveas(gcf, [filename,'_6_TC_ColumnAverage'], 'png');

figure;
plot(0:Parameters.h:Parameters.h*(Parameters.N-1), EC_ColumnAverage(:,32:32:end));
title(['CA Model Stalk Cell Column Average over ',num2str(Parameters.M),' Realizations in increments of 0.2 seconds']);
xlabel('x');
ylabel('E(x,t)');
figurepalette('show')
saveas(gcf, [filename,'_7_EC_ColumnAverage'], 'fig');
saveas(gcf, [filename,'_7_EC_ColumnAverage'], 'png');

figure;
plot(0:Parameters.h:Parameters.h*(Parameters.N-1), Perfused_ColumnAverage(:, 32:32:end));
title(['CA Model Perfused Cell Column Average over ',num2str(Parameters.M),' Realizations in increments of 0.2 seconds']);
xlabel('x');
ylabel('Average Number of Cells');
figurepalette('show');
saveas(gcf, [filename,'_8_PerfusedCells_ColumnAverage'], 'fig');
saveas(gcf, [filename,'_8_PerfusedCells_ColumnAverage'], 'png');

figure;
plot(Time, mean(NumberofTipCellsTime,2))
title(['Number of Tip Cells per Time, Average over ',num2str(Parameters.M),' Realizations'])
xlabel('t')
ylabel('N_{TC}')
figurepalette('show')
saveas(gcf, [filename,'_9_NumberTipCellsVsTime'], 'fig');
saveas(gcf, [filename,'_9_NumberTipCellsVsTime'], 'png');

figure;
plot(Time, mean(NumberSelfLoops,2))
title(['Number of Self Loops that would have formed per Time, Average over ',num2str(Parameters.M),' Realizations'])
xlabel('t')
ylabel('N_{self loops}')
figurepalette('show')
saveas(gcf, [filename,'_10_NumberSelfLoopsVsTime'], 'fig');
saveas(gcf, [filename,'_10_NumberSelfLoopsVsTime'], 'png');

figure;
bar(1:size(AverageSproutLength, 1), AverageSproutLength(:,end))
title(['Average Sprout Length (original + branches) over ',num2str(Parameters.M),' Realizations @ t = ',num2str(Parameters.Tf),'s']);
xlabel('Number of Cells in Sprout');
ylabel('Average Number of Sprouts')
figurepalette('show');
saveas(gcf, [filename,'_11_SproutLengthHistogram'], 'fig');
saveas(gcf, [filename,'_11_SproutLengthHistogram'], 'png');

figure;
bar(1:size(AverageBranchLength,1), AverageBranchLength(:,end))
title(['Average Length of Branches (not original sprout) over ',num2str(Parameters.M),' realizations @ t = ',num2str(Parameters.Tf),'s']);
xlabel('Number of Cells in Branch');
ylabel('Average Number of Branches');
figurepalette('show');
saveas(gcf, [filename,'_12_BranchLengthHistogram'], 'fig');
saveas(gcf, [filename,'_12_BranchLengthHistogram'], 'png');

figure;
histogram(sum(NumberofBranchEvents,1), 'BinMethod', 'integers');
title(['Histogram of Branching Events from ',num2str(Parameters.M),' Realizations'])
xlabel('Number of Branching Events')
ylabel('N_{realizations}')
figurepalette('show')
saveas(gcf, [filename,'_13_BranchEventHistogram'], 'fig');
saveas(gcf, [filename,'_13_BranchEventHistogram'], 'png');

%% Statistics 
fprintf(['The average total number of loops formed over ',num2str(Parameters.Tf),' seconds and ',num2str(Parameters.M),' realizations (no self-loops over 10 cells) is \n']);
fprintf([num2str(mean(TotalLoopsFormedNoSelfLoops)),'\n']);

fprintf('The standard deviation from this mean value is \n');
fprintf([num2str(std(TotalLoopsFormedNoSelfLoops)),'\n']);

fprintf(['The average number of self-loops that would have formed over ',num2str(Parameters.M),' realizations and ',num2str(Parameters.Tf),' seconds is \n']);
fprintf([num2str(mean(sum(NumberSelfLoops, 1))),'\n']);

fprintf('The standard deviation from this mean value is \n');
fprintf([num2str(std(sum(NumberSelfLoops, 1))),'\n']);

fprintf(['The average proportion of loops formed from tip-tip anastomosis events to total loops over ',num2str(Parameters.M),' realizations and ',num2str(Parameters.Tf),' seconds is \n']);
fprintf([num2str(mean(sum(NumberTipTipAnastomoses, 1)./TotalLoopsFormedNoSelfLoops)),'\n']);

fprintf('The standard deviation from this mean value is \n');
fprintf([num2str(std(sum(NumberTipTipAnastomoses, 1)./TotalLoopsFormedNoSelfLoops)),'\n']);

fprintf(['The average proportion of loops formed from tip-sprout anastomosis events to total loops over ',num2str(Parameters.M),' realizations and ',num2str(Parameters.Tf),' seconds is \n']);
fprintf([num2str(mean(sum(NumberTipSproutAnastomoses, 1)./TotalLoopsFormedNoSelfLoops)),'\n']);

fprintf('The standard deviation from this mean value is \n');
fprintf([num2str(std(sum(NumberTipSproutAnastomoses, 1)./TotalLoopsFormedNoSelfLoops)),'\n']);

fprintf(['The average proportion of self-loops from backtracking to the number of self loops over ',num2str(Parameters.M),' realizations and ',num2str(Parameters.Tf),' seconds is \n']);
fprintf([num2str(mean(sum(NumberBacktrackingLoops, 1)./sum(NumberSelfLoops, 1))),'\n']);

fprintf('The standard deviation from this mean value is \n');
fprintf([num2str(std(sum(NumberBacktrackingLoops, 1)./sum(NumberSelfLoops, 1))),'\n']);

SelfLoopsNoBacktracking = NumberSelfLoops-NumberBacktrackingLoops;

fprintf(['The average total number of branches that formed over ',num2str(Parameters.Tf),' seconds and ',num2str(Parameters.M),' realizations is \n']);
fprintf([num2str(mean(sum(NumberofBranchEvents,1))),'\n']);

fprintf('The standard deviation from this mean value is \n');
fprintf([num2str(std(sum(NumberofBranchEvents,1))),'\n']);

fprintf(['The average number of self-loops that would have formed within ',num2str(Parameters.Tf),' seconds and ',num2str(Parameters.M),' realizations that are 11 cells long or more is \n']);
fprintf([num2str(mean(NumberofLargeSelfLoops)),'\n'])

fprintf('The standard deviation from this mean value is \n');
fprintf([num2str(std(double(NumberofLargeSelfLoops))),'\n'])

fprintf('The probability of having at least one sprout survive is \n');
fprintf([num2str(mean(SproutsAliveYesNo)),'\n'])

fprintf('The standard deviation from this mean value is \n');
fprintf([num2str(std(SproutsAliveYesNo)),'\n'])

fprintf(['The probability of reaching the tumor within ',num2str(Parameters.Tf),' s is \n']);
fprintf([num2str(mean(ReachedTumorYesNo)),'\n'])

fprintf('The standard deviation from this mean value is \n');
fprintf([num2str(std(ReachedTumorYesNo)),'\n'])

fprintf(['The average number of tip cells that have survived after ',num2str(Parameters.Tf),' s is \n']);
fprintf([num2str(mean(NumberofTipCellsTime(end,:))),'\n'])

fprintf('The standard deviation from this mean value is \n');
fprintf([num2str(std(double(NumberofTipCellsTime(end,:)))),'\n']);

% Also Output these statistics as a table for easier data collection:
stats_descriptions = {'Number of Total Loops'; ...
    'Proportion of Tip-Tip Anastomosis Events'; ...
    'Proportion of Tip-Sprout Anastomosis Events';...
    'Number of Branch Events';...
    'Number of Tip Cells at Final Time';...
    'Proportion of Realizations with Surviving Sprouts at Final Time';...
    'Proportion of Realizations that Reached the Tumor';...
    'Number of Successful Self Loops > 10 cells';...
    'Number of Failed Self Loops <= 10 cells (repeats possible)';...
    'Proportion of Failed Self Loops Resulting from Backtracking'
    };
Mean = [mean(TotalLoopsFormedNoSelfLoops); ...
    mean(sum(NumberTipTipAnastomoses, 1)./TotalLoopsFormedNoSelfLoops);...
    mean(sum(NumberTipSproutAnastomoses, 1)./TotalLoopsFormedNoSelfLoops);...
    mean(sum(NumberofBranchEvents,1));...
    mean(NumberofTipCellsTime(end,:));...
    mean(SproutsAliveYesNo);...
    mean(ReachedTumorYesNo);...
    mean(NumberofLargeSelfLoops);...
    mean(sum(NumberSelfLoops, 1));...
    mean(sum(NumberBacktrackingLoops, 1)./sum(NumberSelfLoops, 1))
    ];
Standard_Deviation = [std(TotalLoopsFormedNoSelfLoops);...
    std(sum(NumberTipTipAnastomoses, 1)./TotalLoopsFormedNoSelfLoops);...
    std(sum(NumberTipSproutAnastomoses, 1)./TotalLoopsFormedNoSelfLoops);...
    std(sum(NumberofBranchEvents,1));...
    std(double(NumberofTipCellsTime(end,:)));...
    std(SproutsAliveYesNo);...
    std(ReachedTumorYesNo);...
    std(double(NumberofLargeSelfLoops));...
    std(sum(NumberSelfLoops, 1));...
    std(sum(NumberBacktrackingLoops, 1)./sum(NumberSelfLoops, 1))
    ];
SummaryStatistics = table(Mean, ...
    Standard_Deviation, 'RowNames', stats_descriptions);
%% Write CSV Files for use in 2D/1D PDE Models
i = find(Time==0.2);
X = reshape(repmat(linspace(0,1,Parameters.N),Parameters.N,1),Parameters.N^2,1);
Y = repmat(linspace(0,1,Parameters.N)', Parameters.N, 1);
% Save 2D P--ABM results at t = 0.2:
csvwrite([filename,'_Average2DStalkCellNetwork_t2.csv'], [X, Y, reshape(Average2DStalkCellNetwork(:,:,i), Parameters.N^2, 1)]);
csvwrite([filename,'_Average2DTipCellNetwork_t2.csv'], [X, Y, reshape(Average2DTipCellNetwork(:,:,i), Parameters.N^2, 1)]);
for t = 4:2:20
i = 1 + 16*t;
csvwrite([filename,'_Average2DStalkCellNetwork_t',num2str(t),'.csv'], [X, Y, reshape(Average2DStalkCellNetwork(:,:,i), Parameters.N^2, 1)]);
csvwrite([filename,'_Average2DTipCellNetwork_t',num2str(t),'.csv'], [X, Y, reshape(Average2DTipCellNetwork(:,:,i), Parameters.N^2, 1)]);
end

% Save column averaged 2D P--ABM results at t = 0.2, 0.4, ..., 2:
csvwrite([filename, '_1DTipCellNetworkIC_t2.csv'], [linspace(0,1,Parameters.N)', TC_ColumnAverage(:,i)]);
csvwrite([filename, '_1DStalkCellNetworkIC_t2.csv'], [linspace(0,1,Parameters.N)', EC_ColumnAverage(:,i)]);
for t = 4:2:20
i = 1 + 16*t;
csvwrite([filename, '_Average1DTipCellNetwork_t',num2str(t),'.csv'], [linspace(0,1,Parameters.N)', TC_ColumnAverage(:,i)]);
csvwrite([filename, '_Average1DStalkCellNetwork_t',num2str(t),'.csv'], [linspace(0,1,Parameters.N)', EC_ColumnAverage(:,i)]);
end
%% Save Workspace for Future Reference
save([filename,'.mat']);
%% Write Video Files using the Data You've Collected
PillayCA_WriteVideos(filename, Networks, PerfusedNetworks, BranchLengths, SproutLengths, Time, Average2DStalkCellNetwork, Average2DTipCellNetwork, TC_Solution, Parameters);
end % function Pillay_2017_CA_MultipleRealizations.m
%--------------------------------------------------------------------------
%% Subfunctions
function [Time, Sprouts, TC_ColumnAverage, EC_ColumnAverage, Network, NumberofSelfLoops, NumberTipTipAnastomoses, NumberTipSproutAnastomoses, NumberofTipCells, NumberofBranchEvents, SproutLengths, NumberofBackTrackingLoops, BranchLengths, NumberofLargeSelfLoops, PerfusedNetwork, Perfused_ColumnAverage, TC_Matrix, TC_Solution, EC_Solution, LateralMovement, LengthofSprout, RatioofLateralMovement, SurvivingSproutsYesNo, ReachedTumorYesNo] = Pillay_2017_CA(c,Parameters)
% 2D CA Model from Pillay et al.
% Output:     Time, a T x 1 vector containing the time points of the
%               simulation
%             Network, X x Y x T array of cell positions at every unit of
%               time 
%             TC_Matrix, a matrix containing (row,col) location of all tip
%               cells
%             NumberofTipCells, a T x 1 vector containing the number
%               of tip cells at a time point
% Steady State Linear TAF Concentration: c(x,y,t) = x
%--------------------------------------------------------------------------
%----------- Solution/Storage Matrices ----------------------
Time = (0:Parameters.dt:Parameters.Tf)'; % Vector of times
Lt = length(Time);

% Initial Condition
% TC_Matrix = [100, 50; 95, 50]; % Store the index of R tip cells
TC_Matrix = uint16([(2:2:Parameters.N)', ones(length(2:2:Parameters.N),1)]);
% TC_Matrix = [(102:2:Parameters.N)', ones(length(102:2:Parameters.N), 1)];
% TC_Matrix = [150, 2; 50, 2];
Network = zeros(Parameters.N, Parameters.N, Lt, 'uint8'); % 3D array that stores cell occupancy at all time points (TC = EC = 1)
Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(:,1), TC_Matrix(:,2))) = 1; % Initial Condition

TC_Solution = uint8(Network);
EC_Solution = zeros(Parameters.N, Parameters.N, Lt, 'uint8');

NumberofSelfLoops = zeros(Lt, 1, 'uint16'); % Track the total number of self tip-sprout loops that would have been formed at every time step
NumberTipTipAnastomoses = NumberofSelfLoops; % Tracks the total number of tip-tip loop formation at every time step
NumberTipSproutAnastomoses = NumberofSelfLoops; % Tracks successful tip-sprout loop formation at every time step
NumberofBackTrackingLoops = NumberofSelfLoops; % Tracks the number of self loops formed by backtracking
NumberofLargeSelfLoops = 0; % Tracks the number of self loops that are 11 cells long or more
%% 5/24/19: Commented out Sprouts to put in the History matrix instead (code is the same, just different names and matrix sizes)
% Create a matrix containing the (x,y) location of TC Paths for every
% live tip cell
Sprouts = zeros(1000, 2, size(TC_Matrix, 1), 'uint16'); 
% Assumes here that no more than 1000 ECs will be created in the loop and
% that branching events create no more than 50 extra tips

% Create a matrix containing the (row, col) location of the last 10 (or last single) TC
% locations
Path_History = zeros(10, 2, size(TC_Matrix,1),'uint16');

% Create a vector that tracks how many branching points occur vs time
NumberofBranchEvents = NumberofSelfLoops;
% BranchLengths = [];
BranchLengths = zeros(Parameters.Q,Lt, 'uint16');
Branch_EC_Counter = zeros(0,0,'uint16');
BranchIndices = zeros(0,0,'uint16');



% Note: any cell > 0, Empty Site == 0.
EC_Counter = zeros(size(TC_Matrix,1),1, 'uint16'); % Will keep count of how long the sprout is
NumberofTipCells = zeros(Lt, 1, 'uint16'); % Keeps Track of Number of active tip cells
NumberofTipCells(1) = size(TC_Matrix,1);
index_store = zeros(0,0,'uint16');

% SproutLengths = []; % Will store the number of cells that formed a branch before an anastomosis event
SproutLengths = zeros(Parameters.Q, Lt, 'uint16'); % Stores both the original vessels AND branches
SproutLengths(1) = size(TC_Matrix,1);

PerfusedNetwork = zeros(Parameters.N,Parameters.N,Lt, 'logical');
PerfusedNetworkSubs = zeros(0,0,'uint16');

t = 1;
while t<Lt
    % Drop Stalk Cells where TCs used to be
    Network(:,:,t+1) = Network(:,:,t);
    PerfusedNetwork(:,:,t+1) = PerfusedNetwork(:,:,t);
    EC_Solution(:,:,t+1) = EC_Solution(:,:,t) + EC_Solution(:,:,t+1);
    TC_Solution(:,:,t+1) = TC_Solution(:,:,t);
    
    % Included as of 6/4/19
    SproutLengths(:,t+1) = SproutLengths(:,t+1)+SproutLengths(:,t);
    BranchLengths(:,t+1) = BranchLengths(:,t+1)+BranchLengths(:,t);
    
    if size(TC_Matrix,1) ~=0 % As long as you have tip cells,      
    % Choose R TCs at Random with Replacement
    R1 = randi(size(TC_Matrix,1),size(TC_Matrix,1),1);
    for ii = 1:length(R1)
        %% Step 1: Movement/Anastomosis
        % Draw a uniformly distributed random number
        r = rand;
        if r <= Parameters.P_m && ~isnan(R1(ii)) && any(TC_Matrix(R1(ii),1))
            % Record all sprouts
            % Add on another EC to the length of the sprout
            EC_Counter(R1(ii)) = EC_Counter(R1(ii)) + 1; % Storage index for the particular sprout you are at
            % Add on another EC to the length of the branch, if sprout
            % originated from a branching event
            Branch_EC_Counter(BranchIndices==R1(ii)) = Branch_EC_Counter(BranchIndices==R1(ii))+1;
            Sprouts(EC_Counter(R1(ii)), :, R1(ii)) = TC_Matrix(R1(ii),:);
            % 5/24/19: Added the following 2 lines for Path_History array
            % (store the last 10 paths, first row is newest, last row is
            % oldest)
            Path_History(2:end, :, R1(ii)) = Path_History(1:end-1, :, R1(ii));
            Path_History(1, :, R1(ii)) = TC_Matrix(R1(ii),:);
            % Delete the TC from the network, add it to the stalk cell network
            EC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) = EC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1))+1;
            TC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) = TC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1))-1;
            %% Movement
            % Draw a random number S in [0,1]
            S = rand;
            % Determine g_x, g_y for Probabilities
            if TC_Matrix(R1(ii), 2) ~= 1 && TC_Matrix(R1(ii),2)~= Parameters.N && TC_Matrix(R1(ii),1)~=1 && TC_Matrix(R1(ii),1)~=Parameters.N
            g_x = Parameters.k*(c(sub2ind(size(c), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2)+1))-c(sub2ind(size(c), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2)-1)));
            g_y = Parameters.k*(c(sub2ind(size(c), TC_Matrix(R1(ii),1)+1, TC_Matrix(R1(ii),2)))-c(sub2ind(size(c), TC_Matrix(R1(ii),1)-1, TC_Matrix(R1(ii),2))));
            else
                if TC_Matrix(R1(ii), 2) == 1
                    g_x = Parameters.k*(c(sub2ind(size(c), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2)+2))-c(sub2ind(size(c), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2))));
                elseif TC_Matrix(R1(ii),2)== Parameters.N
                    g_x = Parameters.k*(c(sub2ind(size(c), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2)))-c(sub2ind(size(c), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2)-2)));
                elseif TC_Matrix(R1(ii), 2) ~= 1 && TC_Matrix(R1(ii), 2) ~= Parameters.N
                    g_x = Parameters.k*(c(sub2ind(size(c), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2)+1))-c(sub2ind(size(c), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2)-1)));
                end
                if TC_Matrix(R1(ii),1)==1
                    g_y = Parameters.k*(c(sub2ind(size(c), TC_Matrix(R1(ii),1)+2, TC_Matrix(R1(ii),2)))-c(sub2ind(size(c), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2))));
                elseif TC_Matrix(R1(ii),1)==Parameters.N
                    g_y = Parameters.k*(c(sub2ind(size(c), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2)))-c(sub2ind(size(c), TC_Matrix(R1(ii),1)-2, TC_Matrix(R1(ii),2))));
                elseif TC_Matrix(R1(ii),1)~=1 && TC_Matrix(R1(ii),1)~=Parameters.N
                    g_y = Parameters.k*(c(sub2ind(size(c), TC_Matrix(R1(ii),1)+1, TC_Matrix(R1(ii),2)))-c(sub2ind(size(c), TC_Matrix(R1(ii),1)-1, TC_Matrix(R1(ii),2))));
                end
            end % if TC_Matrix
            assert(abs(g_x)<=1+2e-15, 'g_x is greater than 1. Adjust k.');
            assert(abs(g_y)<=1+2e-15, 'g_y is greater than 1. Adjust k.');
%             if abs(g_x)>1
%                 g_x = sign(g_x)*1;
%             end % if abs(g_x)>1
%             if abs(g_y)>1
%                 g_y = sign(g_y)*1;
%             end % if abs(g_y)>1
            if S < (1-g_x)/4
                if TC_Matrix(R1(ii),2)~=1
                % Move Left if you're not at boundary
                TC_Matrix(R1(ii), 2) = TC_Matrix(R1(ii),2)-1;
                assert(TC_Matrix(R1(ii),2)>0, 'Out of range subscript.');
                end
            elseif S>=(1-g_x)/4 && S<0.5
                if TC_Matrix(R1(ii),2)~=Parameters.N
                % Move Right if you're not at boundary
                TC_Matrix(R1(ii), 2) = TC_Matrix(R1(ii),2)+1;
                assert(TC_Matrix(R1(ii),2)<(Parameters.N+1), 'Out of range subscript.');
                end
            elseif S>=0.5 && S<0.5 + (1-g_y)/4
                if TC_Matrix(R1(ii),1)~=1 
                % Move Down if you're not at boundary
                TC_Matrix(R1(ii),1) = TC_Matrix(R1(ii),1)-1;
                assert(TC_Matrix(R1(ii),1)>0, 'Out of range subscript.');
                end
            else
                if TC_Matrix(R1(ii),1)~=Parameters.N
                % Move Up if you're not at boundary
                TC_Matrix(R1(ii),1) = TC_Matrix(R1(ii),1)+1;
                assert(TC_Matrix(R1(ii),1)<(Parameters.N+1), 'Out of range subscript.');
                end
            end % if S
            
            %% Make sure you update the length of the sprout
            SproutLengths(EC_Counter(R1(ii)), t+1) = SproutLengths(EC_Counter(R1(ii)), t+1)-1;
            SproutLengths(EC_Counter(R1(ii))+1, t+1) = SproutLengths(EC_Counter(R1(ii))+1, t+1)+1;
            % Branch_EC_Counter has already been incremented to account
            % for additional EC, but must add 1 in second line so that
            % TC is also accounted for
            BranchLengths(Branch_EC_Counter(BranchIndices==R1(ii)), t+1) = BranchLengths(Branch_EC_Counter(BranchIndices==R1(ii)), t+1)-1;
            BranchLengths(Branch_EC_Counter(BranchIndices==R1(ii))+1, t+1) = BranchLengths(Branch_EC_Counter(BranchIndices==R1(ii))+1, t+1)+1;
            %% Anastomosis    
            % Tip-to-Tip Anastomosis
            % NB: Tip-Tip Anastomosis will take priority over all other
            % anastomosis events
            AnastomosisEventOccured = false;
            Z = ismember(TC_Matrix, TC_Matrix(R1(ii),:), 'rows');
            if Parameters.a_n && sum(Z)>1
                % Record successful tip-tip loop formation
                NumberTipTipAnastomoses(t+1) = NumberTipTipAnastomoses(t+1)+1;
                
                % Find all tip cells that are located at the location of
                % TC_Matrix (In practice the number of tip cells colliding
                % can be greater than 2, since branching can deposit a TC
                % onto a location already occupied by a TC). For
                % simplicity, we have chosen to delete all TCs in this
                % location.
                SproutsIndices = find(Z==1);
                % Remove the TCs (potentially more than 1) that you collide with from the network
                TC_Solution(sub2ind(size(TC_Solution), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) = TC_Solution(sub2ind(size(TC_Solution), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1))-(sum(Z)-1);
                % Create an EC here
                EC_Solution(sub2ind(size(TC_Solution), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) = EC_Solution(sub2ind(size(TC_Solution), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1))+1;
                Sprouts(EC_Counter(R1(ii))+1,:,R1(ii)) = TC_Matrix(R1(ii),:);
                % This preserves indexing of TC Matrix so that you actually
                % drop a perfused cell where the TCs collided
                EC_Counter(R1(ii)) = EC_Counter(R1(ii))+1;
                % Collect the (row, col) indices of the perfused
                % tip/endothelial cells
                % If the number of TCs colliding is
                % greater than 2, then the following line of code captures
                % all sprouts that collided:
                for iq = 1:length(SproutsIndices)
                    PerfusedNetworkSubs = [PerfusedNetworkSubs; Sprouts(1:EC_Counter(SproutsIndices(iq)),:,SproutsIndices(iq))];
                end % for iq                
                
                % If one/both of these vessels was formed from a branch,
                [C, index_deletion] = intersect(BranchIndices,SproutsIndices);
                if ~isempty(C)
                    % Delete BranchIndices
                    Branch_EC_Counter(index_deletion) = NaN;
                    BranchIndices(index_deletion) = NaN;
                end % if ~isempty(C)

                % Delete all TCs that collide
                TC_Matrix(Z,:) = NaN; % Delete all the TCs that are in the location
                EC_Counter(Z) = NaN; % Must stop all sprouts from growing further also
                index_store = [index_store; SproutsIndices]; % This stores which matrices in the 3D array Sprouts that we will have to remove (so that we correctly reindex at the end of a time step)
                Sprouts(end,1,Z) = NaN;

                Path_History(end,1,Z) = NaN;
                R1(R1 == ismember(R1, SproutsIndices)) = NaN; % Delete any further instances of these TCs moving
                AnastomosisEventOccured = true;
            % Tip-to-Sprout Anastomosis
            elseif Parameters.a_e && ~AnastomosisEventOccured && EC_Solution(sub2ind(size(EC_Solution), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) && ~TC_Solution(sub2ind(size(TC_Solution), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1))
%                 if isempty(intersect(TC_Matrix(R1(ii),:), Sprouts(1:EC_Counter(R1(ii)),:,R1(ii)), 'rows'))
%                 if isempty(intersect(TC_Matrix(R1(ii),:), Path_History(:,:,R1(ii)), 'rows'))
                if ~ismember(TC_Matrix(R1(ii),:), Path_History(:,:,R1(ii)), 'rows')
                    % Record successful tip-sprout loop formation
                    NumberTipSproutAnastomoses(t+1) = NumberTipSproutAnastomoses(t+1) + 1;
                    
                    % Record successful self-loop formation that actually
                    % resulted in an anastomosis event (ie., if you're
                    % colliding more than 10 cells away from a TC)
                    if EC_Counter(R1(ii))>=11 && ~ismember(TC_Matrix(R1(ii),:), Sprouts(EC_Counter(R1(ii))-9:EC_Counter(R1(ii)), :, R1(ii)), 'rows')
                        NumberofLargeSelfLoops = NumberofLargeSelfLoops + 1;
                    end % if you successfully collide with your own sprout
                    
                    % Collect (row, col) data of Perfused Network
                        % Find which sprout you intersected with
                    if PerfusedNetwork(sub2ind(size(PerfusedNetwork), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) ~= 1
                        PerfusedNetworkSubs = [PerfusedNetworkSubs; Sprouts(1:EC_Counter(R1(ii)),:,R1(ii)); TC_Matrix(R1(ii),:)];
                        CheckPageNumber = zeros(size(Sprouts,3),1);
                        for iv = 1:size(Sprouts,3)
                            CheckPageNumber(iv) = ismember(TC_Matrix(R1(ii),:), Sprouts(:,:,iv), 'rows');
                        end
                            CheckPageNumber(R1(ii)) = 0;
                            V = find(CheckPageNumber==1, 1);
                            if ~isempty(EC_Counter(V))
                                % The above line checks to make sure you
                                % are colliding with a sprout that has
                                % already anastomosized - everything should
                                % be perfused already, so you don't have to
                                % add any new subscripts to
                                % PerfusedNetworkSubs.
                                
                                % If you collided with a sprout that has
                                % NOT anastomosized yet, the following
                                % lines of code will find the active sprout
                                % you collided with
                                % Find the last occurence of the EC that
                                % you collided with
                                [~, EC_Index] = intersect(flipud(Sprouts(:,:,V)), TC_Matrix(R1(ii),:), 'rows');
                                EC_Index = size(Sprouts,1)+1-EC_Index;
                                % Perfusion is everywhere up to and
                                % including that EC cell
                                PerfusedNetworkSubs = [PerfusedNetworkSubs; Sprouts(1:EC_Index, :, V)];
                            end
                    else
                        PerfusedNetworkSubs = [PerfusedNetworkSubs; Sprouts(1:EC_Counter(R1(ii)),:,R1(ii))];
                    end % if you are colliding with an already perfused network

                    if ismember(R1(ii), BranchIndices)
                        Branch_EC_Counter(BranchIndices==R1(ii)) = NaN;
                        BranchIndices(BranchIndices==R1(ii)) = NaN;
                    end % if ismember
                    % Delete the TC that collides with an EC
                    TC_Matrix(R1(ii),:) = [NaN, NaN];
                    EC_Counter(R1(ii)) = NaN;   
                    index_store = [index_store; R1(ii)];
                    Sprouts(end,1,R1(ii)) = NaN;
                    % 5/24/19: Following line is new code for Path_History
                    % stuff
                    Path_History(end,1,R1(ii)) = NaN;
                    R1(R1 == R1(ii)) = NaN; % Delete any further instances of this TC moving
                else
                    %% If self-loop would be formed, act as if anastomosis didn't happen
                    NumberofSelfLoops(t+1) = NumberofSelfLoops(t+1)+1;
                    
                    % If self-loop occured by backtracking, record it
%                     if EC_Counter(R1(ii))~=1 && ismember(TC_Matrix(R1(ii),:),Sprouts(EC_Counter(R1(ii))-1,:,R1(ii)),'rows')
                    % 5/24/19: Following line is new code for Path_History
                    % stuff
                    if EC_Counter(R1(ii))~=1 && ismember(TC_Matrix(R1(ii),:),Path_History(2,:,R1(ii)),'rows')
                        NumberofBackTrackingLoops(t+1) = NumberofBackTrackingLoops(t+1)+1;
                    end % if EC_Counter(R1(...
                    
                    % If a vessel is over 10 cells long and a self loop is
                    % formed that would be more than 10 cells in diameter,
                    % record it (this is only valid if Path_History isn't
                    % used).
%                     if EC_Counter(R1(ii))>=11 && ~ismember(TC_Matrix(R1(ii),:),Sprouts(EC_Counter(R1(ii))-9:EC_Counter(R1(ii)),:,R1(ii)),'rows')
%                         NumberofLargeSelfLoops = NumberofLargeSelfLoops+1;
%                     end % if EC_Counter(R1(ii))>=11                  
                    
                    % Put back the TC to the network and the TC solution
                    % arrays
                    Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) = Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1))+1;
                    TC_Solution(sub2ind(size(TC_Solution), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) = TC_Solution(sub2ind(size(TC_Solution), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) + 1;
                end % if TC doesn't collide with its own sprout, annihilate
            else
            % Add in the new positions of Tip Cells if anastomosis doesn't
            % occur
            Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) = Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1))+1;
            TC_Solution(sub2ind(size(TC_Solution), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) = TC_Solution(sub2ind(size(TC_Solution), TC_Matrix(R1(ii),1), TC_Matrix(R1(ii),2), t+1)) + 1;
            end % if Parameters.a_n
        end % if r
    end % for ii
    
    % These lines of code delete terminated sprouts from the 3D arrays, so
    % that when you delete all NaN variables in TC_Matrix and EC_Counter,
    % you don't accidentally change the indexing between arrays
%     if any(any(any(isnan(Sprouts))))
    if any(~all(TC_Matrix))
        % Reset indices for Branches, for calculating lengths
        for vv = 1:length(BranchIndices)
            BranchIndices(vv) = BranchIndices(vv)-sum(index_store<BranchIndices(vv));
        end % for vv
        for qq = 1:length(index_store)
            if index_store(qq) ~=1 && index_store(qq) ~= size(TC_Matrix,1)
                Sprouts = cat(3, Sprouts(:,:,1:index_store(qq)-1), Sprouts(:,:,index_store(qq)+1:end));
                Path_History = cat(3, Path_History(:,:,1:index_store(qq)-1), Path_History(:,:,index_store(qq)+1:end));
                index_store(index_store>index_store(qq)) = index_store(index_store>index_store(qq))-1;
            elseif index_store(qq) == 1
                Sprouts = Sprouts(:,:,index_store(qq)+1:end);
                Path_History = Path_History(:,:,index_store(qq)+1:end);
                index_store(index_store>index_store(qq)) = index_store(index_store>index_store(qq))-1;
            else
                Sprouts = Sprouts(:,:,1:index_store(qq)-1);
                Path_History = Path_History(:,:,1:index_store(qq)-1);
            end % if index_store(qq)
        end % for qq
        index_store = zeros(0,0,'uint16'); % Reset index storage vector for next round
    end % if any(any(...

    % Delete all annihilated TCs after movement steps are done
    EC_Counter(TC_Matrix(:,1)==0) = [];
    TC_Matrix(TC_Matrix(:,1)==0,:) = []; 
    Branch_EC_Counter(BranchIndices==0) = []; % Freeze any branch lengths that had an anastomosis event
    BranchIndices(BranchIndices==0) = [];
    assert(size(TC_Matrix,1)==sum(sum(TC_Solution(:,:,t+1))), 'The number of tip cells in the network being visualized is not equal to those being tracked. You have lost count of some tip cells.');
    end % if size(TC_Matrix, 1) ~=0
    %% Step 2: Branching
    if size(TC_Matrix,1)~=0
    % Choose R tip cells at random with replacement
    R2 = randi(size(TC_Matrix,1), size(TC_Matrix,1), 1);
    for jj = 1:length(R2)
        % Chose a uniformly distributed random number
        if ~isnan(R2(jj))
        r = rand;
        % Calculate Branching Probability at the Tip Cell's location
        
        % Original Model: Branching is Dependent on TAF Concentration
        P_b = Parameters.P_p*c(sub2ind(size(c), TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2)));

        % Simpler Model: Branching is Independent of TAF Concentration
%         P_b = Parameters.P_p;
        
        if r <= P_b
           %% Add on endothelial cell where you are, if you add EC after branching
           if Parameters.EC_Branching_YesNo
            Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) = Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1))+1;
            EC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) = EC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1))+1;
            SproutLengths(EC_Counter(R2(jj))+1,t+1) = SproutLengths(EC_Counter(R2(jj))+1,t+1)-1;
            SproutLengths(EC_Counter(R2(jj))+2,t+1) = SproutLengths(EC_Counter(R2(jj))+2,t+1)+1;
           end % if Parameters.EC_Branching_YesNo == 1
           
           %% Record position of Tip Cell
           Sprouts(EC_Counter(R2(jj))+1, :, R2(jj)) = TC_Matrix(R2(jj),:);
           
           %% Record a successful branching event
           NumberofBranchEvents(t+1) = NumberofBranchEvents(t+1) + 1;

           %% Delete tip Cell from Solution Arrays
           Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) = Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1))-1;
           TC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) = TC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1))-1;
           %% First Tip Cell branching
           if TC_Matrix(R2(jj),1) ~= 1
            % Branch in y-direction only
            TC_Matrix(R2(jj), 1) = TC_Matrix(R2(jj), 1) - 1;
            assert(TC_Matrix(R2(jj),1)>0, 'Out of range subscript.');
             % Adjust Solution Matrix
            Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) = Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) + 1;
            TC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) = TC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) + 1;
           else
                % No flux is reflecting - stay where you are
                % Adjust Solution Matrix
           Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) = Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) + 1;
           TC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) = TC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(R2(jj),1), TC_Matrix(R2(jj),2), t+1)) + 1;
           end

           %% New Tip Cell branching/creation
           if TC_Matrix(R2(jj),1)+2<=Parameters.N
               % Branch in y-direction only, top endothelial cell
               TC_Matrix = [TC_Matrix; [TC_Matrix(R2(jj),1) + 2, TC_Matrix(R2(jj),2)]];
               assert(TC_Matrix(R2(jj),1)<(Parameters.N+1), 'Out of range subscript.');
           else
               % Reflecting BC - Stay where you are
               TC_Matrix = [TC_Matrix; [TC_Matrix(R2(jj),1) + 1, TC_Matrix(R2(jj),2)]];
           end
           % Adjust Solution Matrix
           Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(end,1), TC_Matrix(end,2), t+1)) = Network(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(end,1), TC_Matrix(end,2), t+1))+1;
           TC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(end,1), TC_Matrix(end,2), t+1)) = TC_Solution(sub2ind([Parameters.N,Parameters.N,Lt], TC_Matrix(end,1), TC_Matrix(end,2), t+1))+1;

           % Define a new sprout as beginning where the TC cell is, while
           % keeping the complete history of the parent sprout                     
           Sprouts = cat(3, Sprouts, Sprouts(:,:,R2(jj)));
           Path_History = cat(3, Path_History, Path_History(:,:,R2(jj)));
           
           % Include if you add EC after branching
           if Parameters.EC_Branching_YesNo
              Path_History(3:end, :, R2(jj)) = Path_History(1:end-2, :, R2(jj));
              Path_History(2, :, R2(jj)) = Sprouts(EC_Counter(R2(jj))+1, :, R2(jj));
              Path_History(3:end, :, end) = Path_History(3:end, :, R2(jj));
              Path_History(2, :, end) = Path_History(2, :, R2(jj));
              EC_Counter(R2(jj)) = EC_Counter(R2(jj))+1;
           end % if you add EC after Branching
           
           Path_History(1, :, R2(jj)) = TC_Matrix(R2(jj), :);
           Path_History(1, :, end) = TC_Matrix(end, :);
           
           % Add on new Sprout to the EC_Counter
           EC_Counter = [EC_Counter; EC_Counter(R2(jj))];
           % Adjust BranchLengths and SproutLenghts to include new branches
           BranchLengths(1, t+1) = BranchLengths(1, t+1) + 2;           
           SproutLengths(EC_Counter(R2(jj))+1, t+1) = SproutLengths(EC_Counter(R2(jj))+1, t+1)+1;

           % Store the branch indices, so that you can record their length
           % later
           if ~ismember(R2(jj), BranchIndices)
             BranchIndices = [BranchIndices; R2(jj); size(Sprouts,3)];
             Branch_EC_Counter = [Branch_EC_Counter; 0; 0];
           else
             BranchIndices = [BranchIndices; size(Sprouts,3)];
             Branch_EC_Counter(BranchIndices==R2(jj)) = 0;
             Branch_EC_Counter = [Branch_EC_Counter; 0];
           end % if ~ismember
        R2(R2==R2(jj)) = NaN; % Delete another instance of branching from this tip cell
        assert(size(TC_Matrix,1)==sum(sum(TC_Solution(:,:,t+1))), 'The number of tip cells in the network being visualized is not equal to those being tracked. You have lost count of some tip cells.');
        end % if R < P_b
        end % if ~isnan(R2(jj))
    end % for jj
    end % if size(TC_Matrix,1)~=0
    
    if size(PerfusedNetworkSubs,1)~=0
    % Use (row, col) data to fill in Perfused Network structure
    PerfusedNetwork(sub2ind(size(PerfusedNetwork), PerfusedNetworkSubs(:,1), PerfusedNetworkSubs(:,2), (t+1)*ones(size(PerfusedNetworkSubs,1),1))) = 1;
    PerfusedNetworkSubs = [];
    end % if size(PerfusedNetworkSubs...)
    %% Step 3: Update Time
    NumberofTipCells(t+1) = size(TC_Matrix,1);
    if sum(Network(:,end,t))~=0 || size(TC_Matrix,1)==0
        Network(:,:,t+2:end) = repmat(Network(:,:,t+1), 1, 1, length(t+2:Lt));
        PerfusedNetwork(:,:,t+2:end) = repmat(PerfusedNetwork(:,:,t+1), 1, 1, length(t+2:Lt));
        TC_Solution(:,:,t+2:end) = repmat(TC_Solution(:,:,t+1), 1,1, length(t+2:Lt));
        EC_Solution(:,:,t+2:end) = repmat(EC_Solution(:,:,t+1), 1,1, length(t+2:Lt));
        NumberofTipCells(t+2:end) = NumberofTipCells(t+1);
        SproutLengths(:, t+2:end) = repmat(SproutLengths(:, t+1), 1, 1, length(t+2:Lt));
        BranchLengths(:, t+2:end) = repmat(BranchLengths(:, t+1), 1, 1, length(t+2:Lt));
        break
    end % if TC's have reached the tumor
    t = t + 1;
end % while t

%% Collect Final Data
% For all surviving sprouts, record their total length in cells as
% well as the total cell lengths that they moved in the y-dimension, which
% if the probability of moving left is 0, means that it is the total cell
% length minus their final x-position
LateralMovement = zeros(size(Sprouts,3),1);
LengthofSprout = LateralMovement;
RatioofLateralMovement = LateralMovement;
for vi = 1:size(Sprouts, 3)
    % Records the number of lateral steps you have, while the index tells
    % you how many cell lengths in total the sprout is
    LateralMovement(vi) = EC_Counter(vi) + 1 - TC_Matrix(vi, 2);
    LengthofSprout(vi) = EC_Counter(vi)+1;
    RatioofLateralMovement(vi) = LateralMovement(vi)/LengthofSprout(vi);
end % for i

% If at least one TC survives at the end of the simulation, consider that
% as a successful angiogenesis event
if size(TC_Matrix,1)
    SurvivingSproutsYesNo = true;
else 
    SurvivingSproutsYesNo = false;
end
% If you reached the tumor, record it
if sum(Network(:,end,end))
    ReachedTumorYesNo = true;
else
    ReachedTumorYesNo = false;
end
% Take Column Averages
TC_ColumnAverage = squeeze(mean(double(TC_Solution), 1));
EC_ColumnAverage = squeeze(mean(double(EC_Solution), 1));
Perfused_ColumnAverage = squeeze(mean(double(PerfusedNetwork), 1));
end % function Pillay_2017_CA.m 

function PillayCA_WriteVideos(filename, Networks, PerfusedNetworks, BranchLengths, SproutLengths, Time, Average2DStalkCellNetwork, Average2DTipCellNetwork, TC_Solution, Parameters)
v = VideoWriter([filename,'_ExampleNetwork-with-VesselLengths']);
v.Quality = 100;
v.FrameRate = 5; % Want about 5 time steps per second, although you can change this

g = VideoWriter([filename,'_ExampleNetwork-with-BranchLengths']);
g.Quality = 100;
g.FrameRate = 5;

H = VideoWriter([filename,'_ExampleNetworkOnly']);
H.Quality = 100;
H.FrameRate = 5;

K = VideoWriter([filename,'_PerfusedNetworkOnly']);
K.Quality = 100;
K.FrameRate = 5;


l = VideoWriter([filename,'_Average2DTipCellNetwork']);
l.Quality = 100;
l.FrameRate = 5;


m = VideoWriter([filename,'_Average2DStalkCellNetwork']);
m.Quality = 100;
m.FrameRate = 5;

n = VideoWriter([filename,'_AverageTCWithExampleNetwork']);
n.Quality = 100;
n.FrameRate = 5;

r = VideoWriter([filename,'_AverageECWithExampleNetwork']);
r.Quality = 100;
r.FrameRate = 5;

u = VideoWriter([filename,'_NetworkOverlayPerfusion+TC']);
u.Quality = 100;
u.FrameRate = 5;

fig = figure;
    ax1 = subplot(2,1,1);
    graph1 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', NaN);
    xlabel('x');
    ylabel('y');
    axis tight
    caxis([0 max(max(max(Networks)))]);
    colormap(jet);
    clabel = colorbar;
    set(get(clabel, 'Title'), 'String', 'Number of Cells');
%     axis tight
    set(ax1,'nextplot','replacechildren', 'Visible','on');
    title(['Network @ t = ',num2str(Time(1),'%7.5f'),' s']);

ax2 = subplot(2,1,2);
    graph2 = bar(NaN,NaN);
    xlabel('Length of Sprouts, in Number of Cells');
    ylabel('Number of Sprouts');
    xlim([0, Parameters.Q]);
    ylim([0, max(SproutLengths(2:end, end, end))]);
    set(graph2, 'XData', 1:Parameters.Q);
    set(ax2, 'nextplot', 'replacechildren', 'visible', 'on');
    
    
fig2 = figure;
    ax4 = subplot(2,1,1);
    graph4 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', NaN);
    xlabel('x');
    ylabel('y');
    axis tight
    caxis([0 max(max(max(Networks)))]);
    colormap(jet);
    clabel = colorbar;
    set(get(clabel, 'Title'), 'String', 'Number of Cells');
    set(ax4, 'nextplot','replacechildren', 'Visible','on');
    title(['Network @ t = ',num2str(Time(1),'%7.5f'),' s']);

    ax5 = subplot(2,1,2);
    graph5 = bar(NaN, NaN);
    xlabel('Length of Branches, in Number of Cells')
    ylabel('Number of Branches')
    xlim([0, 100]);
    ylim([0, 10]);
    set(graph5, 'XData', 1:Parameters.Q);
    set(ax5, 'nextplot', 'replacechildren', 'visible', 'on');
    
fig3 = figure;
    graph6 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata',NaN);
    xlabel('x')
    ylabel('y');
    axis tight
    caxis([0 max(max(max(Networks)))]);
    colormap(jet);
    clabel = colorbar;
    set(get(clabel, 'Title'), 'String', 'Number of Cells');
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    title(['Network @ t = ',num2str(Time(1),'%7.5f'),' s']);
    ax6 = gca;
    
fig4 = figure;
    graph7 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', NaN);
    xlabel('x')
    ylabel('y');
    axis tight
    caxis([0 1]);
    Map1 = [1 1 1; 1 0 0];
    colormap(Map1);
    clabel = colorbar;
    set(get(clabel, 'Title'), 'String', 'Perfusion Status');
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    title(['Perfused Network @ t = ',num2str(Time(1),'%7.5f'),' s']);
    ax7 = gca;   
    
fig5 = figure;
    graph8 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata',NaN);
    xlabel('x')
    ylabel('y');
    axis tight
    caxis([0 max(max(Average2DTipCellNetwork(:,:,end)))]);
    colormap(jet);
    clabel = colorbar;
    set(get(clabel, 'Title'), 'String', 'Number of Cells');
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    title(['Average TC Network @ t = ',num2str(Time(1),'%7.5f'),' s']);
    ax8 = gca;
    
fig6 = figure;
    graph9 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', NaN);
    xlabel('x')
    ylabel('y');
    axis tight
    caxis([0 max(max(max(Average2DStalkCellNetwork)))]);
    colormap(jet);
    clabel = colorbar;
    set(get(clabel, 'Title'), 'String', 'Number of Cells');
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    title(['Average EC Network @ t = ',num2str(Time(1),'%7.5f'),' s']);
    ax9 = gca;
    
fig7 = figure;
    ax10 = subplot(2,1,1);
    graph10 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', NaN);
    xlabel('x');
    ylabel('y');
    axis tight
    caxis([0 max(max(Average2DTipCellNetwork(:,:,end)))]);    
    colormap(jet);
    clabel = colorbar;
    set(get(clabel, 'Title'), 'String', 'Number of Cells');
%     axis tight
    set(ax10,'nextplot','replacechildren', 'Visible','on');
    title(['Average TC Network @ t = ',num2str(Time(1),'%7.5f'),' s']);
    
    ax11 = subplot(2,1,2);
    graph11 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', NaN);
    xlabel('x');
    ylabel('y');
    axis tight
    caxis([0 max(max(max(Networks)))]);
    colormap(jet);
    clabel = colorbar;
    set(get(clabel, 'Title'), 'String', 'Number of Cells');
%     axis tight
    set(ax11,'nextplot','replacechildren', 'Visible','on');
    title(['Network @ t = ',num2str(Time(1),'%7.5f'),' s']);
    
fig8 = figure;
    ax12 = subplot(2,1,1);
    graph12 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', NaN);
    xlabel('x');
    ylabel('y');
    axis tight
    caxis([0 max(max(Average2DStalkCellNetwork(:,:,end)))]);    
    colormap(jet);
    clabel = colorbar;
    set(get(clabel, 'Title'), 'String', 'Number of Cells');
%     axis tight
    set(ax12,'nextplot','replacechildren', 'Visible','on');
    title(['Average EC Network @ t = ',num2str(Time(1),'%7.5f'),' s']);
    
    ax13 = subplot(2,1,2);
    graph13 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', NaN);
    xlabel('x');
    ylabel('y');
    axis tight
    caxis([0 max(max(max(Networks)))]);
    colormap(jet);
    clabel = colorbar;
    set(get(clabel, 'Title'), 'String', 'Number of Cells');
%     axis tight
    set(ax13,'nextplot','replacechildren', 'Visible','on');
    title(['Network @ t = ',num2str(Time(1),'%7.5f'),' s']);

fig11 = figure;
    graph16 = imagesc('xdata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'ydata', 0:Parameters.h:Parameters.h*(Parameters.N-1), 'cdata', NaN);
    xlabel('x')
    ylabel('y');
    axis tight
    caxis([0 3]);
    % Colormap in RGB: [No Cells, Unperfused Cells, Tip Cells, Perfused Cells] 
    colormap([1 1 1; 218/255, 112/255, 214/255; 0 0 1; 1 0 0]);
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    title(['Example Network w/Overlay @ t = ',num2str(Time(1),'%7.5f'),' s']);
    ax16 = gca;
    
open(v);
open(g);
open(H);
open(K);
open(l);
open(m);
open(n);
open(r);
open(u);
for frame_counter = 1:length(Time)
    set(graph1, 'CData', Networks(:,:,frame_counter));
    set(graph2, 'YData', SproutLengths(:, frame_counter, end));
    set(get(ax1, 'title'), 'String', ['Network @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
%     set(graph3, 'YData', BranchLengths(:,i,end));
    set(graph4, 'CData', Networks(:,:,frame_counter));
    set(graph5, 'YData', BranchLengths(:,frame_counter,end));
    set(get(ax4, 'title'), 'String', ['Network @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
    set(graph6, 'CData', Networks(:,:,frame_counter));
    set(get(ax6, 'title'), 'String', ['Network @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
    set(graph7, 'CData', PerfusedNetworks(:,:,frame_counter));
    set(get(ax7, 'title'), 'String', ['Perfused Network @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
    set(graph8, 'CData', Average2DTipCellNetwork(:,:,frame_counter));
    set(get(ax8, 'title'), 'String', ['Average TC Network @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
    set(graph9, 'CData', Average2DStalkCellNetwork(:,:,frame_counter));
    set(get(ax9, 'title'), 'String', ['Average EC Network @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
    set(graph10, 'CData', Average2DTipCellNetwork(:,:,frame_counter));
    set(get(ax10, 'title'), 'String', ['Average TC Network @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
    set(graph11, 'CData', Networks(:,:,frame_counter));
    set(get(ax11, 'title'), 'String', ['Network @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
    set(graph12, 'CData', Average2DStalkCellNetwork(:,:,frame_counter));
    set(get(ax12, 'title'), 'String', ['Average EC Network @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
    set(graph13, 'CData', Networks(:,:,frame_counter));
    set(get(ax13, 'title'), 'String', ['Network @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
    set(graph16, 'cdata', (Networks(:,:,frame_counter)>0) + 2*(PerfusedNetworks(:,:,frame_counter)>0) + logical(TC_Solution(:,:,frame_counter)));
    set(get(ax16, 'title'), 'string', ['Example Network w/Overlay @ t = ',num2str(Time(frame_counter),'%7.5f'),' s']);
    
    frame = getframe(fig);
    frame2 = getframe(fig2);
    frame3 = getframe(fig3);
    frame4 = getframe(fig4);
    frame5 = getframe(fig5);
    frame6 = getframe(fig6);
    frame7 = getframe(fig7);
    frame8 = getframe(fig8);
    frame11 = getframe(fig11);
    
    writeVideo(v, frame);
    writeVideo(g, frame2);
    writeVideo(H, frame3);
    writeVideo(K, frame4);
    writeVideo(l, frame5);
    writeVideo(m, frame6);
    writeVideo(n, frame7);
    writeVideo(r, frame8);
    writeVideo(u, frame11);
end % for i
close(v);
close(g);
close(H);
close(K);
close(l);
close(m);
close(n);
close(r);
close(u);
end % function PillayCA_WriteVideos.m