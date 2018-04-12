% The following code is inspired by "importTrackMates" and has
% been motified by Chanya Godzich and Anja Kunze to classify
% vesicle motion into five categories, to calculate the total
% length traveled in a specific time window and to calculate the
% caged dimater after ... et al.
% It includes:
%	- XML file import from TrackMate
% 	- Connecting time points
% 	- Filtering short trajectories
%	- Filtering to long trajectories
% 	- Set origin for each trajectory = 0,0
% 	- Calculates the total Trajectory Length (L), Average
% 		Velocity (V) alonge the trajectory, the MSD (MSD), Caged
%		diameter (CD)
% 	- Outputs data into .csv files

% Appended as the xml file name followed by a letter:
%   T (trajectories), L (length), V (velocity), MSD (MSD)

%% CLEAR ALL
clear all;
close all;
clc;

%% IMPORT - XML file import from TrackMate, Copy in path and it can start

file = 'C:\Users\ucheosahor\Desktop\New folder\6DIV-Chi-Norm-nM-DiD_Tracks3.xml';
[pathstr,name,ext] = fileparts(file);
Name = [pathstr, '\',name];
[tracksorigin, md] = importTrackMateTracks(file);

%% GLOBAL PARAMETER SETTINGS

ma2 = msdanalyzer(3, 'microns','seconds');      %ma2 is used for generating the (0,0) centered trajectory plot

CUT = 20;                                      %Maximal time length for evaluation in sec
TMIN = 0.9*CUT;                                 %Minimum time length for evaluation in sec
tmin = round(TMIN/md.frameInterval);                  %Minimum time length for evaluation in frames
tmax = round(CUT/md.frameInterval);                     %Maximum time length for evaluation in frames

%% CREATE CENTERED TRACKS

trackscenter = tracksorigin;                % "tracks" contains original track data, "track2" is used for centering
totalTracks = numel(trackscenter);          % "totalTracks" Counts total number of tracks

for i = 1:totalTracks                     % For all tracks interative extract the following values
    currentTrack = trackscenter{i};
    tValue = currentTrack(1,1);
    xValue = currentTrack(1,2);
    yValue = currentTrack(1,3);
    totalPointsWithinTrack = size(currentTrack,1);
        
    currentTrack(:,1) = ((currentTrack(:,1) - tValue)*md.frameInterval);
    currentTrack(:,2) = currentTrack(:,2) - xValue;
    currentTrack(:,3) = currentTrack(:,3) - yValue;
    currentTrack(:,4) = i;
    
    trackscenter{i} = currentTrack;          % New centered, tresholded data set includes short and long trajectories
end


%Ttrackscenter = transpose(trackscenter);
%% SELECT CENTERED TRACKS
tracksCUT = cell(totalTracks,1);
for i = 1:totalTracks                     % Sets all short trajectories to zero, standard: minimum length tmin = 0.8*CUT
    currentTrack2 = trackscenter{i};
    totalPointsWithinTrack2 = size(currentTrack2,1);
    total = (totalPointsWithinTrack2 <= tmax);
    timestep = totalPointsWithinTrack2*md.frameInterval;
    
    currentTrack3 = currentTrack2;
    
    if timestep < CUT
        if timestep < TMIN
            currentTrack3(:,2:3)=0;
        end
    else
        currentTrack3((tmin+1):end,:) = [];
    end
    
    %BULLLLSHITTTT
    %if (timestep < CUT)
    %
    %
    %             for j = 1 : totalPointsWithinTrack2
    %                if timestep >= TMIN
    %                 currentTrack3(j,1) = currentTrack2(j,1);
    %                 currentTrack3(j,2) = currentTrack2(j,2);
    %                 currentTrack3(j,3) = currentTrack2(j,3);
    %                 currentTrack3(j,4) = currentTrack2(j,4);
    %                     else
    %                     currentTrack3(j,1) = currentTrack2(j,1);
    %                     currentTrack3(j,2) = 0;
    %                     currentTrack3(j,3) = 0;
    %                     currentTrack3(j,4) = currentTrack2(j,4);
    %                end
    %             end
    %         else
    %             n = 1;
    %             while n < limit + 1
    %                 currentTrack3(n,1) = currentTrack2(n,1);
    %                 currentTrack3(n,2) = currentTrack2(n,2);
    %                 currentTrack3(n,3) = currentTrack2(n,3);
    %                 currentTrack3(n,4) = currentTrack2(n,4);
    %                 n = n + 1;
    %             end
    %         end
    %
    tracksCUT{i} = currentTrack3;               % Sets all short trajectories to zero, standard: minimum length tmin = 0.8*CUT
end

%TtracksCUT = transpose(tracksCUT);


%% PLOT CENTERED TRACKS

clf
ma2 = ma2.addAll(tracksCUT);
ma2.plotTracks;
ma2.labelPlotTracks;
axis([-4,4,-4,4]);

%% CALCULATE MSD

M = trackscenter;
% t = M(:,1);
% X = M(:,2);
% Y = M(:,3);

figure('Name', 'MSD_DZ');

MSDcell = [];

for traj_N = 1:length(M)                     % loop over all trajectories
    ts = M{traj_N}(:,1);                     % time points
    x = M{traj_N}(:,2);                      % x and y
    y = M{traj_N}(:,3);
    N = length(ts);
    
    
    %skip short trajectories
    if  N < tmin
        continue;
    end
    
    MSD = zeros(round(N/2),1);
    for tau = 1:round(N/2)
        for k = 1:(N-tau)
            MSD(tau) = MSD(tau) + ( (x(k+tau)-x(k)).^2 + (y(k+tau)-y(k)).^2 );
        end
        MSD(tau) = MSD(tau) / (N-tau);
        MSDcell(end+1,1:3) = [traj_N, tau*md.frameInterval, MSD(tau)]; %#ok<SAGROW>
        
    end
    
    plot( (1:round(N/2))*ts(2), MSD ); hold on;
    
    
end



%% CALCULATE TOTAL LENGTH AND SUM UP TIME POINTS
tracks3 = trackscenter;
%totalTracks = totalTracks(1,1);
%sizeLTracks = zeros(totalTracks,1);     %will output a matrix of total change in time points: dt

%currentTrack2 = trackscenter{i};
%   totalPointsWithinTrack2 = numel(currentTrack2)/4;
%  total = le(totalPointsWithinTrack2,max);
% currentTrack3 = zeros(total(1,1), 4);
%timestep = totalPointsWithinTrack2*md.frameInterval;

L = zeros(totalTracks,1);
TrackTime = L;

for i = 1:totalTracks                  %Iteration up to total number of tracks
    currentTrackL = tracks3{i};
    
    LTracks = currentTrackL(2:end,2:3) - currentTrackL(1:(end-1),2:3);
    sumofsquares = LTracks(:,1).^2 + LTracks(:,2).^2;
    L(i) = sum( sqrt(sumofsquares(:)) ); 
    
    TrackTime(i) = currentTrackL(end,1);
end

%% CALCULATES MEAN VELOCITY
V = L./(TrackTime);

%% CREATE CAGED RADIUS GRAPH RELATIVE TO FIRST POSITION - CELL ARRAY
tracks5 = trackscenter;
totalTracks5 = numel(tracks5);
CagedTracks2 = cell(totalTracks5,1);
for i = 1:totalTracks5
    currentTrack5 = tracks5{i};
    currentTrack6 = currentTrack5;
    totalPointsWithinTrack5 = size(currentTrack5,1);

    sumofsquares = currentTrack5(:,2).^2 + currentTrack5(:,3).^2;
        
    currentTrack6(:,2) = currentTrack5(:,1);
    currentTrack6(:,3) = sqrt( sumofsquares );
    CagedTracks2{i} = currentTrack6;
    
end



figure;
hold on;
avgCD = zeros(1,numel(CagedTracks2));
for i = 1 :numel (CagedTracks2)
    avgCD(i) = mean((CagedTracks2{i}(:,3))); %Column 3 is Caged Diameter in um
    
end

scatter( L , avgCD, 100, 'b'); %AverageCD vs Length
axis([0,30,0,30]);
ylabel 'Average Caged Diameter  (um)';
xlabel 'Total Length (um)';

figure;
scatter(V, avgCD, 100, 'b') %AverageCD vs Velocity
axis([0,5,0,5]);
ylabel 'Average Caged Diameter  (um)';
xlabel 'Average Velocity (um)';


Cat{1} = [.6, 0, 2, 0];      %DOCKED: CD < 0.6 AND L < 2              (COL 3)
Cat{2} = [ 2,.6, 3, 0];      %CAGED: 0.6 < CD < 2 AND L < 3       	(COL 4)
Cat{3} = [.6, 0, inf, 2];      %DOCKED TRANSPORT: CD < 0.6 AND L > 2 	(COL 5)
Cat{4} = [ 2,.6, inf, 3];      %CAGED TRANSPORT: 0.6 < CD < 2 AND L > 3	(COL 6)
Cat{5} = [ inf, 2, inf, 0];      %DIRECTED TRANSPORT CD > 2               (COL 7)

avgCDT = transpose(avgCD);

HOLD = avgCDT*0;

CD = [avgCDT, L, HOLD, HOLD, HOLD, HOLD, HOLD, HOLD];
NumCD = size(CD,1);

for k=1:numel(Cat)
    CD( avgCDT > Cat{k}(2) & avgCDT < Cat{k}(1) & L > Cat{k}(4) & L < Cat{k}(3), 2+k) = 1;
end

CD( :,end ) = 1 - sum(CD(:,3:end), 2);


% 
% for i = 1: NumCD
%     if CD(i,1) > 0
%         if CD(i,1) > CatE(1,2)              %DIRECTED TRANSPORT CD > 2       (COL 7)
%             CD(i,7) = 1;
%             
%             
%         elseif CD(i,1) < CatB(1,1)          % CD < 5
%             if CD(i,2) < CatA(1,3)          % L < 3
%                 if CD(i,1) < CatA(1,1)      % CD < .75
%                     CD(i,3) = 1;            % DOCKED: CD < 0.5 AND L < 5  (COL 3)
%                 end
%             end
%             if CD(i,2) < CatB(1,3)          % L < 5
%                 if CD(i,2) > CatB(1,4)      % L > 3
%                     CD(i,4) = 1;            %CAGED: CD < 5 AND 3 < L < 5	(COL 4)
%                 end
%             end
%             if CD(i,2) > CatC(1,4) % L > 5
%                 if CD(i,1) < CatC(1,1)  % CD < 2
%                     CD(i,5) = 1;    %DOCKED TRANSPORT: CD < 2 AND L > 5 (COL 5)
%                     
%                 elseif CD(i,1) > CatD(1,2)  % CD > 2
%                     CD(i,6) = 1;    %CAGED TRANSPORT: 2 < CD < 5 AND L > 5	(COL 6)
%                 end
%             end
%             
%         else
%             CD(i,8) = 1;
%             
%         end
%     else
%         CD(i,8) = 1;
%     end
% end
CDCatCount = sum(CD);

%% OUTPUTS DATA

nameCutT = strcat(Name, 'CutTracks.csv');          %Designates a name to .csv file to describe contents as CutTracks
csvwrite(nameCutT, tracksCUT);                    %Outputs Tracks data to .csv file

nameCutC = strcat(Name, 'CenterTracks.csv');       %Designates a name to .csv file to describe contents as CutTracks
csvwrite(nameCutC, trackscenter);                  %Outputs Tracks data to .csv file

%nameCutL = strcat(Name, 'CutL.csv');
%csvwrite(nameCutL, L);                            %Outputs CUT Length data to .csv file

%nameCutV = strcat(Name, 'CutV.csv');
%csvwrite(nameCutV, V);                            %Outputs CUT Velocity data to .csv file

nameCutCD = strcat(Name, 'CutCD-L.csv');
csvwrite(nameCutCD, CD);

nameCDCatCount = strcat(Name, 'CDCatCount.csv');
csvwrite(nameCDCatCount, CDCatCount);

nameCutMSD = strcat(Name, 'CutMSD.csv');
csvwrite(nameCutMSD, MSDcell);                     %Outputs CUT MSD data to .csv file

