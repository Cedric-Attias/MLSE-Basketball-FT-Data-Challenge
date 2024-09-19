%Specify the folder containing the trajectory files (x,y,z data)
folder_path = '';  % Replace with your folder path*************************

% Get a list of all .mat files in the folder
files = dir(fullfile(folder_path, '*.mat'));

% Initialize a figure for plotting
figure;

%Hoop
radius = 0.75;
center = [14.25,0,10];
theta = linspace(0, 2*pi, 100);
xc = radius * cos(theta) + center(1);
yc = radius * sin(theta) + center(2);
zc = center(3) * ones(size(theta)); % z-coordinate is constant
plot3(xc, yc, zc, '-', 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980]);
hold on

% FT line
x_line = [0, 0]; % x-coordinates for the line
y_line = [-6, 6]; % y-coordinates for the line
z_line = [0, 0]; % z-coordinates for the line

% Backboard
xB = [15 15 15 15];  % X-coordinates of the corners
yB = [-3 -3 3 3];  % Y-coordinates of the corners
zB = [13 9 9 13];  % Z-coordinates (all at the same height)
fill3(xB, yB, zB, 'black');  % 'r' specifies red color for the rectangle
line(x_line, y_line, z_line, 'Color', 'red', 'LineWidth', 5);
xlabel('x (ft)');
ylabel('y (ft)');
zlabel('z (ft)');
grid on
set(gca, 'YDir', 'reverse'); % Flip the y-axis
hold on;

missCount = 0;
madeCount = 0;

% Loop through each file
for i = 1:length(files)
    % Construct the full file path
    file_path = fullfile(folder_path, files(i).name);
    
    % Load the data from the file
    data = load(file_path);

    x = data.x_traj;
    y = data.y_traj;
    z = data.z_traj;

    if strcmp(data.result, 'missed')
        % Trajectory
        p = plot3(x, y, z);
        p.MarkerSize = 20;
        p.LineStyle = "--";
        p.Color = 'red';
        missCount = missCount+1;
    else
        % Trajectory
        p = plot3(x, y, z);
        p.MarkerSize = 20;
        p.LineStyle = "--";
        p.Color = 'green';
        madeCount = madeCount+1;
    end
end

FTPerc = (madeCount/(madeCount+missCount)*100);
totalShots = madeCount+missCount;

str = sprintf('FT%%:  %.1f (%i / %i)', FTPerc,madeCount, totalShots); % Example string
txt1 = text(15, -2, 15, str,'FontSize', 20, 'Color', 'black');
hold off;