% Define the folder containing the files
folder_path = 'C:\Users\Cedric Attias\OneDrive - University of Waterloo\MLSE Data Challenge\FT Data\basketballFreethrow_data_P0001';  % Replace with your folder path

% Get a list of all files in the folder
fileList = dir(fullfile(folder_path, '*.json'));  % Replace '*.txt' with the file extension you're interested in

% Loop through each file in the folder
for i = 1:numel(fileList)
    % Get the full path of the current file
    filePath = fullfile(folder_path, fileList(i).name);
    
    % Display the file being processed (optional)
    fprintf('Processing file: %s\n', fileList(i).name);
    
    % Call your existing script or function to process the file
    % If your script processes a single file, you can modify it to take
    % the file path as an input argument

    FTTrajectory(filePath);

    clear FTTrajectory;
    
    % Your processing logic here
    % For example, if you had some operations to perform on `data`
    
    % Save or display results as needed
end


function FTTrajectory(InputFile)
    %% Machine Learning Model for Predicting Free-Throw Success using a Physics-based Basektball Trajectory Model
    % Assumption in the model: 
    % Assuming there is back spin, y-axis (Magnus force on x)
    % Only forces present are that of gravity (-Z) and air resitance against the
    % ball (-X)
    % NBA court dimensions  
    % Ball is a smooth sphere
    
    %% STEP 1: Read the JSON data files and process the data to only include known ball coordinate indices
    input = InputFile;
    % testFolder = 'C:\Users\Cedric Attias\OneDrive - University of Waterloo\MLSE Data Challenge\FT Data\TEST';
    
    jsonFile = fileread(input);
    json_str = jsondecode(jsonFile);

    ballData = struct();
    trimmedBallData = struct();
    
    %Extract pertinent information for this model
    for j = 1:length(json_str.tracking)
        ballData.position{j} = json_str.tracking(j).data.ball;
        ballData.result{j} = json_str.result;
        ballData.landingX{j} = json_str.landing_x;
        ballData.landingY{j} = json_str.landing_y;
        ballData.entryAngle{j} = json_str.entry_angle;
        ballData.frames{j} = json_str.tracking(j).frame;
        ballData.time{j} = json_str.tracking(j).time;
    end
    
    %Crop the data based on first and last known ball position
    nanIndex = zeros(1,length(ballData.position));
    pstn = ballData.position;
    
    for k = 1:length(ballData.position)
        if ~isnan(pstn{1,k}(1,1))
            nanIndex(k) = k;
        end
    end
    
    % Ball release occurs at peak elbow extensions (or height)
    newElbow = 0;
    maxElbow = 0;
    maxElbowFrame =0;
    for kk = 1:length(ballData.position)
        newElbow = (json_str.tracking(kk).data.player.R_ELBOW(end,1));
           if newElbow >= maxElbow
               maxElbow =newElbow;
               maxElbowFrame = kk;
           end
    end
    
    nanIndex(nanIndex == 0) = NaN;
    firstNonNanIndex = min(nanIndex);
    lastNonNanIndex = max(nanIndex);
    
    % Store all the known positions and asssociated data for the ball
    trimmedBallData.position = ballData.position(1,firstNonNanIndex:lastNonNanIndex);
    trimmedBallData.result = ballData.result(1);
    trimmedBallData.landingX = ballData.landingX(1);
    trimmedBallData.landingY = ballData.landingY(1);
    trimmedBallData.entryAngle = ballData.entryAngle(1);
    trimmedBallData.frames = ballData.frames(1,firstNonNanIndex:lastNonNanIndex);
    trimmedBallData.time = ballData.time(1,firstNonNanIndex:lastNonNanIndex);
    
    dT = (trimmedBallData.time{1,length(trimmedBallData.position)})*0.001 - (trimmedBallData.time{1,length(trimmedBallData.position)-2})*0.001; %seconds
    
    %Calculate velocity at each timesetp for checking later on:
    for f = 3:length(trimmedBallData.position)
        VX_CONT(f) = (trimmedBallData.position{f}(1,1) - trimmedBallData.position{f-2}(1,1))/dT;
        VY_CONT(f)=(trimmedBallData.position{f}(2,1) - trimmedBallData.position{f-2}(2,1))/dT;
        VZ_CONT(f)=(trimmedBallData.position{f}(3,1) - trimmedBallData.position{f-2}(3,1))/dT;
    end
    
    ballRelease = length(nanIndex) - maxElbowFrame;
    
    %% STEP 2: Physics Model
    
    %Normalize the ball position with respect to the ball travel distance
    
     % 41.5 ft =  distance from half court to the front of the rim. Interested in ball's last known position
    ballTravelDistance_x = 41.5 - trimmedBallData.position{1,length(trimmedBallData.position)}(1,1);
    
    % 13.5 ft = distance from the FT line to the front of the rim
    ballIniitialPoisiton_relative = 13.5 - ballTravelDistance_x;
    
    %find initial velocity using the 3-point central difference method to calculate both on the last
    %known ball position
    v2_x = trimmedBallData.position{1,length(trimmedBallData.position)}(1,1); %ft
    v1_x = trimmedBallData.position{1,((length(trimmedBallData.position))-2)}(1,1); %ft
    
    v2_y = trimmedBallData.position{1,length(trimmedBallData.position)}(2,1); %ft
    v1_y = trimmedBallData.position{1,((length(trimmedBallData.position))-2)}(2,1); %ft
     
    v2_z = trimmedBallData.position{1,length(trimmedBallData.position)}(3,1); %ft
    v1_z = trimmedBallData.position{1,((length(trimmedBallData.position))-2)}(3,1); %ft
    
    t2 = (trimmedBallData.time{1,length(trimmedBallData.position)})*0.001; %seconds
    t1 = (trimmedBallData.time{1,length(trimmedBallData.position)-2})*0.001; %seconds
    
    veloInitial_x = (v2_x-v1_x)/(t2-t1); %ft/s
    veloInitial_y = (v2_y-v1_y)/(t2-t1); %ft/s
    veloInitial_z = (v2_z-v1_z)/(t2-t1); %ft/s
    
    %Find the time of flight using kinematic equations, need acceleration first
    % Only acceleration  in the x direction acting on the ball is due to air resistance (drag)
    % Compute acceleration
    %A Study by McCluskey et al. (2008): This study analyzed the effect of spin on the trajectory of a basketball and reported Magnus forces in the range of 2 to 4 N for typical game situations.
    % 4 N =  0.897 lbf and the drage is approximately 0.3 lbf (Straight Shooter)
    % Tflight Should be between 0.7-1sec according to Straight Shooter By Bob J. Fisher With a Special Chapter: The Physics of the Free Throw by Larry Silverberg
    %ax = 0.897-0.3/ m; % From literature Magnus - Drag
    
    ballMass = 1.38; % Mass of the basketball (lb)
    ballVelocity = veloInitial_x;
    
    dragForce = 0.3; %lbf from literature on FT
    Magnus = 0.897; %lbf
    
    acceleration = ((-dragForce+Magnus)/ballMass)*3.280839895013123; %ft/s^2, 0.897 contribution of Magnus force
    
    %calculate flight time in x-direction by reaarranging kinematic equations and using quadratic formula
    % Should be between 0.7-1sec according to Straight Shooter By Bob J. Fisher With a Special Chapter: The Physics of the Free Throw by Larry Silverberg
    a = 0.5*acceleration;
    b = veloInitial_x;
    c = -1*ballTravelDistance_x;
    
    timeNeg = (-1*b-sqrt((b^2)-4*a*c))/(2*a);
    timePos = (-1*b+sqrt((b^2)-4*a*c))/(2*a);
    
    if timeNeg > 1.5 || timePos > 0
        flightTime = timePos;
    elseif timeNeg < 0
        flightTime = timePos;
    elseif time
    end
    
    %Pick the correct time root and assign it to flightTime
    if timeNeg > timePos || timePos > 0
        flightTime = timePos;
    elseif timePos > timeNeg || timeNeg > 0
        flightTime = timeNeg;
    elseif timeNeg > timePos || timePos < 0
        flightTime = timeNeg;
    elseif timePos > timeNeg || timeNeg < 0
        flightTime = timePos;
    elseif timeNeg < 0
        flightTime = timePos;
    elseif timePos < 0
        flightTime = timeNeg;
    end
    
    %% STEP 3: Trajectory Model
    
    % Timing
    timeSteps = linspace(0,flightTime,100);
    dt = timeSteps(1,3)-timeSteps(1,2);
    
    % Drag force
    g = 32.15; % Acceleration due to gravity (ft/s^2)
    m = 1.38; % Mass of the basketball (lb)
    
    % Initial conditions
    x0 = ballIniitialPoisiton_relative;
    y0 = trimmedBallData.position{1,length(trimmedBallData.position)}(2,1); % Initial y-position (ft)
    z0 = trimmedBallData.position{1,length(trimmedBallData.position)}(3,1); % Initial z-position (ft)
    vx0 = veloInitial_x; % Initial velocity in x-direction (ft/s)
    vy0 = veloInitial_y; % Initial velocity in y-direction (ft/s)
    vz0 = veloInitial_z; % Initial velocity in z-direction (ft/s)
    
    %Final conditions: Assume ball is at rest when it crosses hoop's frame
    x_f = 13.5+((trimmedBallData.landingY{1,1})/12); % Final x position (ft)
    y_f = ((trimmedBallData.landingX{1,1})/12); % Final y position (ft)
    z_f = 10; % Final z position (ft)
    v_xf = 0; % Final velocity in x (ft/s)
    v_yf = 0; % Final velocity in y (ft/s)
    v_zf = 0; % Final velocity in z (ft/s)
    
    if vx0 < 5
        vx0 = 15;
    end
    
    if vz0 <5
        vz0 = 15;
    end
    
    % Time span for the simulation
    t_span = [0 flightTime];
    
    initial_pos = [x0,y0,z0];
    final_pos = [x_f,y_f,z_f];
    initial_vel = [vx0,vy0,vz0];
    final_vel = [v_xf,v_yf,v_zf];
    
    % Optimization options
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
                   'StepTolerance', 1e-6, 'OptimalityTolerance', 1e-6);
    
    % Objective function: minimize the squared distance between final position
    objective = @(initial_vel_guess) objective_function(initial_vel_guess, ...
                                                initial_pos, final_pos, ...
                                                t_span, m, g);
    
    % Set bounds on the initial velocities if needed (e.g., physical limitations)
    lb = [0, -15, 0]; % Lower bound (m/s)
    ub = [50, 50, 50];    % Upper bound (m/s)
    
    % Initial guess for initial velocities (can start with provided initial velocities)
    initial_guess = initial_vel;
    
    % Call fmincon with the velocity coupling constraint
    [optimized_vel, ~] = fmincon(objective, initial_guess, [], [], [], [], lb, ub, ...
        @(v) velocity_coupling_constraint(v), options);
    
    % Display the optimized initial velocities
    disp('Optimized initial velocities (vx0, vy0, vz0):');
    disp(optimized_vel);
    
    % Simulate the optimized trajectory
    [t, sol] = simulate_trajectory(optimized_vel, initial_pos, t_span, m, g);

    
    %% Step 7: Append previously known data
    additionalTime = trimmedBallData.time{1,1}*0.001; %Seconds
    
    additionalBallPosition = trimmedBallData.position(1,1:end);
    
    for ff = 1:length(additionalBallPosition)
        X_CONT(ff) = additionalBallPosition{ff}(1,1);
        Y_CONT(ff)=additionalBallPosition{ff}(2,1);
        Z_CONT(ff)=additionalBallPosition{ff}(3,1);
    end
    
    % Only tracking the values during the shot phase (where ball height is
    % greater than participant height (1.91 m = 6.27 ft)
    
    shotPhaseIndex = max(find(Z_CONT < 6.27));
    
    appendedTrajX = horzcat((13.5-(41.5-X_CONT(1,shotPhaseIndex:end))), sol(:,1)');
    appendedTrajY = horzcat(Y_CONT(1,shotPhaseIndex:end), sol(:,2)');
    appendedTrajZ = horzcat(Z_CONT(1,shotPhaseIndex:end), sol(:,3)');
    
    %% STEP 6: Visualize
    
    x_traj = appendedTrajX;
    y_traj = appendedTrajY;
    z_traj = appendedTrajZ;
    
    % x_traj = sol(:,1);
    % y_traj = sol(:,2);
    % z_traj = sol(:,3);
    
    trajTime = linspace(0,additionalTime+t(end),length(x_traj));
    
    %Hoop
    radius = 0.75;
    center = [14.25,0,10];
    theta = linspace(0, 2*pi, 100);
    xc = radius * cos(theta) + center(1);
    yc = radius * sin(theta) + center(2);
    zc = center(3) * ones(size(theta)); % z-coordinate is constant
    plot3(xc, yc, zc, '-', 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980]);
    hold on
    
    % Trajectory
    p = plot3(x_traj', y_traj', z_traj');
    p.MarkerSize = 10;
    p.LineStyle = ":";
    p.Color = 'blue';
    
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
    
    %Animate
    % Initialize the plot for the animated point
    h = plot3(x_traj(1), y_traj(1), z_traj(1), 'ro', 'MarkerFaceColor',  [0.8500 0.3250 0.0980], 'MarkerSize', 30); 
    
    result = trimmedBallData.result{1,1};
    exitVelocity = norm(initial_vel);
    launchAngle = atan2(vz0, vx0)*(180/pi);
    entryAngle = trimmedBallData.entryAngle{1,1};
    txt = [];
    txt1 = [];
    txt2 = [];
    
    % Animate the point along the trajectory
    for n = 1:length(x_traj)
        view(-50,10);
        set(h, 'XData', x_traj(n), 'YData', y_traj(n), 'ZData', z_traj(n));  % Update position
        pause(0.005);  % Pause to create animation effect
    
    % Clear the previous text if it exists
        if ~isempty(txt) || ~isempty(txt1) || ~isempty(txt2)
            delete(txt);
            delete(txt1);
            delete(txt2);
        end
    
        % Print new text at the current point
        str = sprintf('Time (s): %d', trajTime(n)); % Example string
        txt1 = text(15, -2, 4, str,'FontSize', 20, 'Color', 'black');
    end
    
    if strcmp(result, 'missed')
        txt = text(15, -4, 5, result,'FontSize', 20, 'Color', 'r');
        txt2 = text(x_f, y_f, z_f, 'X', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'red');
        str3 = sprintf('Exit Velocity [x,y,z] (ft/s): %.1f [%.1f, %.1f, %.1f]',exitVelocity, vx0,vy0,vz0);
        str4 = sprintf('Launch/Entry Angles (deg): %.1f / %.1f',launchAngle,entryAngle);
        txt3 = text(15, -4, 4, str3,'FontSize', 20, 'Color', 'black');
        txt4 = text(15, -4, 3, str4,'FontSize', 20, 'Color', 'black');
    else
        txt = text(15, -4,5, result,'FontSize', 20, 'Color', 'green');
        txt2 = text(x_f, y_f, z_f, 'X', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'red');
        str3 = sprintf('Exit Velocity [x,y,z] (ft/s): %.1f [%.1f, %.1f, %.1f]',exitVelocity, vx0,vy0,vz0);
        str4 = sprintf('Launch/Entry Angles (deg): %.1f / %.1f',launchAngle,entryAngle);
        txt3 = text(15, -4, 4, str3,'FontSize', 20, 'Color', 'black');
        txt4 = text(15, -4, 3, str4,'FontSize', 20, 'Color', 'black');
    end
    pause(1);  % Pause to create animation effect
    
    outputDataFileName = input;

    % Define the character you want to keep around
    char_to_keep = 'BB_FT';
    
    % Find the position of the character
    char_index = strfind(outputDataFileName, char_to_keep);
    
    % If the character is found
    if ~isempty(char_index)
        % Keep everything from the character (inclusive)
        new_str = outputDataFileName(char_index:end-5);  % From the character to the end
    end

    %Save Key Variables
    save('C:\Users\Cedric Attias\OneDrive - University of Waterloo\MLSE Data Challenge\Trajectory Data\'+ convertCharsToStrings(new_str),'initial_pos',"final_pos","initial_vel",'x_traj',"y_traj","z_traj","trajTime","result");

    figures = findobj('Type', 'figure');

    % Maximize the figure window (different methods depending on the platform)
    % For Windows:
    set(figures, 'WindowState', 'maximized'); 
    % For older versions or different OS, might need other methods like:
    % set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    % Pause briefly to ensure the window maximizes (if necessary)
    pause(5);
    
    figTitle = convertCharsToStrings(new_str);
    
    % Replace invalid filename characters in title
    validTitle = regexprep(figTitle, '[^a-zA-Z0-9]', '');

    savePlots ='C:\Users\Cedric Attias\OneDrive - University of Waterloo\MLSE Data Challenge\Trajectory Graph';
    
    % Specify the full path where the figure will be saved
    filePath = fullfile(savePlots, [validTitle, '.png']);

    % Save the figure
    saveas(figures, filePath);

    close all;
end

%% FUNCTIONS FOR TRAJECTORY MODEL

function [t, sol] = simulate_trajectory(initial_vel, initial_pos, t_span, m, g)
    % Solve the trajectory using ODE45
    options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    initial_conditions = [initial_pos(1), initial_pos(2), initial_pos(3), ...
                          initial_vel(1), initial_vel(2), initial_vel(3)];
    [t, sol] = ode45(@(t, y) trajectory_dynamics(t, y, m, g), t_span, initial_conditions, options_ode);
end

function dydt = trajectory_dynamics(~,y,m, g)
% Unpack state variables
    x = y(1);
    y_pos = y(2);
    z = y(3);
    vx = y(4);
    vy = y(5);
    vz = y(6);
    
    % Compute acceleration
    %A Study by McCluskey et al. (2008): This study analyzed the effect of spin on the trajectory of a basketball and reported Magnus forces in the range of 2 to 4 N for typical game situations.
    % 4 N =  0.897 lbf and the drage is approximately 0.3 lbf (Straight Shooter)
    % Tflight Should be between 0.7-1sec according to Straight Shooter By Bob J. Fisher With a Special Chapter: The Physics of the Free Throw by Larry Silverberg
    ax = 0.897-0.3/ m; % From literature Magnus - Drag
    ay =  0;
    az =  -g;

    % Return derivatives
    dydt = [vx; vy; vz; ax; ay; az];
end

function error = objective_function(initial_vel_guess, initial_pos, final_pos, t_span,m, g,final_vel)
    % Simulate trajectory with the guessed initial velocities
    [~, sol] = simulate_trajectory(initial_vel_guess, initial_pos, t_span,m,g);
    
    % Calculate final position and velocity from the trajectory solution
    final_sim_pos = sol(end, 1:3);
    final_sim_vel = sol(end, 4:6);
    
    % Position error (distance between simulated and target final position)
    position_error = norm(final_sim_pos - final_pos);
    
    % Final z-coordinate error (emphasize z-position tracking)
    z_position_error = abs(final_sim_pos(3) - final_pos(3));
    
    % Penalize for z going below zero (ground level)
    below_ground_penalty = sum(sol(:, 3) < 0);  % Count the number of points below zero in z

    %beyond_backboard_penalty
    beyond_backboard_penalty = sum(sol(:, 1) > 15);  % Count the number of points below zero in z
    
    % Penalize high initial z-velocity
    z_velocity_penalty = initial_vel_guess(3)^2;  % Penalize high initial vertical velocity
    
    % Penalize extreme apex heights (stronger penalty if initial z-velocity is high)
    max_z = max(sol(:, 3));
    apex_penalty = max(0, max_z - final_pos(3));  % Penalize if apex goes too high
    
    % Objective: minimize position error, with penalties on z-position tracking, initial z-velocity, and apex height
    error =100000*position_error^2 + 2* beyond_backboard_penalty^2 + 1 * z_position_error^2 + 2 * below_ground_penalty^2 + 10 * apex_penalty^2 + 50 * z_velocity_penalty;  % Adjust weights as needed
end

function [c, ceq] = velocity_coupling_constraint(initial_vel_guess)
    % Extract the initial x- and z-velocities
    x_velocity = initial_vel_guess(1);
    z_velocity = initial_vel_guess(3);
    
    % Ensure that the ratio z/x is below the specified threshold
    velocity_ratio = abs(x_velocity) / abs(z_velocity);
    
    % Inequality constraint: z/x must be less than max_ratio
    c = velocity_ratio - 1.5;
    
    % No equality constraints
    ceq = [];
end
