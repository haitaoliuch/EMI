function [flatten_time,transmission] = Coronavirus(varargin)
% DESCRIPTION
% This function simulates disease transmission amongst a set of n-carriers,
% in a confined space, with some proportion of carriers socially isolating
% themselves. A simple multibody physics model (elastic collision between
% two equal-mass particles) determines carrier trajectory.
%
% OPTIONAL ARGUMENTS
%       'lim': spatial limits [scalar > 0]
%       'n': number of carriers [scalar > 0]
%       'w_loc': vertical quarantine wall x-location (as proportion of lim, set to 1 if no quarantine wall is desired) [scalar > 0]
%       'w_loc_h': horizontal wall y-location (as proportion of lim, set to 1 if no quarantine wall is desired)
%       'w_init': date when quarantine wall opens [in days] [scalar > 0]
%       'w_speed': gate opening speed [scalar > 0]
%       'w_open' : gate opening width (as a percent of the spatial dimension) [0 < scalar < 1]
%       'rad': infection radius in units [scalar > 0]
%       'speed': initial carrier speed in units/day [scalar > 0]
%       'iso': percentage of social isolation [0 <= scalar < 1]
%       't_tot': total simulation time in days [scalar > 0]
%       'delT': simulation time increment in days [scalar > 0]
%       't_recovery': mean recovery time in days [scalar > 0]
%       'p_init': probability of carriers initially infected [0 <= scalar < 1]
%       'p_trans': probability of disease transmission [0 <= scalar < 1]
%       'p_mort': mortality rate [0 <= scalar < 1]
%       'video_save': boolean to save video [true/false]
%       'verbose': boolean to print simulation results to command [true/false]
%       'filename': filename to save video [string]
%
% OUTPUT
%       transmission = a cell array containing the time vector, percentage
%       of unaffected carriers, percentage of infected carriers, percentage
%       of recovered carriers, and percentage of deceased carriers
%
% EXAMPLE USAGE
%       transmission = simulitis('lim',200,'n',200,'rad',1,'speed',10,...
%           'iso',.5,'t_tot',40,'delT',0.1,'t_recovery',14,'p_init',.05,...
%           'p_trans',.95,'p_mort',.05,'video_save',true,'verbose',true,...
%           'filename','video_1')
%

% If running outside the function, parse inputs
if nargin
    p = inputParser;
    
    % Inline functions to check parameter validity
    validParam = @(x) isnumeric(x) && isscalar(x) && (x>=0);
    validParamWhole = @(x) validParam(x) && ~rem(x,1);
    validProportion = @(x) validParam(x) && x<=1;
     
    addOptional(p, 'lim', 200, validParamWhole);
    addOptional(p, 'n', 200, validParamWhole);
    addOptional(p, 'w_loc', .4, validProportion);
    addOptional(p, 'w_loc_h', .4, validProportion);
    addOptional(p, 'w_init', 5, validParam);
    addOptional(p, 'w_speed', 10, validParam);
    addOptional(p, 'w_open', .5, validProportion);
    addOptional(p, 'rad', 1.5, validParam);
    addOptional(p, 'speed', 10, validParam);
    addOptional(p, 'iso', 0.5,  validProportion);
    addOptional(p, 't_tot', 14, validParam);
    addOptional(p, 'delT', 0.1, validParam);
    addOptional(p, 't_recovery', 14, validParam);
    addOptional(p, 'p_init', 0.05, validProportion);
    addOptional(p, 'p_trans', 0.95, validProportion);
    addOptional(p, 'p_mort', 0.03, validProportion);
    addOptional(p, 'video_save', true, @(x) islogical(x));
    addOptional(p, 'verbose', true, @(x) islogical(x));
    addOptional(p, 'filename', 'file', @(x) ischar(x));
    parse(p,varargin{:})
    
    % parse the inputparser object results
    lim = p.Results.lim;
    n = p.Results.n;
    w_loc = p.Results.w_loc;
    w_loc_h = p.Results.w_loc_h;
    w_init = p.Results.w_init;
    w_speed = p.Results.w_speed;
    w_open = p.Results.w_open;
    rad = p.Results.rad;
    speed = p.Results.speed;
    iso = p.Results.iso;
    t_tot = p.Results.t_tot;
    delT = p.Results.delT;
    t_recovery = p.Results.t_recovery;
    p_init = p.Results.p_init;
    p_trans = p.Results.p_trans;
    p_mort = p.Results.p_mort;
    video_save = p.Results.video_save;
    verbose = p.Results.verbose;
    filename = p.Results.filename;
    
    % If running from within the function, you can define model inputs here
else
    % --------MODEL INPUTS----------
    lim = 100;          % Spatial limits
    n = 100;            % Number of carriers
    w_loc = 0.8;         % Vertical quarantine Wall location (as a proportion of lim). If you don't want a wall, set to 1.
    w_loc_h = .2;       % Horizontal quarantine Wall location (as a proportion of lim). If you don't want a wall, set to 1.
    w_init = 2;         % Wall open start time (in days)
    w_speed = 3;        % Wall opening speed (in units/day)
    w_open = .5;        % Total wall opening width (as a fraction of spatial dimension)
    rad = 1;            % Infection radius
    speed = 10;         % Initial speed of carriers (in units/day)
    t_recovery = 14;    % Mean recovery time (in days)
    t_tot = 40;         % Total simulation time (in days)
    delT = 0.02;        % Time increment (in days)
    iso = 0.3;          % Social isolation proportion
    p_init = 0.1;       % Percent of carriers initially infected
    p_trans = 0.99;     % Disease transmission rate
    p_mort = 0.03;      % Mortality rate
    video_save = true;  % Save video?
    verbose = true;     % Print simulation results to command prompt
    filename = 'video.avi'; % Video filename
end
% --------VECTOR ALLOCATION----------
pos = rand(n,2).*(lim);                 % Initial Position vector
isolate = (rand(n,1)<iso);              % Initial Isolation vector
th = rand(n,1).*2.*pi;                  % Initial Direction vector
v = [speed.*cos(th) speed.*sin(th)];    % Velocity Vector
% Depending on choice of quarantine wall locations, either have gates open
% from where the horizontal and vertical walls meet, or if only one wall is
% used, have gates open at the midpoint of that wall.
% If no vertical wall is used, put it outside of frame, and have the
% horizontal wall open from the middle
if w_loc==1                           
    w_loc=1.5;
    wall_x2 = 0.5*lim;          % Left quarantine gate x-coords
    wall_x1 = 0.5*lim;          % Right quarantine gate x-coords
else
    wall_x2 = lim*w_loc;        % Left quarantine gate x-coords
    wall_x1 = lim*w_loc;        % Right quarantine gate x-coords
end
% If no horizontal wall is used, put it outside of frame, and have the
% vertical wall open from the middle.
if w_loc_h==1
    w_loc_h=1.5;
    wall_y1 = 0.5*lim;          % Top quarantine gate y-coords
    wall_y2 = 0.5*lim;          % Bottom quarantine gate y-coords
else
    wall_y1 = lim*w_loc_h;      % Top quarantine gate y-coords
    wall_y2 = lim*w_loc_h;      % Bottom quarantine gate y-coords
end
wall_x = lim*w_loc;         % x-location of vertical wall
wall_y = lim*w_loc_h;       % y-location of horizontal wall.
% Put all initially infected in lower-right quadrant
infected = zeros(n,1);  
if w_loc<1 && w_loc_h<1
    quar_ind = find(pos(:,1)<wall_x1 & pos(:,2)<wall_y1);   % Find positions of quarantined carriers
    n_quar = size(quar_ind,1);                              % Number of quarantined carriers
    infected(quar_ind) = rand(n_quar,1)<p_init;             % Only infect those in quarantine
    
% If no walls, infect population uniformly
else
    infected = rand(n,1)<p_init;
end
healthy = ~infected;                    % Logical array of healthy carriers
recovered = zeros(n,1);                 % Logical array of recovered carriers
dead = zeros(n,1);                      % Logical array of deaths
p_death = rand(n,1)<p_mort;             % Probability that a person will die
t_rec = ceil(t_recovery.*ones(n,1)+4.*randn(n,1))./delT; % Time-to-recover vector (randomized about mean)
t = linspace(0,t_tot,t_tot./delT);      % Time vector
collision = zeros(n,n);                 % Matrix to keep track of recent collisions
n_delay = ceil(1/(delT));                 % Collision delay
% Allocate space for solution vectors
inf_sum = zeros(t_tot/delT,0);          % Percentage of infected carriers
hea_sum = inf_sum;                      % Percentage of unaffected carriers
rec_sum = inf_sum;                      % Percentage of recovered carriers
dead_sum = inf_sum;                     % Percentage of dead carriers
cumulative_sum = inf_sum;               % Percentage of cumulative disease cases
collision_gate = zeros(n,2);            % Gate collision delay matrix
n_delay_gate = 2;                       % Gate collision delay parameter
% Initialize videowriter
if video_save
    vid = VideoWriter(filename,'Uncompressed AVI');
    open(vid);
end
for i = 1:(t_tot/delT)
    
    % Decrement collision delay
    collision = collision-ones(n,n);
    collision(collision<0)=0;
    
    % Update carrier position
    pos_new = pos+v.*(~repmat(isolate,1,2)).*delT;
    
    % Step through each carrier
    for k = 1:n
        
        % If recovery time is up, carrier is either recovered or dead
        if infected(k)&&t_rec(k)<=0
            
            % If recovery time is up and carrier is dead, well, it's dead.
            % Zero it's velocity
            if p_death(k)
                dead(k) = 1;
                v(k,:) = [0 0];
                recovered(k)=0;
                infected(k)=0;
                healthy(k)=0;
            else
                
                % If recovery time is up and carrier is not dead, recover
                % that carrier
                recovered(k)=1;
                infected(k)=0;
                healthy(k)=0;
                dead(k)=0;
            end
            
            % If carrier is infected and not recovered, decrement recovery time
        elseif (infected(k))
            t_rec(k) = t_rec(k)-1;
        end
        
        % Step through all other carriers, looking for collisions, and if
        % so, transmit disease and recalculate trajectory
        for j = 1:1:n
            
            if j~=k
                
                % Get positions of carriers j and k
                pos1 = pos_new(k,:);
                pos2 = pos_new(j,:);
                
                % If collision between two living specimens, re-calcuate
                % direction and transmit virus (but don't check the same
                % two carriers twice)
                if norm(pos1-pos2)<=(2*rad) && ~collision(k,j) && ~collision(j,k)
                    
                    % Create collision delay (i.e. if carrier j and k have
                    % recently collided, don't recompute collisions for a
                    % n_delay time steps in case they're still close in proximity,
                    % otherwise they might just keep orbiting eachother)
                    collision(k,j) = n_delay;
                    collision(j,k) = n_delay;
                    
                    % Compute New Velocities
                    phi = atan2((pos2(2)-pos1(2)),(pos2(1)-pos1(1)));
                    
                    % if one carrier is isolated, treat it like a wall and
                    % bounce the other carrier off it
                    if isolate(j)||dead(j)
                        
                        % Get normal direction vector of 'virtual wall'
                        phi_wall = -phi+pi/2;
                        n_wall = [sin(phi_wall) cos(phi_wall)];
                        dot = v(k,:)*n_wall';
                        
                        % Redirect non-isolated carrier
                        v(k,1) = v(k,1)-2*dot*n_wall(1);
                        v(k,2) = v(k,2)-2*dot*n_wall(2);
                        v(j,1) = 0;
                        v(j,2) = 0;
                        
                    elseif isolate(k)||dead(k)
                        
                        % Get normal direction vector of 'virtual wall'
                        phi_wall = -phi+pi/2;
                        n_wall = [sin(phi_wall) cos(phi_wall)];
                        dot = v(j,:)*n_wall';
                        
                        % Redirect non-isolated carrier
                        v(j,1) = v(j,1)-2*dot*n_wall(1);
                        v(j,2) = v(j,2)-2*dot*n_wall(2);
                        v(k,1) = 0;
                        v(k,2) = 0;
                        
                        % Otherwise, transfer momentum between carriers
                    else
                        
                        % Get velocity magnitudes
                        v1_mag = sqrt(v(k,1)^2+v(k,2)^2);
                        v2_mag = sqrt(v(j,1)^2+v(j,2)^2);
                        
                        % Get directions
                        th1 = atan2(v(k,2),v(k,1));
                        th2 = atan2(v(j,2),v(j,1));
                        
                        % Compute new velocities
                        v(k,1) = v2_mag*cos(th2-phi)*cos(phi)+v1_mag*sin(th1-phi)*cos(phi+pi/2);
                        v(k,2) = v2_mag*cos(th2-phi)*sin(phi)+v1_mag*sin(th1-phi)*sin(phi+pi/2);
                        v(j,1) = v1_mag*cos(th1-phi)*cos(phi)+v2_mag*sin(th2-phi)*cos(phi+pi/2);
                        v(j,2) = v1_mag*cos(th1-phi)*sin(phi)+v2_mag*sin(th2-phi)*sin(phi+pi/2);
                        
                    end
                    
                    % If either is infected and not dead...
                    if((infected(j)||infected(k)))&&((~dead(k)||~dead(j)))
                        
                        % If either is recovered, no transmission
                        if recovered(k)
                            infected(k)=0;
                        elseif recovered(j)
                            infected(j)=0;
                            
                            % Otherwise, transmit virus
                        else
                            transmission = rand(1)<p_trans;
                            if transmission
                                infected(j)=1;
                                infected(k)=1;
                                healthy(j)=0;
                                healthy(k)=0;
                            end
                        end
                    end
                end
            end
        end
        
        % Look for collisions with outer walls and re-direct
        
        % Left Wall
        if pos_new(k,1)<=rad
            if v(k,1)<0
                v(k,1)=-v(k,1);
            end
            
            % Right wall
        elseif pos_new(k,1)>=lim-rad
            if v(k,1)>0
                v(k,1)=-v(k,1);
            end
        end
        
        % Bottom Wall
        if pos_new(k,2) <=rad
            if v(k,2)<0
                v(k,2)=-v(k,2);
            end
            
            % Top Wall
        elseif pos_new(k,2) >=lim-rad
            if v(k,2)>0
                v(k,2)=-v(k,2);
            end
        end
        
        % If there is a quarantine wall, look for collisions with the wall
        % and re-direct
        if w_loc<1 || w_loc_h<1 
            
            
            % Decrement gate collision delay
            collision_gate = collision_gate-ones(n,2);
            collision_gate(collision_gate<0)=0;
    
            % Bottom Gate
            if (pos_new(k,1) > (wall_x-rad) && pos(k,1) <= (wall_x-rad))&& (pos_new(k,2)<wall_y1) && ~collision_gate(k,1)
                v(k,1)=-v(k,1);
                collision_gate(k,1) = n_delay_gate;
                
            % Top Gate
            elseif (pos_new(k,1) > (wall_x-rad) && pos(k,1) <= (wall_x-rad))&& (pos_new(k,2)>wall_y2) && ~collision_gate(k,2)
                v(k,1) = -v(k,1);
                collision_gate(k,2) = n_delay_gate;
                
            % Bottom Gate
            elseif (pos_new(k,1) < (wall_x+rad) && pos(k,1) >= (wall_x+rad))&& (pos_new(k,2)<wall_y1) && ~collision_gate(k,1)
                v(k,1)=-v(k,1);
                collision_gate(k,1) = n_delay_gate;
                
            % Top Gate
            elseif (pos_new(k,1) < (wall_x+rad) && pos(k,1) >= (wall_x+rad))&& (pos_new(k,2)>wall_y2) && ~collision_gate(k,2)
                v(k,1) = -v(k,1);
                collision_gate(k,2) = n_delay_gate;
                
            % Left Gate
            elseif (pos_new(k,2) > (wall_y-rad) && pos(k,2) <= (wall_y-rad))&& (pos_new(k,1)<wall_x1) && ~collision_gate(k,1)
                v(k,2)=-v(k,2);
                collision_gate(k,2) = n_delay_gate;
                
            % Right Gate
            elseif (pos_new(k,2) > (wall_y-rad) && pos(k,2) <= (wall_y-rad))&& (pos_new(k,1)>wall_x2) && ~collision_gate(k,2)
                v(k,2) = -v(k,2);
                collision_gate(k,2) = n_delay_gate;
                
            % Left Gate
            elseif (pos_new(k,2) < (wall_y+rad) && pos(k,2) >= (wall_y+rad))&& (pos_new(k,1)<wall_x1) && ~collision_gate(k,1)
                v(k,2)=-v(k,2);
                collision_gate(k,1) = n_delay_gate;
                
            % Right Gate
            elseif (pos_new(k,2) < (wall_y+rad) && pos(k,2) >= (wall_y+rad))&& (pos_new(k,1)>wall_x2) && ~collision_gate(k,2)
                v(k,2) = -v(k,2);
                collision_gate(k,2) = n_delay_gate;
            end
        end
    end
    
    % Open gate after alotted quarantine time
    if i*delT>w_init && (((wall_y2-wall_y1)<w_open*lim)||((wall_x2-wall_x1)<w_open*lim))
        wall_y1 = wall_y1-w_speed*delT;
        wall_y2 = wall_y2+w_speed*delT;
        wall_x1 = wall_x1-w_speed*delT;
        wall_x2 = wall_x2+w_speed*delT;
    end
    
    % Update color vector
    color = [infected healthy recovered].*(1-dead);
    
    % Update solution vectors
    inf_sum(i) = sum(infected)*100/n;
    hea_sum(i) = sum(healthy)*100/n;
    rec_sum(i) = sum(recovered)*100/n;
    dead_sum(i) = sum(dead)*100/n;
    cumulative_sum(i) = 100-hea_sum(i);
    
    if verbose
    % Initialize plots on first loop iteration
    if i==1
        % Plot transmission simulation
        figure(1);
        lineWidth = 2;
        markerSize = 16;
        set(gcf,'DefaultFigureRendererMode','auto');
        set(gcf,'Position',[100 100 800 500]);
        subplot(2,3,[1 2 4 5]);
        h = scatter(pos_new(:,1),pos_new(:,2),markerSize,color,'filled','MarkerEdgeColor','k'); hold on;
        ha = plot([wall_x wall_x],[0 wall_y1],'k-','LineWidth',lineWidth);hold on;
        hb = plot([wall_x wall_x],[wall_y2 lim],'k-','LineWidth',lineWidth);hold on;
        hc = plot([0 wall_x1],[wall_y wall_y],'k-','LineWidth',lineWidth);hold on;
        hd = plot([wall_x2 lim],[wall_y wall_y],'k-','LineWidth',lineWidth);hold off;
        xlim([0,lim]);
        ylim([0,lim]);
        axis square;
        grid on;
        box on;
        set(gca,'YTickLabel',[],'FontName', 'Times New Roman');
        set(gca,'XTickLabel',[]);
        titlestring = strcat('Percent: Unaffected= ',num2str(hea_sum(i)),', Infected= ',...
            num2str(inf_sum(i)),', Recovered=', num2str(rec_sum(i)),...
            ', Death=', num2str(dead_sum(i)));
        title(titlestring);
        
        % Resize markers to match infection radius
        currentunits = get(gca,'Units');
        set(gca, 'Units', 'Points');
        axpos = get(gca,'Position');
        set(gca, 'Units', currentunits);
        markerWidth = 2*rad/diff(xlim)*axpos(3); % Calculate Marker width in points
        lineWidth = 0.5*markerWidth;
        set(h, 'SizeData', markerWidth^2);
        markerSize = markerWidth^2;
        
        % Plot infection rates vs. time
        subplot(2,3,[3 6]);
        h2 = plot(t(1:i),hea_sum(1:i),'g','LineWidth',2);hold on;
        h3 = plot(t(1:i),inf_sum(1:i),'r','LineWidth',2);hold on;
        h4 = plot(t(1:i),rec_sum(1:i),'b','LineWidth',2);hold on; 
        h5 = plot(t(1:i),dead_sum(1:i),'k','LineWidth',2);hold off;
        legend('Unaffected','Infected','Recovered','Death');
        xlabel('Days');
        ylabel('Percent of Population');
        xlim([0,t_tot]);
        ylim([0,100]);
        set(gcf,'Color','w'); set(gca,'FontName', 'Times New Roman','FontSize',10)
        
        % Update data on subesequent iterations
    else
        subplot(2,3,[1 2 4 5]);
        set(h,'XData',pos_new(:,1));
        set(h,'YData',pos_new(:,2));
        set(h,'CData',color);
        set(ha,'XData',[wall_x wall_x]);
        set(ha,'YData',[0 wall_y1]);
        set(hb,'XData',[wall_x wall_x]);
        set(hb,'YData',[wall_y2 lim]);
        set(hc,'XData',[0 wall_x1]);
        set(hc,'YData',[wall_y wall_y]);
        set(hd,'XData',[wall_x2 lim]);
        set(hd,'YData',[wall_y wall_y]);
        
                % Update title
        titlestring = strcat('Percent: Unaffected= ',num2str(hea_sum(i)),', Infected= ',...
            num2str(inf_sum(i)),', Recovered=', num2str(rec_sum(i)),...
            ', Death=', num2str(dead_sum(i)));
        title(titlestring);
        
        subplot(2,3,[3 6]);
        set(h2,'XData',t(1:i)); set(h2,'YData',hea_sum(1:i));
        set(h3,'XData',t(1:i)); set(h3,'YData',inf_sum(1:i));
        set(h4,'XData',t(1:i)); set(h4,'YData',rec_sum(1:i));
        set(h5,'XData',t(1:i)); set(h5,'YData',dead_sum(1:i));
        
    end 
    drawnow;
    end
    % Print daily results to command prompt
    if ~mod(i,1/delT)&&verbose
        fprintf('Day %i: %.2f percent infected, %.2f percent deceased\n',i*delT,    (i),dead_sum(i));
        if i==t_tot/delT
            fprintf('Simulation Complete\n');
            [~,flatten_time] = max(inf_sum);
            flatten_time = flatten_time*delT;
            fprintf('Curve flattened after %.1f  days\n',flatten_time);
        end
    end
            [~,flatten_time] = max(inf_sum);
            flatten_time = flatten_time*delT;
    % Save video
    if video_save
        frame = getframe(gcf);
        writeVideo(vid,frame);
    end
    
    pos = pos_new;      % Update position
end
if video_save
    close(vid);
end
% Save results to a cell array
transmission = {'time',t};
transmission(2,:) = {'unaffected',hea_sum};
transmission(3,:) = {'infected',inf_sum};
transmission(4,:) = {'recovered',rec_sum};
transmission(5,:) = {'Death',dead_sum};
end