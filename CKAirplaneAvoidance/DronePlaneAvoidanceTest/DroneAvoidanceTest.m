clc, clear, close all

%% Initialize Path
starting_location = [0, 5];
ending_location = [30, 5];
num_iterations = 100;
x_points = linspace(starting_location(1), ending_location(1), num_iterations);
y_points = linspace(starting_location(2), ending_location(2), num_iterations);

%% Initialize Drone
drone_location = [28, 5];
drone_max_speed = 2;


%% Initialize Plane (Uses Mathematica Ported Program)
plane_location = [0, 5];
plane_speed = 5;
avoidance_radius = 2;
squeeze_zone = 1;
angle = deg2rad(0);
translation = [0; 5];

%% Applying Safety Constraints to user variables
squeeze_zone = max(0.01, squeeze_zone);
if (plane_speed <= drone_max_speed)
    drone_max_speed = 0.999 * plane_speed;
end

%% Calculating Neccessary Variables
alpha = atan2(drone_max_speed, plane_speed);
beta = (pi / 2) - alpha;
w = avoidance_radius / sin(alpha);

%% Calculate Distance Between Drone and Plane
distance = sqrt(power((drone_location(1) - plane_location(1)), 2) + power((drone_location(2) - plane_location(2)), 2));
% Given that the object will be moving towards the left of the drone
object_center_point = plane_location;
object_rightmost_point = [plane_location(1) + w, plane_location(2)];
object_leftmost_point = [plane_location(1) - w, plane_location(2)];
% Perpendicular Avoidance must be a function of the three above variables
% with the height of the drone to be returned

% Function Outline:
% Inputs: Drone Current Position, Plane Center Position, Plane Leftmost
% Position, Plane Rightmost Position
% Steps: 
% 1. Divide the plot into these regions: Plane already passed, Plane
% Descending Area, Plane Midpoint (tip of triangle), Plane Ascending Area, 
% and Plane still coming (distance greater than w)
% 2. Determine which area the drone is in based on the given inputs
% If:
%   - (Plane is in ascending area)
%       a. Calculate ascending slope
%       b. Find the drone's position relative to the slope
%           - Represent this as a percentage between rightmost point and
%           center
%       c. Set drone height to be (slope * percentage) + offset



%% Plot Path and Drone
plot(x_points, y_points);
hold on
scatter(drone_location(1), drone_location(2));
scatter(object_rightmost_point(1), object_rightmost_point(2))
scatter(object_leftmost_point(1), object_leftmost_point(2))


%% Creating Evacuation Zone Polygon
% Need to apply a rotation matrix to each set of verticies in the polygon
% to properly rotate each evacuation zone
% Test Angle
rotation_matrix = [cos(angle), -sin(angle); sin(angle), cos(angle)];

%   Evacuation Zone 1
evac_zone_1_x = [0, avoidance_radius * cos(beta), w, avoidance_radius * cos(beta)];
evac_zone_1_y = [0, avoidance_radius * sin(beta), 0, avoidance_radius * -sin(beta)];

%   Evacuation Zone 2
evac_zone_2_x = [0, -avoidance_radius * cos(beta), -w, -avoidance_radius * cos(beta)];
evac_zone_2_y = [0, -avoidance_radius * sin(beta), 0, -avoidance_radius * -sin(beta)];

evac_zone_1_x_rot = zeros([1 4]);
evac_zone_1_y_rot = zeros([1 4]);
evac_zone_2_x_rot = zeros([1 4]);
evac_zone_2_y_rot = zeros([1 4]);

for i = 1:4
    original_verticies = [evac_zone_1_x(i); evac_zone_1_y(i)];
    rotated_verticies = (rotation_matrix * original_verticies) + translation;
    evac_zone_1_x_rot(i) = rotated_verticies(1);
    evac_zone_1_y_rot(i) = rotated_verticies(2);
end

for i = 1:4
    original_verticies = [evac_zone_2_x(i); evac_zone_2_y(i)];
    rotated_verticies = (rotation_matrix * original_verticies) + translation;
    evac_zone_2_x_rot(i) = rotated_verticies(1);
    evac_zone_2_y_rot(i) = rotated_verticies(2);
end

%% Creating Squeeze Zone Rectangle
% Redo Rectangle so it uses polyshape instead of rectangle
rectangle_x_coordinates = [-w, -w, w, w];
rectangle_y_coordinates = [avoidance_radius + squeeze_zone, -avoidance_radius - squeeze_zone, -avoidance_radius - squeeze_zone, avoidance_radius + squeeze_zone];

rectangle_starting_location = [-w, -avoidance_radius - squeeze_zone];
rectangle_width = 2 * w;
rectangle_height = 2 * (avoidance_radius + squeeze_zone);

rectangle_x_coordinates_rot = zeros([1 4]);
rectangle_y_coordinates_rot = zeros([1 4]);

for i = 1:4
    original_verticies = [rectangle_x_coordinates(i); rectangle_y_coordinates(i)];
    rotated_verticies = (rotation_matrix * original_verticies) + translation;
    rectangle_x_coordinates_rot(i) = rotated_verticies(1);
    rectangle_y_coordinates_rot(i) = rotated_verticies(2);
end


my_poly = polyshape(rectangle_x_coordinates_rot, rectangle_y_coordinates_rot);

%% Plotting Polygons
airplane_location = [0, 0];
scatter(translation(1, 1), translation(2, 1), 'filled', c="blue")
%rectangle('Position', [rectangle_starting_location(1), rectangle_starting_location(2), rectangle_width, rectangle_height], FaceColor='y'); 
plot(my_poly)


fill(evac_zone_1_x_rot, evac_zone_1_y_rot,'g');
fill(evac_zone_2_x_rot, evac_zone_2_y_rot,'g');

%% Plot Avoidance Radius!
% Viscircles can be used as a debug method for the rectangle-based circle!
%viscircles(airplane_location, avoidance_radius);
pos = [-avoidance_radius + translation(1, 1), -avoidance_radius + translation(2, 1), avoidance_radius * 2, avoidance_radius * 2];
rectangle('Position',pos,'Curvature',[1 1], FaceColor='r', EdgeColor='r')

%% Creating Squeeze Zone Trendlines
wid = rectangle_width;    % Half-width of outer boundary
L = 0.1;

imax = 8;
for i = 1:imax
    row = (i - 1) / (imax - 1);

    % Start with left horizontal line
    x_left = [-wid, -w];
    y_left = [row * (avoidance_radius + squeeze_zone), row * (avoidance_radius + squeeze_zone)];

    % Arc part (Beta_v from pi - beta to beta in -beta/12 steps)
    beta_v = linspace(pi - beta, beta, 24);  % Approximate step count

    x_arc = avoidance_radius * cos(beta_v);
    y_arc = avoidance_radius * sin(beta_v) + row * (squeeze_zone + avoidance_radius - avoidance_radius * sin(beta_v));

    % Add incline line
    x_inc = [x_left(end), x_arc(1)];
    y_inc = [y_left(1), y_arc(1)];

    % Right horizontal line
    x_right = [w, wid];
    y_right = [row * (avoidance_radius + squeeze_zone), row * (avoidance_radius + squeeze_zone)];

    % Add decline line
    x_dec = [x_arc(end), x_right(1)];
    y_dec = [y_arc(end), y_right(1)];

    % Combine full line path
    x_full = [x_inc, x_arc, x_dec];
    y_full = [y_inc, y_arc, y_dec];

    x_full_rot = zeros([1 28]);
    y_full_rot = zeros([1 28]);

    for j = 1:28
        original_verticies = [x_full(j); y_full(j)];
        rotated_verticies = (rotation_matrix * original_verticies) + translation;
        x_full_rot(j) = rotated_verticies(1);
        y_full_rot(j) = rotated_verticies(2);
    end
    % Draw the line
    line(x_full_rot, y_full_rot);  % 'k' = black line

end

for i = 1:imax
    row = -(i - 1) / (imax - 1);

    % Start with left horizontal line
    x_left = [-wid, -w];
    y_left = [row * (avoidance_radius + squeeze_zone), row * (avoidance_radius + squeeze_zone)];


    % Arc part (Beta_v from pi - beta to beta in -beta/12 steps)
    beta_v = linspace(pi - beta, beta, 24);  % Approximate step count

    x_arc = avoidance_radius * cos(beta_v);
    y_arc = -avoidance_radius * sin(beta_v) + row * (squeeze_zone + avoidance_radius - avoidance_radius * sin(beta_v));

    % Add incline line
    x_inc = [x_left(end), x_arc(1)];
    y_inc = [y_left(1), y_arc(1)];

    % Right horizontal line
    x_right = [w, wid];
    y_right = [row * (avoidance_radius + squeeze_zone), row * (avoidance_radius + squeeze_zone)];

    % Add decline line
    x_dec = [x_arc(end), x_right(1)];
    y_dec = [y_arc(end), y_right(1)];

    % Combine full line path
    x_full = [x_inc, x_arc, x_dec];
    y_full = [y_inc, y_arc, y_dec];

    x_full_rot = zeros([1 28]);
    y_full_rot = zeros([1 28]);

    for j = 1:28
        original_verticies = [x_full(j); y_full(j)];
        rotated_verticies = (rotation_matrix * original_verticies) + translation;
        x_full_rot(j) = rotated_verticies(1);
        y_full_rot(j) = rotated_verticies(2);
    end

    % Draw the line
    line(x_full_rot, y_full_rot);  % 'k' = black line
end

xlim([-5, 35])
ylim([0, 20])
xlabel("X Position [m]")
ylabel("Y Position [m]")
axis equal;