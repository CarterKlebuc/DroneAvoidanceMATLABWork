clc, clear, close all

num_iterations = 50;


number_of_paths = randi([1, 10]);
random_paths(number_of_paths) = RandomPath;

for i = 1:number_of_paths
    random_paths(i) = RandomPath();
end
    
figure();
hold on

for i = 1:number_of_paths
    accessed_path = random_paths(i);
    accessed_path.generate_random_values
    plot(accessed_path.x_values, accessed_path.y_values)

end

point_array = gobjects([1, number_of_paths]);

for i = 1:number_of_paths
    point = plot(random_paths(i).x_points, random_paths(i).y_points, '-o', MarkerFaceColor='red');
    point_array(i) = point;
end

%% Animate Plane moving tangent to the line
% Initialize Airplane Starting Location
for i = 1:num_iterations
    for j = 1:number_of_paths
        point_array(j).MarkerIndices = i;
    end
    exportgraphics(gca, "TestAnim.gif",Append=true)
end

% for i = 1:num_iterations
%     point_array(1).MarkerIndices = i;
%     exportgraphics(gca,"TestAnim.gif",Append=true)
% end
