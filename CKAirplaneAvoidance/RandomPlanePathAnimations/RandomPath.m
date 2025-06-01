classdef RandomPath
    properties
        x_values
        y_values
        x_points
        y_points
        num_iterations;
        starting_location
        ending_location
    end
    
    methods
        function obj = RandomPath()
            obj.x_values = [0, 0];
            obj.y_values = [0, 0];
            obj.x_points = [];
            obj.y_points = [];
            obj.num_iterations = 50;
            obj.starting_location = [0, 0];
            obj.ending_location = [0, 0];
            
            obj = obj.generate_random_values();
        end
        
        function obj = generate_random_values(obj)
            obj.x_values = randi([-10, 10], 1, 2);
            obj.y_values = randi([-10, 10], 1, 2);
            obj.starting_location = [obj.x_values(1), obj.y_values(1)];
            obj.ending_location = [obj.x_values(2), obj.y_values(2)];
            obj.x_points = linspace(obj.starting_location(1), obj.ending_location(1), obj.num_iterations);
            obj.y_points = linspace(obj.starting_location(2), obj.ending_location(2), obj.num_iterations);
        end
    end
end

