classdef Base < handle
    % This is a common base class from which all the other calculation
    % classes inherit from (Gaussian, Indo, Ampac, etc)
    % It defines the common API by which Controller can use them all
    
    properties
        dataPath
        template
        params
        filename
    end
    
    methods
        function obj = Base(dataPath, template, params)
           obj.dataPath = dataPath;
           obj.template = template;
           obj.params = params;

           obj.filename = obj.template;
           f = fieldnames(obj.params);
           for i=1:length(f)
               x = f{i};
               obj.params.(x){1} = num2str(obj.params.(x){1}, '%.1f');
               if obj.params.(x){2}
                   obj.filename = [obj.filename, '_', obj.params.(x){1}];
               end
           end 
        end
        function run(obj)
        end
        function parse(obj)
        end
    end
    
end

