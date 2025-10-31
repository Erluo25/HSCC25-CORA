classdef DensityInfo < handle
    properties
        init_density
        max_density
    end

    methods
        function density_info = DensityInfo(d)
            density_info.init_density = d;
            density_info.max_density = d;
        end
        
        function update_init(density_info, newD)
            density_info.init_density = newD;
        end
        
        function update_max(density_info, newD)
            density_info.max_density = newD;
        end
        
        function result = get_init(density_info)
            result = density_info.init_density;
        end

        function result = get_max(density_info)
            result = density_info.max_density;
        end
    end

end
