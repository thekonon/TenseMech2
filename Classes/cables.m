classdef cables

    properties
        count
        specific_stiffness
        specific_dumpings
        stiffness
        dumpings
        variable_cables_indexes
        fixed_cables_indexes
        from_to
        free_length
        free_lengths_cables %jen pomocná vlastnost - pro generátor tensegrit
        connectivity_matrix
    end
    
    methods
        function obj = cables()

        end
    end
end

