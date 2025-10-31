% Defeine the memory tracker 
classdef MemTracker < handle
    properties
        val
    end
    methods
        function mem_track = MemTracker(v)
            mem_track.val = v;
        end

        function add(mem_track, newV)
            mem_track.val = mem_track.val + newV;
        end

        function result = get(mem_track)
            result = mem_track.val;
        end
    end
end