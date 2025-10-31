function [result, mem, density_info] = improved_benchmark(pZ, hs, splits)
    % Since we will do maximization, so need to set the dir with -1
    dir = -hs.c;
    d = hs.d;
    
    % We are always doing the maximization in the benchmark implementation.
    mem = 0;
    pZsplit = {dir' * pZ};
    
    % Record the number of generators, name it density
    [~, temp_density] = size(pZ.E);
    density_info.init_density = temp_density;
    density_info.max_density = temp_density;
   
    for i = 0:(splits+1)
        % preinit to avoid copying
        qZnew = cell(2*length(pZsplit), 1);
        c = 0; % counter
        
        if isempty(pZsplit)
            % Empty queue, finished 
            result = 0;
            return
        else
            if i == (splits+1)
                result = 1;
                return
            end
        end
          
        for j = 1:length(pZsplit) 
            if i == 0
                % check input without splitting (efficient for zonotopic input)
                res = pZsplit(j);
            else
                res = splitLongestGen(pZsplit{j});
            end
            
            for k = 1:length(res)
                res_k = res{k};
                
                % Update the memory usage information and the density
                % information
                mem = mem + (size(res_k.G,1) * size(res_k.G,2)) + (size(res_k.E,1) * size(res_k.E,2));
                [~, temp_density] = size(res_k.E);
                
                % Update the density information and have the print outs
                if temp_density > density_info.max_density
                    density_info.max_density = temp_density;
                    %fprintf("New max density is %d\n", density_info.max_density);
                end
                
                % compute support function for enclosing zonotope
                [max_k, ~, ~] = supportFunc(zonotope(res_k),1);

                if (-1*max_k) <= d
                    % add new set to queue which need splitting
                    c = c + 1;
                    qZnew{c} = res_k; 
                end
            end
        end
        % update remaining splitted sets
        pZsplit = qZnew(1:c);
    end
end