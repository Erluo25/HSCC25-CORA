function [result, mem] = improved_benchmark(pZ, hs, splits)
    % Since we will do maximization, so need to set the dir with -1
    dir = -hs.c;
    d = hs.d;
    
    % We are always doing the maximization in the benchmark implementation.
    mem = 0;
    pZsplit = {dir' * pZ};

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
                mem = mem + (size(res_k.G,1) * size(res_k.G,2)) + (size(res_k.E,1) * size(res_k.E,2));
    
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


 
%{
function [val, mem] = aux_support_benchmark(pZ, dir, splits, d)
    % We are always doing the maximization in the benchmark implementation.
    mem = 0;
    if isa(pZ, 'cell')
        pZsplit = pZ;
    else
        pZsplit = {pZ};
    end

    % project the polynomial zonotope onto the direction
    for i=1:length(pZsplit)
        pZsplit{i} = dir' * pZsplit{i};
    end

    % goal is to find a tight upper bound of all sets in pZsplit
    
    % determine lower bound of upper bound
    % is used to exclude entire sets with no need to split them further
    minUpperBound = -inf;
    
    % the result will be the smallest upper bound of the upper bound
    maxUpperBound = -inf; % result for empty set
    
    % for optimizations, we also store the largest exact upper bound
    maxExactUpperBound = -inf;
    
    % split the polynomial zonotope multiple times to obtain a better 
    % over-approximation of the real shape
    
    for i = 0:splits
        % preinit to avoid copying
        qZnew = cell(2*length(pZsplit), 1);
        c = 0; % counter
    
        % reset to only consider leftover splitted subsets
        maxUpperBound = maxExactUpperBound;
    
        for j = 1:length(pZsplit) 
            if i == 0
                % check input without splitting (efficient for zonotopic input)
                res = pZsplit(j);
            else
                res = splitLongestGen(pZsplit{j});
            end
            
            for k = 1:length(res)
                res_k = res{k};
                mem = mem + (size(res_k.G,1) * size(res_k.G,2)) + (size(res_k.E,1) * size(res_k.E,2));
    
                % compute support function for enclosing zonotope
                [max_k,~,alpha] = supportFunc(zonotope(res_k),1);
                
                % update upper and lower bound
                maxUpperBound = max(maxUpperBound, max_k);
                
                if max_k >= minUpperBound
                    % update min upper bound by 'most critical' point
                    % aka largest point in zonotope subset
    
                    % exract zonotopic generators from E
                    ind1 = sum(res_k.E,1) == 1;
                    ind2 = sum(res_k.E(:,ind1),2) == 1;
                    alpha_ = zeros(size(res_k.E,1),1);
                    alpha_(ind2) = alpha(ind1);
                    
                    % use result from zonotope supportFunc
                    minMax_k = res_k.c + ...
                        sum(res_k.G .* prod(alpha_.^res_k.E,1)); 
                    
                    if ~isempty(res_k.GI)
                        % same for GI
                        beta = alpha(size(res_k.E,2)+1:end);
                        minMax_k = minMax_k + res_k.GI*beta;
                    end
    
                    minUpperBound = max(minUpperBound, minMax_k);
    
                    if withinTol(minMax_k, max_k)
                        % found exact upper bound for current set
                        maxExactUpperBound = max(maxExactUpperBound, max_k);
                        continue;
                    end
                    
                    % TODO other min upper bound? 
                    % update min upper bound by largest random point (slow)
                    % minMax_k = max(res_k.randPoint(11));
                    % minUpperBound = max(minUpperBound, minMax_k);
    
                    % add new set to queue
                    c = c + 1;
                    qZnew{c} = res_k; 
                end
            end
        end
    
        if withinTol(minUpperBound, maxUpperBound)
            % exact upper bound is found
            val = maxUpperBound;
            return
        end
        
        % update remaining splitted sets
        pZsplit = qZnew(1:c);
    end
    
    % return result
    val = maxUpperBound;
end
%}