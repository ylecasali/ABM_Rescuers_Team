function target = selectTarget( nodes, excluded_nodes )
    
    available_nodes = setdiff(nodes, excluded_nodes);
    target = available_nodes(randi(length(available_nodes)));
end

