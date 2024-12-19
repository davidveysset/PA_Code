function concentrations = create_concentrations(unique_concentrations, num_iter);
    concentrations = 0;
    for i = 1:size(unique_concentrations,1)
        row = unique_concentrations(i,:);
        tiled_row = repmat(row, num_iter, 1);
        concentrations = [concentrations; tiled_row];
    end
end
