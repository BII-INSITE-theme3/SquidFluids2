function [pre_found,index, stks,F,dUflowy,Uflowx,Uflowy] = check_existing_runs()

    % Get the current parameter file contents
    path0 = 'parameters.m';
    content0 = fileread(path0); 
    
    % Booliean to check if a precalculated version has been found
    pre_found = false;
    
    % Define directory path and pattern (e.g., all .m files)
    pattern = 'outputs/parameters_*.m'; % Define the parameter save format.
    
    % get the number of prior runs for which the data has been pre-run.
    import matlab.buildtool.io.FileCollection
    fc1 = FileCollection.fromPaths(pattern);
    
    for i = 1:length(fc1.paths)
    
        % Set file name to get contents of
        path1 = ['outputs/parameters_' num2str(i) '.m'];
        content1 = fileread(path1);
    
        % Compare content (ignoring whitespaces)
        if strcmp(strtrim(content0), strtrim(content1))
            pre_found = true;
            index = (fc1.length) + 1;
            load(['outputs/output_' num2str(i) '.mat'],"stks","F","dUflowy","Uflowx","Uflowy");
        end
    
    end

    if ~pre_found

        index = 0;

        stks = [];
        F = [];
        dUflowy = [];
        Uflowx = [];
        Uflowy = [];

    end

end

