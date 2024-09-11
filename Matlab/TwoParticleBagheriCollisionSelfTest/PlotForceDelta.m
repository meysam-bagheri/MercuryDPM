hold all;
fstat_path ='/nishome/meysam/MercuryDPM/MercuryBuild/Drivers/SelfTests/Interactions/test5/'; %<-- Change this to the path of saved fstat files

% List all files matching the pattern
file_list = dir([fstat_path '*.fstat.*']);

% Extract the numeric part using a regular expression and convert it to a number
file_numbers = cellfun(@(x) str2double(regexp(x, '\d+$', 'match')), {file_list.name});

% Get the maximum number
last_number = max(file_numbers);

% Display the last number
disp(['Number of fstat files: ', num2str(last_number)]);

num_fstat_files = last_number;

% make array for Force and Overlap (delta)
force = zeros(num_fstat_files,1);
overlap = zeros(num_fstat_files,1);

for i = 1:num_fstat_files
    if i<10
        fstat = read_fstat([fstat_path 'TwoParticleBagheriCollisionSelfTest.fstat.000',num2str(i)]);
    elseif i<100
        fstat = read_fstat([fstat_path 'TwoParticleBagheriCollisionSelfTest.fstat.00',num2str(i)]);
    elseif i<1000
        fstat = read_fstat([fstat_path 'TwoParticleBagheriCollisionSelfTest.fstat.0',num2str(i)]);
    else
        fstat = read_fstat([fstat_path 'TwoParticleBagheriCollisionSelfTest.fstat.',num2str(i)]);
    end
    
    
    if (isempty(fstat.deltan)~= 1)
        loc = find(fstat.PJ >= 0);
        if (isempty(loc)~= 1)
            force(i) = fstat.forcen(min(loc));
            overlap(i) = fstat.deltan(min(loc));
        end
    end
end

% Plotting data where liquid bridge force is non-zero and particles are separated (negative overlap)
loc = find(force~=0 & overlap < 0);
scatter(overlap(loc(1:1:end)), force(loc(1:1:end)), 90, 'd', 'filled', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'y');

% Set x-axis and y-axis limits
xlim([-2.5e-4 0]);  % Set x-axis limits [xmin, xmax]
ylim([-3e-4 0]);        % Set y-axis limits [ymin, ymax]

% Add label to axes
xlabel('Separation [m]', 'FontSize', 12);
ylabel('Force [N]', 'FontSize', 12);

% Add a title for the plot
title('Force vs. Separation', 'FontSize', 14);