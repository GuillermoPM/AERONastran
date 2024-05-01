%% AERONastran
% Grupo 14

clear all;
close all;

%% === Files and format ===

% input and output file
entry_file = 'Alas.bdf';
output_file = 'archivo_ordenado.bdf';
matrices_file  = 'matrix.bdf';

% format line for the .bdf
format_line = '\n$...1...|...2...|...3...|...4...|...5...|...6...|...7...|...8...|...9...|..10... \n';

% Heading line for the matrix.bdf
heading_line = '';

%% === Geometry ===

% set wing span and chord
span = 2.5;
chord = 0.5;
angle = 0;

wingspecs = [span, chord, angle];

% FEM material
material = struct('E', 10000, 'dens', 2850, 'poisson', 0.3);


%% === Set analysis parameters ===
analysis = 2; % select the analyisis type (1 : flutter, 2 : divergence)
aeroparam = struct('M', 0.0, 'Q', 1500, 'AOA', 0.0174533);

nnodes = 45;
total_nodes = nnodes/5;

% set panel span and chord division
% span_division = [0.0000, 0.0625, 0.1250, 0.1875, 0.2500, 0.3167, 0.3834, ...
%                  0.4483, 0.5117, 0.56, 0.5725, 0.6317, 0.6892, 0.7450, 0.7992, ...
%                  0.8525, 0.9042, 0.9534, 1.0000];
span_division = linspace(0,1,11);
chord_division = [0.0000, 0.1091, 0.2182, 0.3886, 0.5614, 0.7273, 0.8705, 1.0000];


%% === Analyisis case ===

switch analysis
    case 1
        disp('Flutter')
    case 2
        disp('Divergence')
        create_bdf_file(output_file);

        write_solaeroelasticity_section(output_file);

        [X,Y,Z] = FEM(wingspecs, total_nodes);
        nodes = get_nodes(X, Y, Z, nnodes);
        meshplot(X,Y,Z, nodes, true);
        nodedump(nodes, format_line, output_file);

        shell_elements = gen_shell(nodes);
        shelldump(output_file, shell_elements, format_line);
        set_material(output_file, material, format_line);

        matrices(matrices_file);
        shelldump(matrices_file, shell_elements, format_line);
        set_material(matrices_file, material, format_line);
        nodedump(nodes, format_line, matrices_file);
        end_file(matrices_file);
        
        [XP, YP, ZP] = dlmpanels(wingspecs, span_division, chord_division);
        panels = get_panels(XP,YP,ZP);

        panelplot(XP,YP,ZP, panels);
        method(output_file, format_line);
        paneldivision(output_file, format_line, chord_division, span_division);
        wingstructcoupling(output_file, nodes, format_line,panels);
        set_aeros(output_file, chord, span, format_line);
        set_trim(output_file, aeroparam, format_line);
        end_file(output_file);



 
end

%% === Geometry definition ===
function paneldivision(file, format_line, chord_division, span_division)
% Sets the panel division for the DLM and writes it into the input file for
% Nastran to execute the analysis.
%   file : input file
%   chor_division : chord division array of points
%   span_division : span division array of points

% Open file in append mode
fileopen_entry = fopen(file, 'a');

% Check if open == ERROR
if fileopen_entry == -1
    error('No se pudo abrir el archivo de entrada');
end

% Write format line
fprintf(fileopen_entry, format_line);

% Write span division
fprintf(fileopen_entry, 'AEFACT  101    ');
write_array(fileopen_entry, span_division);

fprintf(fileopen_entry, format_line);
% Write chord division
fprintf(fileopen_entry, 'AEFACT  151    ');
write_array(fileopen_entry, chord_division);


% Close file
fclose(fileopen_entry);

end

%% === Doublet lattice method ===
function [X,Y,Z] = dlmpanels(wingspecs, span_division, chord_division)
    span = wingspecs(1);
    chord = wingspecs(2);
    angle = wingspecs(3);

    % Calculate the number of nodes for each section
    nnodes_section = ceil(length(span_division) / 2);

    % Generate nodes for the first section
    span_nodes1 = span * span_division(1:nnodes_section);
    chord_nodes1 = chord * chord_division - chord/2;

    % Calculate the slope for the first section
    theta1 = deg2rad(angle);
    m1 = tan(theta1);

    % Generate nodes for the second section
    % Start at the end of the first section
    span_nodes2 = span * span_division(nnodes_section:end); % Include one more node for smooth connection

    % Initialize arrays to store coordinates of all nodes
    X1 = zeros(nnodes_section, length(chord_division));
    Y1 = zeros(nnodes_section, length(chord_division));
    Z1 = zeros(nnodes_section, length(chord_division));
    X2 = zeros(nnodes_section, length(chord_division)); % Correct size initialization
    Y2 = zeros(nnodes_section, length(chord_division)); % Correct size initialization
    Z2 = zeros(nnodes_section, length(chord_division)); % Correct size initialization

    % Generate coordinates for the first section
    for i = 1:nnodes_section
        for j = 1:length(chord_division)
            X1(i, j) = chord_nodes1(j) + m1 * span_nodes1(i);
            Y1(i, j) = span_nodes1(i);
            Z1(i, j) = 0;
        end
    end
    span_nodes2 = flip(span_nodes2);
    % Generate coordinates for the second section
    for i = 1:length(span_nodes2)
        for j = 1:length(chord_division)
            X2(i, j) = chord_nodes1(j) + m1 * span_nodes1(i); % Use chord_nodes1
            Y2(i, j) = span_nodes2(i);
            Z2(i, j) = 0;
        end
    end

    X2 = flip(X2, 1);
    Y2 = flip(Y2, 1);


    X = vertcat(X1, X2);
    Y = vertcat(Y1, Y2);
    Z = vertcat(Z1, Z2);
end

function panels = get_panels(X, Y, Z)
    elim_row = size(X)/2;
    X(elim_row(1),:) = [];
    Y(elim_row(1),:) = [];
    Z(elim_row(1),:) = [];

    [rows, cols] = size(X);
    panels = struct('id', [], 'x', [], 'y', [], 'z', []);
    panel_id = 101001; % Starting ID for panels
    
    for i = 1:rows-1
        for j = 1:cols-1
            panel = struct();
            panel.id = panel_id;
            panel.x = [X(i, j), X(i+1, j), X(i+1, j+1), X(i, j+1)];
            panel.y = [Y(i, j), Y(i+1, j), Y(i+1, j+1), Y(i, j+1)];
            panel.z = [Z(i, j), Z(i+1, j), Z(i+1, j+1), Z(i, j+1)];
            
            panels(end+1) = panel;
            panel_id = panel_id + 1;
        end
    end

end



function method(filename, format_line)
    % Open the file for appending
    fid = fopen(filename, 'a');
    if fid == -1
        error('Could not open the file');
    end

    str = 'CAERO1  101001  100000                          101     151     1       +';
    str_nodes = '0.25    0.0     0.0     0.5     0.25     2.5     0.0     0.5';
    % Append the string to the file
    fprintf(fid, '$');
    fprintf(fid, format_line);
    fprintf(fid, '%s\n', str);
    fprintf(fid, '+       ');
    fprintf(fid, '%s', str_nodes);
    fprintf(fid, format_line);
    fprintf(fid, 'PAERO1  100000');
    % Close the file
    fclose(fid);
end

function wingstructcoupling(file, nodes, format_line, panels)
    % Open the file for appending
    fid = fopen(file, 'a');
    if fid == -1
        error('Could not open the file');
    end
    
    % Append the format line
    fprintf(fid, format_line);
    fprintf(fid, 'SPLINE1 101001  101001  101001  %i  101001\n',panels(end).id);
    
    fprintf(fid, 'SET1    101001');
    write_set(fid, [nodes(:).id])

    fclose(fid);
end


function set_aeros(file, chord, span, format_line)
% Used in sol 144 to give the reference coordinates system and values. The
% reference values are taken from the chord and span, and the reference
% coordinate system is the default.

fid = fopen(file, 'a');
if fid == -1
    error('Could not open the file');
end

% Write the AEROS line
aerosinfo = '\n$AEROS  ACSID   RCSID   REFC    REFB    REFS    SYMXZ   SYMXY';

fprintf(fid,aerosinfo);
fprintf(fid,format_line);
fprintf(fid, 'AEROS                   %.1f     %.1f     %.2f     \n', chord, span, span*chord);

% Close the file
fclose(fid);

end

function set_trim(file, aeroparam, format_line)
format_str = '$TRIM   ID      MACH    Q       LABEL1  UX1     LABEL2  UX2     AEQR    +\n';
format2_str = '$+      LABEL3  UX3     ...';

fid = fopen(file, 'a');
if fid == -1
    error('Could not open the file');
end
fprintf(fid,format_str);
fprintf(fid, format2_str);
fprintf(fid, format_line);
fprintf(fid, 'TRIM    6       %.1f     %.1e  ANGLEA %.6f                %.1f',aeroparam.M, aeroparam.Q, aeroparam.AOA,1);
fprintf(fid, '\nAESTAT        61  ANGLEA');


end

%% === Set FEM properties ===
function set_material(file, material, format_line)
    fid = fopen(file, 'a');
    if fid == -1
        error('Could not open the file');
    end
    
    fprintf(fid, format_line);
    fprintf(fid, 'MAT1    1       %.1e         %.1f     %.1f\n', material.E, material.poisson, material.dens);
    
    fclose(fid);

end

%% === FEM Grid generation ===
function [X,Y,Z] = FEM(wingspecs, total_nodes)
    span = wingspecs(1);
    chord = wingspecs(2);
    angle = wingspecs(3);

    
    % Calculate the number of nodes for each section
    nnodes_section = round(total_nodes / 2);
    
    % Generate nodes for the first section
    span_nodes1 = linspace(0, span/2, nnodes_section);
    chord_nodes = linspace(chord/2, -chord/2, nnodes_section);
    
    % Calculate the slope for the first section
    theta1 = deg2rad(angle);
    m1 = tan(theta1);
    
    % Generate nodes for the second section
    % Start at the end of the first section
    span_offset = span_nodes1(end); % Offset for the start of the second section
    span_nodes2 = linspace(0, span/2, nnodes_section) + span_offset; % Include one more node for smooth connection
   
   
    % Initialize arrays to store coordinates of all nodes
    X1 = zeros(nnodes_section, nnodes_section);
    Y1 = zeros(nnodes_section, nnodes_section);
    Z1 = zeros(nnodes_section, nnodes_section);

    X2 = X1; Y2 = Y1; Z2 = Z1;
    
    % Generate coordinates for the first section
    for i = 1:nnodes_section
        for j = 1:nnodes_section
            X1(i, j) = chord_nodes(j) + m1 * span_nodes1(i);
            Y1(i, j) = span_nodes1(i);
            Z1(i, j) = 0;
        end
    end
    span_nodes2 = flip(span_nodes2);
    % Generate coordinates for the second section
    for i = 1:(nnodes_section) % Include one more node for smooth connection
        for j = 1:nnodes_section
            X2(i, j) = chord_nodes(j) + m1 * span_nodes1(i);
            Y2(i, j) = span_nodes2(i);
            Z2(i, j) = 0;
        end
    end

    X2 = flip(X2,1);
    Y2 = flip(Y2,1);

    X = vertcat(X1, X2); 
    Y = vertcat(Y1, Y2);
    Z = vertcat(Z1, Z2);


end

function shell_elements = gen_shell(nodes)
% Generates the shell elements that conform the wing
% Shell elements are defined in the .bdf file with:
%   - index
%   - pshell pointer
%   - id's of the corner nodes

shell_elements = struct('id', [], 'pshell', [], 'node1', [], ...
    'node2',[], 'node3',[], 'node4',[]);
idstart = 20000;

count = 0;
% Iterate over nodes
for i = 1:length(nodes)
    % Check if y-coordinate is zero
    if nodes(i).y == 0
        % Increment count
        count = count + 1;
    end
end

nchord = count;
nspan = length(nodes)/5;

for j = 1:nchord-1
    for i = 1:nspan-1
    nodeid = 1000;
    idstart = idstart +1;

    shell_elements(end+1) = struct('id', idstart, 'pshell', 1, ...
        'node1', nodeid*j + i, ...
        'node2', nodeid*j+i+1, ...
        'node3', nodeid*(j+1) + i+1, ...
        'node4', nodeid*(j+1) +i);
    end

end
shell_elements(1) = [];

end

function nodes = get_nodes(X, Y, Z, nnodes)
    % Initialize array to store node coordinates
    node_coordinates = struct('id', [], 'x', [], 'y', [], 'z', [], 'bc', []);
    
    % Get the number of nodes and dimensions
    [nnodes_row, nnodes_col] = size(X); nnodes_row = nnodes_row-1;
    X(round(nnodes/10),:) = [];
    Y(round(nnodes/10),:) = [];
    Z(round(nnodes/10),:) = [];
    
    % Loop through the coordinates and assign boundary conditions
    for j = 1:nnodes_col
        start_id = 1000*j;
        for i = 1:nnodes_row
            x = X(i, j);
            y = Y(i, j);
            z = Z(i, j);
            id = start_id + i;
            if mod(id, 1000) == 1
            % Update the 'bc' attribute for these nodes
            bc = 123456; % Assuming 'bc' is set to 1 for these nodes
            else
                bc = [];
            end
            % Assign a unique ID for each node
            
            node_coordinates(end+1) = struct('id', id, 'x', x, 'y', y, 'z', z, 'bc', bc);
        end
    end
    
    % Remove the first empty entry
    node_coordinates(1) = [];
    
    nodes = node_coordinates;

end

function set_boundarycond(nodes)
% Sets the boundary conditions in the designed nodes.
    for i = 1:numel(nodes)
        % Check if the node id follows the pattern (e.g., 1001, 2001, 3001, ...)
        if mod(nodes(i).id, 1000) == 1
            % Update the 'bc' attribute for these nodes
            nodes(i).bc = 123456; % Assuming 'bc' is set to 1 for these nodes
        end
    end

end

%% === FEM stiffness and mass matrices ===
% To obtain the FEM stiffness and mass matrices another .bdf file is going
% to be created with the punch (.pch) file output.

function matrices(matrixfile)
    fid = fopen(matrixfile, 'w');
    if fid == -1
        error('Could not open the file');
    end
    fprintf(fid, '$ Direct Text Input for Nastran System Cell Section\n');
    fprintf(fid, '$ Direct Text Input for File Management Section\n');
    fprintf(fid, '$ Direct Text Input for Executive Control\n');
    fprintf(fid, '$ Normal Modes Analysis, Database\n');
    fprintf(fid, 'SOL 103\nCEND\n');
    fprintf(fid, '$ Direct Text Input for Global Case Control Data\n');
    fprintf(fid, 'TITLE = MSC.Nastran job created on %s\n', datestr(now, 'dd-mmm-yy at HH:MM:SS'));
    fprintf(fid, 'ECHO = NONE\n');
    fprintf(fid, 'RESVEC = NO\n');
    fprintf(fid, '$\n');  % End of section
    fprintf(fid, '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n');
    fprintf(fid, 'EXTSEOUT(STIF,DAMP,MASS,EXTID=1,DMIGPCH)\n');
    fprintf(fid, '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n');
    fprintf(fid, 'SUBCASE 1\n');
    fprintf(fid, '$ Subcase name : Default\n');
    fprintf(fid, '   SUBTITLE=Default\n');
    fprintf(fid, '   METHOD = 1\n');
    fprintf(fid, '   VECTOR(SORT1,REAL)=ALL\n');
    fprintf(fid, '   SPCFORCES(SORT1,REAL)=ALL\n');
    fprintf(fid, 'BEGIN BULK\n');
    fprintf(fid, '$ Direct Text Input for Bulk Data\n');
    fprintf(fid, 'PARAM    POST    0\n');
    fprintf(fid, 'PARAM   PRTMAXIM YES\n');
    fprintf(fid, 'EIGRL    1                       10      0\n');

    fclose(fid);
end

%% === File dump functions ===

function nodedump(nodes, format_line, output_file)
    % Dumps nodes into the .bdf file in order, adding the format line to
    % give the columns.

    % Order nodes by ID
    [~, idx] = sort([nodes.id]);
    sorted_nodes = nodes(idx);

    % Open the output file
    foutput = fopen(output_file, 'a');
    if foutput == -1
        error('Error opening the output file');
    end

    % Write the format line to the output file
    fprintf(foutput, format_line);

    % Write nodes to the output file
    for i = 1:numel(sorted_nodes)
        node = sorted_nodes(i);
        fprintf(foutput, 'GRID    %d            %.4f  %.4f  %.4f          %i\n', node.id, node.x, node.y, node.z, node.bc);
    end

    % Close the output file
    fclose(foutput);
end

function shelldump(output_file, shell_elements, format_line)
    % Open the output file for writing
    fid = fopen(output_file, 'a');
    if fid == -1
        error('Unable to open the output file.');
    end

    fprintf(fid,format_line);
    fprintf(fid,'PSHELL  1       1       .008    1               1\n\n');
    % Write shell elements to the file
    for i = 1:length(shell_elements)
        fprintf(fid, 'CQUAD4  %d   1       %d    %d    %d    %d\n', ...
            shell_elements(i).id, shell_elements(i).node1, shell_elements(i).node2, ...
            shell_elements(i).node3, shell_elements(i).node4);
    end
    
    % Close the output file
    fclose(fid);
end

%% === Plot functions ===
function meshplot(X, Y, Z, nodes, nodeplotflag)
    % Plotting the mesh
    figure;
    mesh(X, Y, Z);
    xlabel('Chord');
    ylabel('Span');
    zlabel('Z');
    title('Full Wing Mesh with Double Angle');
    grid on;

    % Plotting the node IDs if nodeplotflag is true
    if nodeplotflag
        hold on;
        for i = 1:numel(nodes)
            text(nodes(i).x, nodes(i).y, nodes(i).z, num2str(nodes(i).id), 'Color', 'red');
        end
        hold off;
    end
end

function panelplot(X, Y, Z, panels)
    % Plotting the mesh
    figure;
    mesh(X, Y, Z);
    hold on;
    
    % Plot panel IDs in the center
    for i = 1:numel(panels)
        panel = panels(i);
        x_center = mean(panel.x);
        y_center = mean(panel.y);
        z_center = mean(panel.z);
        text(x_center, y_center, z_center+0.02, num2str(panel.id), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
    end
    ylim([-0.5,3]);
    xlim([-0.5,0.5]);
    xlabel('Chord');
    ylabel('Span');
    zlabel('Z');
    title('Wing panel distribution');
    grid on;
    hold off;
end

function mesh_plot(filename)
    % Open the file for reading
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open the file');
    end
    
    % Initialize arrays to store node numbers and coordinates
    node_numbers = [];
    coordinates = [];
    
    % Read each line of the file
    tline = fgetl(fid);
    while ischar(tline)
        % Check if the line starts with 'GRID'
        if startsWith(tline, 'GRID')
            % Split the line into parts
            parts = strsplit(tline);
            % Extract node number and coordinates
            node_number = str2double(parts{2});
            x_coord = str2double(parts{3});
            y_coord = str2double(parts{4});
            z_coord = str2double(parts{5});
            % Store node number and coordinates
            node_numbers(end+1) = node_number;
            coordinates(end+1, :) = [x_coord, y_coord, z_coord];
        end
        % Read the next line
        tline = fgetl(fid);
    end
    
    % Close the file
    fclose(fid);
    
    % Plot the nodes in 3D
    figure;
    scatter3(coordinates(:, 1), coordinates(:, 2), coordinates(:, 3), 'filled', 'SizeData', 50);
    xlim([0,1])
    ylim([0,3])
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Nodes');
    grid on;
end

function plot_node_ids(nodes)
    % Extract x, y, z coordinates from the nodes struct
    x_coords = [nodes.x];
    y_coords = [nodes.y];
    z_coords = [nodes.z];

    % Plot node IDs
    scatter3(x_coords, y_coords, z_coords, 'filled');
    hold on;
    for i = 1:numel(nodes)
        text(x_coords(i), y_coords(i), z_coords(i), num2str(nodes(i).id), 'Color', 'red');
    end
    hold off;
    xlabel('Chord');
    ylabel('Span');
    zlabel('Z');
    title('Node IDs Plot');
    grid on;
end

function plot_edges(edge_nodes)
    % Open a new figure
    figure;
    
    % Plot the edges formed by edge nodes
    plot3(edge_nodes(:, 1), edge_nodes(:, 2), edge_nodes(:, 3), 'b-', 'LineWidth', 2);
    
    % Set labels and title
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Edges of the Wing');
    
    % Adjust the aspect ratio of the plot
    axis equal;
    
    % Show grid
    grid on;
end

%% === Helpers ===
function write_solaeroelasticity_section(filename)
    content = generate_solaeroelasticity_section();
    fid = fopen(filename, 'a');
    if fid == -1
        error('Unable to open file for writing');
    end
    fprintf(fid, content);
    fclose(fid);
end

function content = generate_solaeroelasticity_section()
    content = '';
    content = strcat(content, '$ ======================================================================\n');
    content = strcat(content, '$                       SOL 144 - Static Aeroelasticity\n');
    content = strcat(content, '$ ======================================================================\n');
    content = strcat(content, '$\n');
    content = strcat(content, '$ ------------------   File Management Section   -----------------------\n');
    content = strcat(content, '$ ASSIGN OUTPUT4 = ''btb.f50'',UNIT = 50,FORM=FORMATTED\n');
    content = strcat(content, '$\n');
    content = strcat(content, '$ ------------------   Executive Control Deck    -----------------------\n');
    content = strcat(content, '$\n');
    content = strcat(content, 'ID                  MSC.NASTRAN NORMAL MODES ANALYSIS\n');
    content = strcat(content, 'SOL                 144\n');
    content = strcat(content, 'TIME                5000\n');
    content = strcat(content, '$\n');
    content = strcat(content, '$ ------------------   Executive Control Deck    -----------------------\n');
    content = strcat(content, '$\n');
    content = strcat(content, '$ INCLUDE ''/homesun/dinamica/nast/2005R2/alt145_mod_interp.v2005R2''\n');
    content = strcat(content, '$\n');
    content = strcat(content, 'CEND\n');
    content = strcat(content, '$\n');
    content = strcat(content, '$ ------------------   Case Control Deck   -----------------------------\n');
    content = strcat(content, '$\n');
    content = strcat(content, 'TITLE               = BTB AEROELASTIC MODEL\n');
    content = strcat(content, 'SUBTITLE            = (SOL103)\n');
    content = strcat(content, '$\n');
    content = strcat(content, 'ECHO                = NONE\n');
    content = strcat(content, 'SEALL               = ALL\n');
    content = strcat(content, 'TRIM                = 6\n');
    content = strcat(content, 'RESVEC              = NO\n');
    content = strcat(content, 'DISPLACEMENT        = ALL\n');
    content = strcat(content, 'DIVERG              = 1\n');
    content = strcat(content, 'CMETHOD             = 10\n');
    content = strcat(content, '$\n');
    content = strcat(content, 'PARAM   POST    -1\n');
    content = strcat(content, 'PARAM   GRDPNT  0\n');
    content = strcat(content, 'PARAM   OPPHIPA 1\n');
    content = strcat(content, '$\n');
    content = strcat(content, 'DIVERG   1       2        0.0\n');
    content = strcat(content, 'EIGC     10      CLAN                                     2\n');
end

function create_bdf_file(filename)
    % Add the .bdf extension to the filename if not already present
    if ~endsWith(filename, '.bdf')
        filename = strcat(filename, '.bdf');
    end
    
    % Open the file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file for writing.');
    end
        
    % Close the file
    fclose(fid);
    
    disp(['File "', filename, '" created successfully.']);
end


function end_file(output_file)
    fid = fopen(output_file, 'a');
    if fid == -1
        error('Unable to open the output file.');
    end

    fprintf(fid, '\n$\nENDDATA');
    fclose(fid);
end

function write_array(fid, array)
    % Write the first line with 7 slots
    for i = 1:min(7, numel(array))
        fprintf(fid, ' %.4f ', array(i));
    end
    fprintf(fid,' +\n+      ');
    % % Print lines with an 8-slot limit
    for i = 8:numel(array)
        if mod(i, 8) == 0 && i ~= 8
            fprintf(fid, ' +\n+      ');
        end
        fprintf(fid, ' %.4f ', array(i));
    end
    fprintf(fid, '\n');
end

function write_set(fid, array)
    % Define the maximum number of elements per line
    
    % Write the first line with 7 slots
    for i = 1:min(7, numel(array))
        fprintf(fid, '  %i  ', array(i));
    end
    fprintf(fid,'  +\n+     ');
    % % Print lines with an 8-slot limit
    for i = 8:numel(array)
        if mod(i, 8) == 0 && i ~= 8
            fprintf(fid, '  +\n+     ');
        end
        fprintf(fid, '  %i  ', array(i));
    end
    fprintf(fid, '\n');


end