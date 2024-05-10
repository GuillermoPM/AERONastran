%% AERONastran
% Helper code to perform aeroelasticity simulations in NASTRAN.

clear all;
close all;

%% === Files and format ===

% input and output file

f06_file = 'archivo_ordenado.f06';

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
material = struct('E', 7.E10, 'dens', 2750, 'poisson', 0.3);



%% === Set analysis parameters ===
analysis = 2; % select the analyisis type (1 : flutter, 2 : divergence, 3 : stiffness and mass matrix)
flightparam = struct('M', 0.0, 'Q', 1500, 'AOA', 0.0174533);
meshdiv = [9, 35];

% set panel span and chord division
panspandiv = linspace(0,1,11);
panchorddiv = [0.0000, 0.1091, 0.2182, 0.3886, 0.5614, 0.7273, 0.8705, 1.0000];

dens = linspace(0.3,1.225,20);
M = 0.2*ones(1,length(dens));
vel = round(linspace(50,150,20),2);

flutterparam = [dens;M;vel];
%% DEBUG
main(flutterparam, wingspecs, material, analysis, flightparam, meshdiv, panspandiv, panchorddiv, format_line)


%% === Main ===
function main(flutterparam, wingspecs, material, analysis, flightparam, meshdiv, panspandiv, panchorddiv,format_line)
    dens = flutterparam(1,:); M = flutterparam(2,:); vel = flutterparam(3,:);
    if wingspecs(3) == 0
        angleflag = 0;
    else
        angleflag = 1;
    end
    % FEM model generation
    [X,Y,Z] = FEM(wingspecs, meshdiv(1), meshdiv(2));
    nodes = get_nodes(X, Y, Z);
    meshplot(X,Y,Z, nodes, true);
    shell_elements = gen_shell(nodes);
    
    % DLM model generation
    [XP, YP, ZP] = dlmpanels(wingspecs, panspandiv, panchorddiv);
    panels = get_panels(XP,YP,ZP);
    [spandiv, ~] = size(X);
    
    switch analysis
        case 1
            disp('Flutter')
            output_file = create_bdf_file(wingspecs, 'flutter');
            write_flutter_analysis_file(output_file);

            nodedump(nodes, format_line, output_file);
            shelldump(output_file, shell_elements, format_line);
            set_material(output_file, material, format_line);
    
            panelplot(XP,YP,ZP, panels);
            method(output_file, format_line, angleflag, nodes, spandiv, panels);
            paneldivision(output_file, format_line, panchorddiv, panspandiv);
            wingstructcoupling(output_file, nodes, format_line,panels, angleflag);
            set_eigen(output_file, format_line, 0, 15,10);
            set_flightcond(output_file, format_line, 0.2, 1500);
            set_flutter(output_file, format_line);
            set_mkaero(output_file, format_line);
            set_flfacts(output_file, dens, M, vel, format_line);
            end_file(output_file);
       
        case 2
            disp('Divergence')
            output_file = create_bdf_file(wingspecs, 'divergence');
            write_solaeroelasticity_section(output_file);
    
            nodedump(nodes, format_line, output_file);

            shelldump(output_file, shell_elements, format_line);
            set_material(output_file, material, format_line);
            
            panelplot(XP,YP,ZP, panels);
            method(output_file, format_line, angleflag, nodes, spandiv, panels);
            paneldivision(output_file, format_line, panchorddiv, panspandiv);
            wingstructcoupling(output_file, nodes, format_line,panels, angleflag);
            set_aeros(output_file, wingspecs, format_line);
            set_trim(output_file, flightparam, format_line);
            end_file(output_file);
        case 3
            disp('Matrices');
            matrices_file  = 'matrix.bdf';
            matrices(matrices_file);

    end

    runNastran(output_file);


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
    % Double lattice panel generation
    % Input
    %   - wingspecs : span, chord, angle [1x3]
    %   - span_division : span number of divisions
    %   - chord_division : chord number of divisions
    % Output
    %   - [X,Y,Z] : coordinates of the panel nodes
    
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
    % flip the matrices to fit the wing
    X2 = flip(X2, 1); 
    Y2 = flip(Y2, 1);

    % concatenate the matrices and create the full mesh
    X = vertcat(X1, X2);
    Y = vertcat(Y1, Y2);
    Z = vertcat(Z1, Z2);
end

function panels = get_panels(X, Y, Z)
    % Generates the DLM panels from the point mesh that gives the location
    % of each node. The panel generation is only for plot purpposes, as in
    % the .bdf only the division is needed.
    % Input:
    %   - [X,Y,Z] : mesh coordinates
    % Output:
    %   - panels : panel struct

    elim_row = size(X)/2;
    X(elim_row(1),:) = [];
    Y(elim_row(1),:) = [];
    Z(elim_row(1),:) = [];

    [rows, cols] = size(X);
    panels = struct('id', [], 'x', [], 'y', [], 'z', []);
    panel_id = 40001; % Starting ID for panels
    
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
    panels(1) = [];

end

function method(filename, format_line,angleflag, nodes, spandiv, panels)
    % Sets the CAERO method for the .bdf
    % Input
    %   - filename : bdf file
    %   - format_line : division line
    %   - angleflag : checks if wing has some angle
    %   - nodes : node array
    %   - spandiv : span division


    % Open the file for appending
    fid = fopen(filename, 'a');
    if fid == -1
        error('Could not open the file');
    end
    index = round(spandiv/2);

    if angleflag == 0
        str = 'CAERO1  40001   10000                           101     151     1       +';
        str_nodes = '-0.25   0.0     0.0     0.5     -0.25   2.5     0.0     0.5';
        % Append the string to the file
        fprintf(fid, '$');
        fprintf(fid, format_line);
        fprintf(fid, '%s\n', str);
        fprintf(fid, '+       ');
        fprintf(fid, '%s', str_nodes);
        fprintf(fid, format_line);
        fprintf(fid, 'PAERO1  10000');
        % Close the file
        fclose(fid);
    else
        str_w1 = 'CAERO1  40001   10000                           101     151     1       +';
       
        midxcoord = 2*nodes(1).x - nodes(index).x;
        midycoord = nodes(index).y;
        chord = 2*nodes(1).x;

        fprintf(fid, '$');
        fprintf(fid, format_line);
        fprintf(fid, '%s\n', str_w1);
        fprintf(fid, '+       ');
        fprintf(fid, '-%.2f   0.0     0.0     %.1f     -%.4f %.2f    0.0     %.1f',nodes(1).x,chord,midxcoord,midycoord,chord);
        fprintf(fid, format_line);
        fprintf(fid, 'CAERO1  %i   10000                           101     151     1       +\n',panels(end).id+1);
        fprintf(fid, '+       ');
        fprintf(fid, '-%.4f %.2f    0.0     %.1f     %.2f   %.4f  0.0     %.1f',midxcoord,midycoord, chord, nodes(end).x, nodes(end).y, chord);
        fprintf(fid, format_line);
        fprintf(fid, 'PAERO1  10000');


    end


end

function wingstructcoupling(file, nodes, format_line, panels, angleflag)
    % Open the file for appending wingstructcoupling.
    % Input:
    %   - file
    %   - nodes : node struct array
    %   - format_line
    %   - panels : panel struct array
    %   - angleflag : flag that activates if wing has angle
    fid = fopen(file, 'a');
    if fid == -1
        error('Could not open the file');
    end
    pindex = length(panels)/2;
    lastid = num2str(nodes(end).id);
    nindex = round(str2double(lastid(end-1:end))/2);
    wing1 = [];
    wing2 = [];

    for j = 1:numel(nodes)
        if nodes(j).wing == 1
            wing1(end+1) = nodes(j).id;
        elseif nodes(j).wing == 2
            wing2(end+1) = nodes(j).id;
        end

    end

    fprintf(fid, '$          EID    CAERO    BOX1    BOX2    SETG');
    fprintf(fid, format_line);
    if angleflag == 0
        % Append the format line
        fprintf(fid, 'SPLINE1 50001   40001   40001   %i   901\n',panels(end).id);
        fprintf(fid, 'SET1    901   ');
        write_set(fid, [nodes(:).id])

    else
        fprintf(fid, 'SPLINE1 50001   40001   40001   %i   901\n',panels(end).id);
        fprintf(fid, 'SPLINE1 50002   %i   %i   %i   902',panels(end).id+1,panels(end).id+1, mod(panels(end).id,100)*2+40000);
        fprintf(fid, '\nSET1    901   ');
        write_set(fid, wing1)
        fprintf(fid, '\nSET1    902   ');
        write_set(fid, wing2)
    end
   

    fclose(fid);
end

%% === Solving divergence ===
function set_aeros(file, wingspecs, format_line)
% Used in sol 144 to give the reference coordinates system and values. The
% reference values are taken from the chord and span, and the reference
% coordinate system is the default.

chord = wingspecs(1); span = wingspecs(2);

fid = fopen(file, 'a');
if fid == -1
    error('Could not open the file');
end

% Write the AEROS line
aerosinfo = '\n$AEROS  ACSID   RCSID   REFC    REFB    REFS    SYMXZ   SYMXY';

fprintf(fid,aerosinfo);
fprintf(fid,format_line);
fprintf(fid, 'AEROS                   %.1f     %.1f     %.2f    +1\n', chord, span, span*chord);

% Close the file
fclose(fid);

end

function set_trim(file, aeroparam, format_line)
    % Sets TRIM entry for static analysis.
    % Input:
    %   - file : .bdf file
    %   - aeroparam : trim conditions
    %   - format_line
    
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
%% === Solving flutter ===

function set_mkaero(outputfile, format_line)
    % Sets MKAERO1 entry for flutter analysis.
    % Input:
    %   - outputfile : .bdf file
    %   - format_line

    fid = fopen(outputfile, 'a');
    if fid == -1
        error('Could not open the file');
    end
    fprintf(fid, format_line);
    fprintf(fid, 'MKAERO1 0.2000\n');
    fprintf(fid, '+       1.0000  0.9000  0.8000  0.7000  0.6000  0.5000  0.4000  0.3000\n');
    fprintf(fid, 'MKAERO1 0.2000\n');
    fprintf(fid, '+       0.2000  0.1000  0.0100');
    fclose(fid);
end

function set_eigen(outputfile, format_line, V1, V2, nroots)
    % Sets eigenvalue extraction method and prints it in the .bdf file.
    % Inputs:
    %   - outputfile : .bdf file
    %   - format_line
    %   - V1 : lower frequency limit
    %   - V2 : higher frequency limit
    %   - nroots : number of modes
    
    fid = fopen(outputfile, 'a');
    if fid == -1
        error('Could not open the file');
    end
    fprintf(fid, '$\n$ Modal Analysis\n$ ');
    fprintf(fid, format_line);
    fprintf(fid, '$          SID     V1      V2      ND    MSGLVL  MAXSET  SHFSCL   NORM\n');
    fprintf(fid, 'EIGRL   200     %.2f    %.2f   %i                              MAX', V1,V2, nroots);
    fclose(fid);
end  

function set_flightcond(outputfile, format_line, M, Q)
    % Sets flight conditions
    % Input:
    %   - outputfile : .bdf file
    %   - format_line
    %   - M : Mach number
    %   - Q : dynamic pressure

    fid = fopen(outputfile, 'a');
    if fid == -1
        error('Could not open the file');
    end
    fprintf(fid, format_line);
    fprintf(fid, 'PARAM   MACH    %.4f\n', M);
    fprintf(fid, 'PARAM   Q       %.3f', Q);
    fprintf(fid, format_line);
    fprintf(fid, '$               VELOCITY  REFC  RHOREF  SIMXZ\n');
    fprintf(fid, 'AERO            1.00      0.5   1.225   +1\n');
    fclose(fid);
end

function set_flutter(outputfile, format_line)
    % Sets flutter method and parameters
    % Input:
    %   - outputfile : .bdf file
    %   - format_line

    fid = fopen(outputfile, 'a');
    if fid == -1
        error('Could not open the file');
    end
    fprintf(fid, '$\n$ Flutter Analysis\n$\n');
    fprintf(fid, 'PARAM   VREF    0.514444');
    fprintf(fid, format_line);
    fprintf(fid, '$               METHOD  DENS    MACH    VEL     IMETH\n');
    fprintf(fid, 'FLUTTER 300     PKNL    1       2       3       L');
    fclose(fid);
    

end

function set_flfacts(outputfile, dens, Mach, Vel, format_line)
    % Sets FLFACTS entry
    % Input:
    %   - 
    fid = fopen(outputfile, 'a');
    if fid == -1
        error('Could not open the file');
    end
    fprintf(fid,'\n$    SID     F1      F2      F3      F4      F5      F6      F7\n');
    fprintf(fid, format_line);
    fprintf(fid,'FLFACT  1      ');
    write_array(fid, dens);
    fprintf(fid,'FLFACT  2      ');
    write_array(fid, Mach);
    fprintf(fid,'FLFACT  3      ');
    write_flutter_array(fid, Vel);
    fclose(fid);

end

%% === Set FEM properties ===
function set_material(file, material, format_line)
    % Sets material entry in the .bdf file
    % Input:
    %   - file
    %   - material : material struct [E, dens, poisson]
    fid = fopen(file, 'a');
    if fid == -1
        error('Could not open the file');
    end
    
    fprintf(fid, format_line);
    fprintf(fid, '$          MID      E       G      NU      RHO      A     TREF     GE\n');
    fprintf(fid, 'MAT1    1       %.1e         %.1f     %.1f                  0.05\n', material.E, material.poisson, material.dens);
    fclose(fid);

end

%% === FEM Grid generation ===
function [X,Y,Z] = FEM(wingspecs, chordnodes,spannodes)
    % Generates FEM grid and outputs the X,Y,Z coordinates of the nodes
    % Input:
    %   - wingspecs : wing specs array [1x3]
    %   - chordnodes : chord node division
    %   - spannodes : span node division
    % Output:
    %   - [X,Y,Z] : node coordinates mesh

    % get wing specs
    span = wingspecs(1); chord = wingspecs(2); angle = wingspecs(3);

    % generate nodes for the first section
    span_nodes1 = linspace(0, span/2, round(spannodes/2));
    chord_nodes = linspace(chord/2, -chord/2, chordnodes);
    
    % calculate the slope for the first section
    theta1 = deg2rad(angle);
    m1 = tan(theta1);
    
    % generate nodes for the second section
    % start at the end of the first section
    span_offset = span_nodes1(end); % Offset for the start of the second section
    span_nodes2 = linspace(0, span/2, round(spannodes/2)) + span_offset; % Include one more node for smooth connection
    span_nodes2(1) = [];
   
    % Initialize arrays to store coordinates of all nodes
    X1 = zeros(round(spannodes/2), chordnodes);
    Y1 = zeros(round(spannodes/2), chordnodes);
    Z1 = zeros(round(spannodes/2), chordnodes);

    X2 = zeros(round(spannodes/2)-1, chordnodes);
    Y2 = zeros(round(spannodes/2)-1, chordnodes);
    Z2 = zeros(round(spannodes/2)-1, chordnodes);
    
    % Generate coordinates for the first section
    for i = 1:round(spannodes/2)
        for j = 1:chordnodes
            X1(i, j) = chord_nodes(j) + m1 * span_nodes1(i);
            Y1(i, j) = span_nodes1(i);
            Z1(i, j) = 0;
        end
    end
    span_nodes2 = flip(span_nodes2);
    % Generate coordinates for the second section
    for i = 1:round(spannodes/2)-1
        for j = 1:chordnodes
            X2(i, j) = chord_nodes(j) + m1 * span_nodes1(i);
            Y2(i, j) = span_nodes2(i);
            Z2(i, j) = 0;
        end
    end

    X2 = flip(X2,1);
    Y2 = flip(Y2,1);

    % X,Y,Z => [spandiv, chordiv]
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
nspan = length(nodes)/nchord;

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

function nodes = get_nodes(X, Y, Z)
    % Initialize array to store node coordinates
    node_coordinates = struct('id', [], 'x', [], 'y', [], 'z', [], 'bc', [],'wing',[]);
    
    % Get the number of nodes and dimensions
    [nnodes_row, nnodes_col] = size(X);
    
    % Loop through the coordinates and assign boundary conditions
    for j = 1:nnodes_col
        start_id = 1000*j;
        for i = 1:nnodes_row
            x = X(i, j);
            y = Y(i, j);
            z = Z(i, j);
            id = start_id + i; % Assign a unique ID for each node
            if mod(id, 1000) == 1
            % Update the 'bc' attribute for these nodes
                bc = 123456; % Assuming 'bc' is set to 1 for these nodes
            else
                bc = [];
            end
            if mod(id, 100) <= round(nnodes_row/2)
               wing = 1;
            else
               wing = 2;
            end
            
            node_coordinates(end+1) = struct('id', id, 'x', x, 'y', y, 'z', z, 'bc', bc,'wing',wing);
        end
    end
    
    % Remove the first empty entry
    node_coordinates(1) = [];
    
    nodes = node_coordinates;

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

%% === File comprehension functions ===
function flutter_extract(filename)
    % Extracts the flutter results from the .f06 file.
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open the file');
    end
    
    nPointTables = 0; % Initialize nPointTables
    nModes = 0; % Initialize nModes
    nLines = 0; % Initialize nLines
    
    while true
        tline = fgetl(fid);
        
        if ~ischar(tline), break, end;
        
        if length(tline) >= 12 && strcmp(tline(8:12), 'POINT')
            nPointTables = nPointTables + 1;
            nModes = str2num(tline(18:19));
        end
        nLines = nLines + 1;
    end
    
    frewind(fid);
    iModeBreak = 0;
    
    while true
        tline = fgetl(fid);
        if ~ischar(tline), break, end;
        values = strsplit(tline);
        if length(tline) >= 12 && strcmp(tline(8:12), 'POINT')
            nPointTables = nPointTables + 1;
            iMode = str2num(tline(18:19));
            for j = 1:3
                fgetl(fid);
            end
            tline(1:1) = ' ';
            if iMode ~= iModeBreak
                iLine = 0;
            end
            while tline(1:1) == ' '
                tline = fgetl(fid);
                if strcmp(tline(1:1), '1')
                    iModeBreak = iMode;
                    break;
                end
                if strcmp(tline(1:4), ' ***') || strcmp(tline(1:3), '***')
                    break;
                end
                iLine = iLine + 1;
                INVk(iLine, iMode) = str2num(tline(17:28));
                DENS(iLine, iMode) = str2num(tline(31:41));
                VELO(iLine, iMode) = str2num(tline(59:69));
                EAS(iLine, iMode) = VELO(iLine, iMode) * sqrt(DENS(iLine, iMode) / 1.225);
                DAMP(iLine, iMode) = str2num(tline(72:83));
                FREQ(iLine, iMode) = str2num(tline(87:97));
                REVA(iLine, iMode) = str2num(tline(100:111));
                IMVA(iLine, iMode) = str2num(tline(114:125));
            end
        end
    end
    
    fclose(fid);
    figure('Position', [50 50 650 600]);
    
    subplot(2, 1, 1);
    hold on;
    for iMode = 1:nModes
        plot(EAS(:, iMode), -DAMP(:, iMode), '.-');
    end
    grid;
    ylim([-0.10 0.5]);
    ylabel('Damping [g]');
    xlabel('Flight Speed [KTAS]');
    
    subplot(2, 1, 2);
    hold on;
    for iMode = 1:nModes
        plot(EAS(:, iMode), FREQ(:, iMode), '.-');
    end
    ylim([0 30]);
    ylabel('Frequency [Hz]');
    grid;

end


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
    fprintf(fid,'$          PID     MID1    T      MID2             MID3    TS/T    NSM\n');
    fprintf(fid,'PSHELL  1       1       .008    1               1\n');
    fprintf(fid,format_line);

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

function plot_displaced_nodes(f06_file, original_nodes)
    % Read the F06 file and extract displacement data
    fid = fopen(f06_file, 'r');
    if fid == -1
        error('Error: Could not open F06 file');
    end
    
    displacement_data = [];
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, '          ') && contains(line, ' G ')
            % Extract displacement data for each node
            displacement_values = sscanf(line, '%d %*s %f %f %f %f %f %f');
            displacement_data = [displacement_data; displacement_values'];
        end
    end
    fclose(fid);
    
    
    % Check if displacement data matches the number of nodes
    if size(displacement_data, 1) ~= numel(original_nodes)
        error('Error: Number of nodes in F06 file does not match original node data');
    end
    displcmnt = 1*displacement_data(:,2:4);
    rot = displacement_data(:,5:7);
    
    % Extract original coordinates of nodes
    original_x = [original_nodes.x];
    original_y = [original_nodes.y];
    original_z = [original_nodes.z];
    
    % Extract displacement vectors
    displacement_x = displcmnt(:, 1);
    displacement_y = displcmnt(:, 2);
    displacement_z = displcmnt(:, 3);

    % Extract rotation data
    rotation_x = rot(:, 1);
    rotation_y = rot(:, 2);
    rotation_z = rot(:, 3);

    % Calculate new coordinates after displacement
    displaced_x = original_x + displacement_x';
    displaced_y = original_y + displacement_y';
    displaced_z = original_z + displacement_z';
    
    % Apply rotations to displaced points
    rotated_displaced_points = [];
    for i = 1:size(displacement_data, 1)
        % Define rotation matrices for each axis
        Rx = [1 0 0; 0 cosd(rotation_x(i)) -sind(rotation_x(i)); 0 sind(rotation_x(i)) cosd(rotation_x(i))];
        Ry = [cosd(rotation_y(i)) 0 sind(rotation_y(i)); 0 1 0; -sind(rotation_y(i)) 0 cosd(rotation_y(i))];
        Rz = [cosd(rotation_z(i)) -sind(rotation_z(i)) 0; sind(rotation_z(i)) cosd(rotation_z(i)) 0; 0 0 1];
        
        % Apply rotations
        rotated_point = Rz * Ry * Rx * [displaced_x(i); displaced_y(i); displaced_z(i)];
        rotated_displaced_points(i, :) = rotated_point';
    end

       
    % Plot mesh using displaced nodes
    figure;
    plot3(rotated_displaced_points(:, 1), rotated_displaced_points(:, 2), rotated_displaced_points(:, 3), '-o', 'LineWidth', 1.5);    

    ylim([-0.5,6]);
    xlim([-1,1]);
    zlim([-0.1,3]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Displaced Nodes and Deformed Wing Mesh');
    
    % Show grid
    grid on;
    
    % Show plot

end

function meshplot(X, Y, Z, nodes, nodeplotflag)
    % Plotting the mesh
    figure;
    mesh(X, Y, Z, 'EdgeColor', 'blue', 'FaceAlpha',0);
    xlabel('Chord');
    ylabel('Span');
    zlabel('Z');
    title('Wing structural mesh');
    grid on;

    % Plotting the node IDs if nodeplotflag is true
    if nodeplotflag
        hold on;
        for i = 1:numel(nodes)
            plot3(nodes(i).x, nodes(i).y, nodes(i).z, 'Marker', 'o', 'MarkerFaceColor', 'red');
            text(nodes(i).x, nodes(i).y, nodes(i).z-0.05, num2str(nodes(i).id), 'Color', 'red', 'FontSize',8);
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
    ylim([-0.5,6]);
    xlim([-1,1]);
    xlabel('Chord');
    ylabel('Span');
    zlabel('Z');
    title('Wing panel distribution');
    grid on;
    hold off;
end

%% === Run NASTRAN ===
function runNastran(input_file)
    % build nastran path
    nastran_path = '"C:/Program Files/MSC.Software/NaPa_SE/20231/Nastran/bin/nastranw.exe"';
    
    % full input path
    input_path = fullfile(pwd, input_file);
    
    % execute nastran and input bdf file
    status = system([nastran_path ' input.bdf < ' input_path]);
    
    % verify status
    if status == 0
        disp('NASTRAN executed.');
    else
        disp('EXE error check path.');
    end
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

function write_flutter_analysis_file(filename)
    fid = fopen(filename, 'a');
    if fid == -1
        error('Error: Unable to open file for writing');
    end
    
    % Write the header
    fprintf(fid, '$ ======================================================================\n');
    fprintf(fid, '$                       SOL 145 - flutter Analysis \n');
    fprintf(fid, '$ ======================================================================\n');
    fprintf(fid, '$\n');
    fprintf(fid, '$ ------------------   File Management Section   -----------------------\n');
    fprintf(fid, '$ ASSIGN OUTPUT4 = ''btb.f50'',UNIT = 50,FORM=FORMATTED        \n');
    fprintf(fid, '$\n');
    fprintf(fid, '$ ------------------   Executive Control Deck    -----------------------\n');
    fprintf(fid, '$\n');
    fprintf(fid, 'ID                  MSC.NASTRAN NORMAL MODES ANALYSIS\n');
    fprintf(fid, 'SOL                 145\n');
    fprintf(fid, 'TIME                5000\n');
    fprintf(fid, '$\n');
    fprintf(fid, '$ ------------------   Executive Control Deck    -----------------------\n');
    fprintf(fid, '$\n');
    fprintf(fid, 'CEND\n');
    fprintf(fid, '$\n');
    fprintf(fid, '$ ------------------   Case Control Deck   -----------------------------\n');
    fprintf(fid, '$\n');
    fprintf(fid, 'TITLE               = BTB AEROELASTIC MODEL\n');
    fprintf(fid, 'SUBTITLE            = (SOL103)\n');
    fprintf(fid, '$\n');
    fprintf(fid, 'ECHO                = NONE\n');
    fprintf(fid, 'SEALL               = ALL\n');
    fprintf(fid, '$ SPC                 = 100\n');
    fprintf(fid, 'METHOD              = 200\n');
    fprintf(fid, 'FMETHOD             = 300\n');
    fprintf(fid, 'RESVEC              = NO\n');
    fprintf(fid, 'DISPLACEMENT        = ALL\n');
    
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

function outputfile = create_bdf_file(wingspecs, analysis)
    % Extract angle from wingspecs
    angle = wingspecs(3);
    
    % Convert angle to a string
    if angle >= 0
        angle_str = num2str(angle);
    else
        angle_str = ['neg_', num2str(abs(angle))];
    end

    % Construct the filename
    filename = ['AERONastran_angle_', angle_str, '_', analysis, '.bdf'];
    
    % Open the file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file for writing.');
    end
    outputfile = filename;
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

function write_flutter_array(fid, array)
    % Calculate the maximum number of digits in the array
    max_digits = max(floor(log10(abs(array))) + 1);
    
    % Write the first line with appropriate spacing
    for i = 1:min(7, numel(array))
        fprintf(fid, '%*.*f ', max_digits + 4, 2, array(i));
    end
    fprintf(fid, ' +\n+      ');
    
    % Print lines with an 8-slot limit and appropriate spacing
    for i = 8:numel(array)
        if mod(i, 8) == 0 && i ~= 8
            fprintf(fid, ' +\n+      ');
        end
        fprintf(fid, '%*.*f ', max_digits + 4, 2, array(i));
    end
    
    fprintf(fid, '\n');
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