%% AERONastran
% Helper code to perform aeroelasticity simulations in NASTRAN.

%% === Wing class ===
classdef aeronastran < handle
    properties
        version = '2024-05-16';
        geom = wing_geom();
        format = file_format();
        model = build_model();
        oper = oper_cond();
        flutter = flutter_cond();
    end
    
    methods
        function W = aeronastran()
            % Constructor

        end
        function set_flightcond(W)
            aero_param(W)
        end
        function solve(W)
            if (W.oper.analysis == 1)
                % Solve flutter
                solve_flutter(W);
            elseif (W.oper.analysis == 2)
               % Solve divergence
                solve_divergence(W);
            else
            
            end
            runNastran(W.format.savefile);

            if (W.oper.analysis == 1)
                output_file_f06 = W.format.savefile;
                output_file_f06(end-3:end) = '.f06';
                flutter_extract(output_file_f06)
            end

        end
        % Build mesh
        function gen_fem(W)
            build_mesh(W)
        end
        % Build lattice
        function lattice(W)
            build_dlm(W)
        end
        % Ploting
        function panel_plot(W)
            panelplot(W);
        end
        function mesh_plot(W)
            meshplot(W);
        end


    end
end

%% ===== STRUCTURE ======
function S = oper_cond()
    % Operating conditions
    S.analysis = 2;
    S.h = 0;
    S.dens = [];
    S.vel = [];
    S.M = [];
    S.flightparam = struct('M', 0.0, 'Q', 1500, 'AOA', 0.0174533);

end

function flut = flutter_cond()
    flut.freq = size(2,1);
    flut.nroots = 0;
end

function F = file_format()
     % File and format
     F.format_line = '\n$...1...|...2...|...3...|...4...|...5...|...6...|...7...|...8...|...9...|..10... \n'; 
     F.heading_line = '';
     F.savefile = '';
end

function M = build_model()
    % Model
    M.chorddiv = 7;
    M.spandiv = 45;
    M.panspandiv = 0;
    M.panchorddiv = 0;
    M.panspan = 21;
    M.panchord = 7;

    M.xpanelcoords = 0;
    M.ypanelcoords = 0;
    M.zpanelcoords = 0;
    
    M.xnodecoords = 0;
    M.ynodecoords = 0;
    M.znodecoords = 0;
    M.nodeplotflag = true;
    M.shell = [];
    M.nodes = [];
    M.panels = [];

end

function S = wing_geom()
    % Wing geometry
    S.angle = 0;
    S.chord = 0.5;
    S.span = 2.5;
    S.angleflag = 0;
    S.material = struct('E', 7.E10, 'dens', 2750, 'poisson', 0.3);

end

%% === Main ===
function build_mesh(W)
    % Generates mesh for FEM
    % Input:
    %   - W : wing
    % Output:
    %   - W.model.nodes
    %   - W.model.shell

    % FEM model generation
    FEM(W);
    nodes = get_nodes(W);
    shell_elements = gen_shell(nodes);

    W.model.nodes = nodes;
    W.model.shell = shell_elements;
end

function build_dlm(W)
    % Generates doubblet lattice method for the wing
    % Input:
    %   - W : wing class
    % Sets:
    %   - W.model.panels : panels for dlm
    %   - W.model.panelcoords : panel coordinates 


    W.model.panspandiv= linspace(0,1,W.model.panspan);
    W.model.panchorddiv = 0.5 + 0.5*cos(linspace(pi,0,W.model.panchord));

    [XP, YP, ZP] = dlmpanels(W);
    panels = get_panels(XP,YP,ZP);

    W.model.panels = panels;
    W.model.xpanelcoords = XP;
    W.model.ypanelcoords = YP;
    W.model.zpanelcoords = ZP;

end

function solve_flutter(W)
    % Generates and solves bdf file for sol 145 (flutter).
    % Input:
    %   - W : wing

    disp('Flutter')
    create_bdf_file(W);
    write_flutter_analysis_file(W);

    nodedump(W);
    shelldump(W);
    set_material(W);

    method(W);

    write_paneldivision(W);
    wingstructcoupling(W);
    set_eigen(W);
    set_flightcond(W);

    set_flutter(W);
    set_mkaero(W);
    set_flfacts(W);
    end_file(W.format.savefile);


end

function solve_divergence(W)
    % Generates and solves bdf file for sol 144 (divergence).
    % Input:
    %   - W : wing

    disp('Divergence')
    create_bdf_file(W);
    write_divergence_analysis_file(W);

    nodedump(W);
    shelldump(W);
    set_material(W);

    method(W);
    write_paneldivision(W);
    wingstructcoupling(W);
    set_aeros(W);
    set_trim(W);
    end_file(W.format.savefile);

end

function fligthcond(W)
    h = W.oper.h;
end


%% === Geometry definition ===
function write_paneldivision(W)
% Sets the panel division for the DLM and writes it into the input file for
% Nastran to execute the analyisis.
% 
% Input:
%   - W : aeronastran wing
% Output:
%   - 

file = W.format.savefile;
format_line = W.format.format_line;
chord_division = W.model.panchorddiv;
span_division = W.model.panspandiv;

% Open file in append mode
fileopen_entry = fopen(file, 'a');

% Check if open == ERROR
if fileopen_entry == -1
    error('Cannot open entry file');
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
function [X,Y,Z] = dlmpanels(W)
    % Double lattice panel generation
    % Input
    %   - W : wing
    % Output
    %   - [X,Y,Z] : coordinates of the panel nodes 

    span = W.geom.span;
    chord = W.geom.chord;
    angle = W.geom.angle;
    span_division = W.model.panspandiv;
    chord_division = W.model.panchorddiv;

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

function method(W)
    % filename, format_line,angleflag, nodes, spandiv, panels
    % Sets the CAERO method for the .bdf
    % Input
    %   - filename : bdf file
    %   - format_line : division line
    %   - angleflag : checks if wing has some angle
    %   - nodes : node array
    %   - spandiv : span division
    
    filename = W.format.savefile;
    format_line = W.format.format_line;
    angle = W.geom.angle;
    nodes = W.model.nodes;
    spandiv = W.model.spandiv;
    panels = W.model.panels;

    % Open the file for appending
    fid = fopen(filename, 'a');
    if fid == -1
        error('Could not open the file');
    end
    index = round(spandiv/2);

    if angle == 0
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

function wingstructcoupling(W)
    % Open the file for appending wingstructcoupling.
    % Input:
    %   - W : aeronastran wing

    file = W.format.savefile;
    nodes = W.model.nodes;
    format_line = W.format.format_line;
    panels = W.model.panels;
    angle = W.geom.angle;

    fid = fopen(file, 'a');
    if fid == -1
        error('Could not open the file');
    end
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
    if angle == 0
        % Append the format line
        fprintf(fid, 'SPLINE1 50001   40001   40001   %i   901\n',panels(end).id);
        fprintf(fid, 'SET1    901     ');
        write_set(fid, [nodes(:).id])

    else
        fprintf(fid, 'SPLINE1 50001   40001   40001   %i   901\n',panels(end).id);
        fprintf(fid, 'SPLINE1 50002   %i   %i   %i   902',panels(end).id+1,panels(end).id+1, mod(panels(end).id,1000)*2+40000);
        fprintf(fid, '\nSET1    901     ');
        write_set(fid, wing1)
        fprintf(fid, '\nSET1    902     ');
        write_set(fid, wing2)
    end
   

    fclose(fid);
end

%% === Solving divergence ===
function set_aeros(W)
% Used in sol 144 to give the reference coordinates system and values. The
% reference values are taken from the chord and span, and the reference
% coordinate system is the default.

file = W.format.savefile;
chord = W.geom.chord;
span = W.geom.span;
format_line = W.format.format_line;


fid = fopen(file, 'a');
if fid == -1
    error('Could not open the file');
end

% Write the AEROS line
aerosinfo = '\n$AEROS  ACSID   RCSID   REFC    REFB    REFS    SYMXZ   SYMXY';

fprintf(fid,aerosinfo);
fprintf(fid,format_line);
fprintf(fid, 'AEROS                   %.1f     %.1f     %.2f    +1\n', span, chord, span*chord);

% Close the file
fclose(fid);

end

function set_trim(W)
    % Sets TRIM entry for static analysis.
    % Input:
    %   - W : aeronastran wing

    file = W.format.savefile;
    aeroparam = W.oper.flightparam;
    format_line = W.format.format_line;

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
function set_mkaero(W)
    % Sets MKAERO1 entry for flutter analysis.
    % Input:
    %   - outputfile : .bdf file
    %   - format_line

    outputfile = W.format.savefile;
    format_line = W.format.format_line;
    M = W.oper.M(1);

    fid = fopen(outputfile, 'a');
    if fid == -1
        error('Could not open the file');
    end
    fprintf(fid, format_line);
    fprintf(fid, 'MKAERO1 %.4f\n', M);
    fprintf(fid, '+       1.0000  0.9000  0.8000  0.7000  0.6000  0.5000  0.4000  0.3000\n');
    fprintf(fid, 'MKAERO1 %.4f\n', M);
    fprintf(fid, '+       0.2000  0.1000  0.0100');
    fclose(fid);
end

function set_eigen(W)
    % Sets eigenvalue extraction method and prints it in the .bdf file.
    % Inputs:
    %   - outputfile : .bdf file
    %   - format_line
    %   - V1 : lower frequency limit
    %   - V2 : higher frequency limit
    %   - nroots : number of modes

    outputfile = W.format.savefile;
    format_line = W.format.format_line;
    V1 = W.flutter.freq(1);
    V2 = W.flutter.freq(2);
    nroots = W.flutter.nroots;
    
    fid = fopen(outputfile, 'a');
    if fid == -1
        error('Could not open the file');
    end
    fprintf(fid, '$\n$ Modal Analysis\n$ ');
    fprintf(fid, format_line);
    fprintf(fid, '$          SID     V1      V2      ND    MSGLVL  MAXSET  SHFSCL   NORM\n');
    fprintf(fid, 'EIGRL   200     %.2f    %.2f   %i                               MAX', V1,V2, nroots);
    fclose(fid);
end  

function set_flightcond(W)
    % Sets flight conditions
    % Input:
    %   - outputfile : .bdf file
    %   - format_line
    %   - M : Mach number
    %   - Q : dynamic pressure
    
    outputfile = W.format.savefile;
    format_line = W.format.format_line;
    M = W.oper.M(1);
    Q = 1500;

    fid = fopen(outputfile, 'a');
    if fid == -1
        error('Could not open the file');
    end
    fprintf(fid, format_line);
    fprintf(fid, 'PARAM   MACH    %.4f\n', M);
    fprintf(fid, 'PARAM   Q       %.3f', Q);
    fprintf(fid, format_line);
    fprintf(fid, '$               VELOCITY  REFC  RHOREF  SIMXZ\n');
    fprintf(fid, 'AERO            9999.00   0.5   1.225   +1\n');
    fclose(fid);
end

function set_flutter(W)
    % Sets flutter method and parameters
    % Input:
    %   - outputfile : .bdf file
    %   - format_line

    outputfile = W.format.savefile;
    format_line = W.format.format_line;

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

function set_flfacts(W)
    % Sets FLFACTS entry
    % Input:
    %   - 

    outputfile = W.format.savefile;
    format_line = W.format.format_line;
    dens = W.oper.dens;
    Mach = W.oper.M;
    Vel = W.oper.vel;

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
function set_material(W)
    % Sets material entry in the .bdf file
    % Input:
    %   - W : aeronastran wing

    format_line = W.format.format_line;
    material = W.geom.material;
    file = W.format.savefile;

    fid = fopen(file, 'a');
    if fid == -1
        error('Could not open the file');
    end
    
    fprintf(fid, format_line);
    fprintf(fid, '$          MID      E       G      NU      RHO      A     TREF     GE\n');
    fprintf(fid, 'MAT1    1       %.1e         %.1f     %.1f                  0.03\n', material.E, material.poisson, material.dens);
    fclose(fid);

end

%% === FEM Grid generation ===
function FEM(W)
    % Generates FEM grid and the X,Y,Z coordinates of the nodes
    % Input:
    %   - W : aeronastran wing
    % Output:
    %   - W.model.nodecoords : node coordinates (X,Y,Z) [3 x n]

    % get wing specs
    span = W.geom.span; chord = W.geom.chord; angle = W.geom.angle;

    % generate nodes for the first section
    span_nodes1 = linspace(0, span/2, round(W.model.spandiv/2));
    chord_nodes = linspace(chord/2, -chord/2, W.model.chorddiv);
    
    % calculate the slope for the first section
    theta1 = deg2rad(angle);
    m1 = tan(theta1);
    
    % generate nodes for the second section
    % start at the end of the first section
    span_offset = span_nodes1(end); % Offset for the start of the second section
    span_nodes2 = linspace(0, span/2, round(W.model.spandiv/2)) + span_offset; % Include one more node for smooth connection
    span_nodes2(1) = [];
   
    % Initialize arrays to store coordinates of all nodes
    X1 = zeros(round(W.model.spandiv/2), W.model.chorddiv);
    Y1 = zeros(round(W.model.spandiv/2), W.model.chorddiv);
    Z1 = zeros(round(W.model.spandiv/2), W.model.chorddiv);

    X2 = zeros(round(W.model.spandiv/2)-1, W.model.chorddiv);
    Y2 = zeros(round(W.model.spandiv/2)-1, W.model.chorddiv);
    Z2 = zeros(round(W.model.spandiv/2)-1, W.model.chorddiv);
    
    % Generate coordinates for the first section
    for i = 1:round(W.model.spandiv/2)
        for j = 1:W.model.chorddiv
            X1(i, j) = chord_nodes(j) + m1 * span_nodes1(i);
            Y1(i, j) = span_nodes1(i);
            Z1(i, j) = 0;
        end
    end
    span_nodes2 = flip(span_nodes2);
    % Generate coordinates for the second section
    for i = 1:round(W.model.spandiv/2)-1
        for j = 1:W.model.chorddiv
            X2(i, j) = chord_nodes(j) + m1 * span_nodes1(i);
            Y2(i, j) = span_nodes2(i);
            Z2(i, j) = 0;
        end
    end

    X2 = flip(X2,1);
    Y2 = flip(Y2,1);

    % X,Y,Z => [spandiv, chordiv]
    W.model.xnodecoords = vertcat(X1, X2); 
    W.model.ynodecoords = vertcat(Y1, Y2);
    W.model.znodecoords = vertcat(Z1, Z2);


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

function nodes = get_nodes(W)
    % Get node struct from node coordinates
    % Input:
    %   - W : aeronastran wing
    % Output
    %   - W.model.nodes : FEM nodes

    % Initialize array to store node coordinates
    node_coordinates = struct('id', [], 'x', [], 'y', [], 'z', [], 'bc', [],'wing',[]);
    
    X = W.model.xnodecoords;
    Y = W.model.ynodecoords;
    Z = W.model.znodecoords;

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
    fprintf(fid, 'EIGRL    1                       3       0\n');

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
            nModes = str2double(tline(18:19));
        end
        nLines = nLines + 1;
    end
    
    frewind(fid);
    iModeBreak = 0;
    modeCounter = 1; % Initialize mode counter
    
    while true
        tline = fgetl(fid);
        if ~ischar(tline), break, end;
        values = strsplit(tline);
        if length(tline) >= 12 && strcmp(tline(8:12), 'POINT')
            nPointTables = nPointTables + 1;
            iMode = str2double(tline(18:19));
            for j = 1:3
                fgetl(fid);
            end
            tline(1:1) = ' ';
            if iMode ~= iModeBreak
                iLine = 0;
                iModeBreak = iMode;
            end
            while tline(1:1) == ' '
                tline = fgetl(fid);
                if strcmp(tline(1:1), '1')
                    break;
                end
                if strcmp(tline(1:4), ' ***') || strcmp(tline(1:3), '***')
                    break;
                end
                iLine = iLine + 1;
                INVk(iLine, iMode) = str2double(tline(17:28));
                DENS(iLine, iMode) = str2double(tline(31:41));
                VELO(iLine, iMode) = str2double(tline(59:69));
                EAS(iLine, iMode) = VELO(iLine, iMode) * sqrt(DENS(iLine, iMode) / 1.225);
                DAMP(iLine, iMode) = str2double(tline(72:83));
                FREQ(iLine, iMode) = str2double(tline(87:97));
                REVA(iLine, iMode) = str2double(tline(100:111));
                IMVA(iLine, iMode) = str2double(tline(114:125));
            end
        end
    end
    
    fclose(fid);
    figure('Position', [50 50 650 600]);
    
    subplot(2, 1, 1);
    hold on;
    for iMode = 1:nModes
        plot(EAS(:, iMode), DAMP(:, iMode), '.-');
    end
    grid;
    ylim([-0.10 0.5]);
    ylabel('Damping [g]');
    xlabel('Flight Speed [KTAS]');
    legend(arrayfun(@(x) ['Mode ' num2str(x)], 1:nModes, 'UniformOutput', false));
    
    subplot(2, 1, 2);
    hold on;
    for iMode = 1:nModes
        plot(EAS(:, iMode), FREQ(:, iMode), '.-');
    end
    ylim([0 50]);
    ylabel('Frequency [Hz]');
    legend(arrayfun(@(x) ['Mode ' num2str(x)], 1:nModes, 'UniformOutput', false));
    grid;
end



function nodedump(W)
    % nodes, format_line, output_file
    % Dumps nodes into the .bdf file in order, adding the format line to
    % give the columns.

    % Order nodes by ID
    nodes = W.model.nodes;
    format_line = W.format.format_line;
    output_file = W.format.savefile;


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

function shelldump(W)
    % output_file, shell_elements, format_line
    % Open the output file for writing

    output_file = W.format.savefile;
    shell_elements = W.model.shell;
    format_line = W.format.format_line;

    fid = fopen(output_file, 'a');
    if fid == -1
        error('Unable to open the output file.');
    end

    fprintf(fid,format_line);
    fprintf(fid,'$          PID     MID1    T      MID2             MID3    TS/T    NSM\n');
    fprintf(fid,'PSHELL  1       1       .008     1               1\n');
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

function meshplot(W)
    % Mesh plot
    % Input:
    %   - W : aeronastran wing
    %   - nodeflag : whether to plot the node id or not

    X = W.model.xnodecoords;
    Y = W.model.ynodecoords;
    Z = W.model.znodecoords;
    nodes = W.model.nodes;

    
    % Plotting the mesh
    figure;

    mesh(X, Y, Z, 'EdgeColor', 'blue', 'FaceAlpha',0);
    xlabel('Chord');
    ylabel('Span');
    zlabel('Z');
    title('Wing structural mesh');
    grid on;

    % Plotting the node IDs if nodeplotflag is true
    if W.model.nodeplotflag
        hold on;
        for i = 1:numel(nodes)
            plot3(nodes(i).x, nodes(i).y, nodes(i).z, 'Marker', 'o', 'MarkerFaceColor', 'red');
            text(nodes(i).x, nodes(i).y, nodes(i).z-0.05, num2str(nodes(i).id), 'Color', 'red', 'FontSize',8);
        end
        hold off;
    end
end

function panelplot(W)
    % Plotting the mesh
    figure;

    X = W.model.xpanelcoords;
    Y = W.model.ypanelcoords;
    Z = W.model.zpanelcoords;
    panels = W.model.panels;

    mesh(X, Y, Z,'FaceAlpha',0);
    hold on;
    
    % Plot panel IDs in the center
    for i = 1:numel(panels)
        panel = panels(i);
        x_center = mean(panel.x);
        y_center = mean(panel.y);
        z_center = mean(panel.z);
        text(x_center, y_center, z_center-0.02, num2str(panel.id-40000), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red', 'FontSize',5);
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

%% === Run NASTRAN ===
function runNastran(input_file)
    % build nastran path
    nastran_path = '"C:/Program Files/MSC.Software/NaPa_SE/20231/Nastran/bin/nastranw.exe"';
        
    % execute nastran and input bdf file
    [~,status] = system([nastran_path ' input.bdf < ' input_file]);

    resultsfile = char(input_file);
    resultsfile(end-3:end) = '.f06';
        
    progress = 0;
    
    while true
        pause(1.5);
        progress = progress + 10;
        if progress > 100
            progress = 0;  
            
            fprintf('Generating: [%s%s] %d%%\r', repmat('#', 1, progress/10), repmat('-', 1, (100-progress)/10), progress);    
            if exist(resultsfile, 'file')
                break;
            end
        end
        % verify status
    end
    if status == 0
        disp('NASTRAN executed.');
    else
        disp('EXE error check path.');
    end
    
end

%% === Helpers ===

function aero_param(W)
    % Input:
    %   - h : altitude
    % Output:
    %   - dens : calculated density at certain height
    %   - vel : velocity
    %   - M : Mach number
    h = W.oper.h;
    
    dens_rel = (1-22.558*10^(-6)*h)^4.2559; % calculate density
    T = 288.15 - 6.5*h/1000; % calculate temperature
    
    R_air = 287;
    a = sqrt(1.4*R_air*T);  % sound speed

    % set oper_cond in wing class
    W.oper.dens = dens_rel*ones(1,60);
    W.oper.vel = round(linspace(10,150,60),2);
    W.oper.M = W.oper.vel/a;

end

function write_divergence_analysis_file(W)
    filename = W.format.savefile;
    fid = fopen(filename, 'a');
    if fid == -1
        error('Unable to open file for writing');
    end
    
    fprintf(fid, '$ ======================================================================\n');
    fprintf(fid, '$                       SOL 144 - Static Aeroelasticity\n');
    fprintf(fid, '$ ======================================================================\n');
    fprintf(fid, '$\n');
    fprintf(fid, '$ ------------------   File Management Section   -----------------------\n');
    fprintf(fid, '$ ASSIGN OUTPUT4 = ''btb.f50'',UNIT = 50,FORM=FORMATTED\n');
    fprintf(fid, '$\n');
    fprintf(fid, '$ ------------------   Executive Control Deck    -----------------------\n');
    fprintf(fid, '$\n');
    fprintf(fid, 'ID                  MSC.NASTRAN NORMAL MODES ANALYSIS\n');
    fprintf(fid, 'SOL                 144\n');
    fprintf(fid, 'TIME                5000\n');
    fprintf(fid, '$\n');
    fprintf(fid, '$ ------------------   Executive Control Deck    -----------------------\n');
    fprintf(fid, '$\n');
    fprintf(fid, '$ INCLUDE ''/homesun/dinamica/nast/2005R2/alt145_mod_interp.v2005R2''\n');
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
    fprintf(fid, 'TRIM                = 6\n');
    fprintf(fid, 'RESVEC              = NO\n');
    fprintf(fid, 'DISPLACEMENT        = ALL\n');
    fprintf(fid, 'DIVERG              = 1\n');
    fprintf(fid, 'CMETHOD             = 10\n');
    fprintf(fid, '$\n');
    fprintf(fid, 'PARAM   POST    -1\n');
    fprintf(fid, 'PARAM   GRDPNT  0\n');
    fprintf(fid, 'PARAM   OPPHIPA 1\n');
    fprintf(fid, '$\n');
    fprintf(fid, 'DIVERG   1       5        0.0\n');
    fprintf(fid, 'EIGC     10      CLAN                                     5\n');

    fclose(fid);
end

function write_flutter_analysis_file(W)
    
    filename = W.format.savefile;
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


function create_bdf_file(W)
% wingspecs, analysis
    % Extract angle from wingspecs
    angle = W.geom.angle;

    if W.oper.analysis == 1
    
        str = 'flutter';
    else
        str = 'div';
    end
    
    % Convert angle to a string
    if angle >= 0
        angle_str = num2str(angle);
    else
        angle_str = ['neg_', num2str(abs(angle))];
    end

    % Construct the filename
    filename = ['AERONastran_angle_', angle_str, '_', str, '.bdf'];
    % Open the file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file for writing.');
    end
    % Close the file
    fclose(fid);

    W.format.savefile = char(filename);
    
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
    % Determine the maximum number of digits in the node IDs
    max_digits = numel(num2str(max(array)));

    % Write the first line with appropriate spacing
    for i = 1:min(7, numel(array))
        fprintf(fid, ['%-' num2str(max_digits) 'i    '], array(i));
    end
    fprintf(fid, '  +\n+       ');

    % Print subsequent lines with adjusted spacing based on the number of digits
    for i = 8:numel(array)
        if mod(i, 8) == 0 && i ~= 8
            fprintf(fid, '  +\n+       ');
        end
        fprintf(fid, ['%-' num2str(max_digits) 'i    '], array(i));
    end
    fprintf(fid, '\n');
end

