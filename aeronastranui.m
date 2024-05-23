
function aeronastranui
    % Create the main figure window
    W = aeronastran;
    W.gen_fem;
    W.lattice;
   
    fig = uifigure('Name', sprintf('aeronastran (version %s)',W.version), 'Position', [100, 100, 800, 600]);
    % Geometry
    % Create labels and input fields for wingspan, chord, angle, and case
    lblWingspan = uilabel(fig, 'Text', 'Wingspan:', 'Position', [50, 500, 100, 22]);
    
    % Create a numeric input field for the wingspan
    txtWingspan = uieditfield(fig, 'numeric', 'Position', [150, 500, 100, 22]);
    
    % Set the initial value of the input field
    txtWingspan.Value = W.geom.span;

    txtWingspan.ValueChangedFcn = @(txt, event) set_span(W, txt.Value);

    lblChord = uilabel(fig, 'Text', 'Chord:', 'Position', [50, 450, 100, 22]);
    txtChord = uieditfield(fig, 'numeric', 'Position', [150, 450, 100, 22]);

    txtChord.Value = W.geom.chord;
    txtChord.ValueChangedFcn = @(txt, event) set_chord(W, txt.Value);

    lblAngle = uilabel(fig, 'Text', 'Angle:', 'Position', [50, 400, 100, 22]);
    txtAngle = uieditfield(fig, 'numeric', 'Position', [150, 400, 100, 22]);

    txtAngle.Value = W.geom.angle;
    txtAngle.ValueChangedFcn = @(txt, event) set_angle(W, txt.Value);
    
    % FEM nodes
    % Create labels and input fields for wingspan, chord, angle, and case
    lblchordnodes = uilabel(fig, 'Text', 'Chord nodes:', 'Position', [300, 500, 100, 22]);
    
    % Create a numeric input field for the wingspan
    txtchordnodes = uieditfield(fig, 'numeric', 'Position', [400, 500, 100, 22]);
    
    % Set the initial value of the input field
    txtchordnodes.Value = W.model.chorddiv;

    txtchordnodes.ValueChangedFcn = @(txt, event) set_chordnodes(W, txt.Value);

    lblspannodes = uilabel(fig, 'Text', 'Span nodes:', 'Position', [300, 450, 100, 22]);
    txtspannodes = uieditfield(fig, 'numeric', 'Position', [400, 450, 100, 22]);

    txtspannodes.Value = W.model.spandiv;
    txtspannodes.ValueChangedFcn = @(txt, event) set_spannodes(W, txt.Value);

    lblspanpanels = uilabel(fig, 'Text', 'Span Panels:', 'Position', [300, 400, 100, 22]);
    txtspanpanels = uieditfield(fig, 'numeric', 'Position', [400, 400, 100, 22]);

    txtspanpanels.Value = W.model.panspan;
    txtspanpanels.ValueChangedFcn = @(txt, event) set_spanpanels(W, txt.Value);

    lblchordpanels = uilabel(fig, 'Text', 'Span Panels:', 'Position', [300, 350, 100, 22]);
    txtchordpanels = uieditfield(fig, 'numeric', 'Position', [400, 350, 100, 22]);

    txtchordpanels.Value = W.model.panchord;
    txtchordpanels.ValueChangedFcn = @(txt, event) set_chordpanels(W, txt.Value);

    lblCase = uilabel(fig, 'Text', 'Case:', 'Position', [550, 350, 100, 22]);
    dropdownCase = uidropdown(fig, 'Items', {'flutter', 'divergence'}, 'Position', [650, 350, 100, 22]);

    dropdownCase.ValueChangedFcn = @(dd, event) set_case(W, dd.Value);

    % Create a "Calculate" button
    btnCalculate = uibutton(fig, 'Text', 'Calculate', 'Position', [50, 300, 100, 22], ...
        'ButtonPushedFcn', @(btnCalculate, event) calculateButtonPushed(txtWingspan.Value, txtChord.Value, txtAngle.Value, ax));

    % Create a "Close" button
    btnClose = uibutton(fig, 'Text', 'Close', 'Position', [200, 300, 100, 22], ...
        'ButtonPushedFcn', @(btnClose, event) closeButtonPushed(fig));
    
    
    % Create "Show Nodes" button
    btnShowNodes = uibutton(fig, 'Text', 'Show Nodes', 'Position', [50, 250, 100, 22], ...
        'ButtonPushedFcn', @(btnShowNodes, event) showNodesButtonPushed(W));

    % Create "Show Panels" button
    btnShowPanels = uibutton(fig, 'Text', 'Show Panels', 'Position', [200, 250, 100, 22], ...
        'ButtonPushedFcn', @(btnShowPanels, event) showPanelsButtonPushed(W));
end

function showNodesButtonPushed(W)
    W.mesh_plot;

end

function set_span(W, txt)
     W.geom.span = txt;
     W.gen_fem;

end

function set_spanpanels(W,txt)
    W.model.panspan = txt;
    W.lattice;
    
end

function set_chordpanels(W, txt)
    W.model.panchord = txt;
    W.lattice;
end

function set_chordnodes(W, txt)
    W.model.chorddiv = txt;
    W.gen_fem;

end

function set_spannodes(W, txt)
    W.model.spandiv = txt;
    W.gen_fem;
end

function set_chord(W, txt)
    W.geom.chord = txt;
    W.gen_fem;
end

function set_angle(W, txt)
    W.geom.angle = txt;
    W.gen_fem;
end

function set_case(W, selectedValue)
    switch selectedValue
        case 'flutter'
            W.model.analysis = 1; % Update this to the appropriate parameter in your class
        case 'divergence'
            W.model.analysis = 2; % Update this to the appropriate parameter in your class
    end
end