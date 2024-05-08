simpleGUI()

function simpleGUI()
    % Create the main figure window
    fig = uifigure('Name', 'Simple GUI Example', 'Position', [100, 100, 800, 600]);

    % Create labels and input fields for wingspan, chord, and angle
    lblWingspan = uilabel(fig, 'Text', 'Wingspan:', 'Position', [50, 500, 100, 22]);
    txtWingspan = uieditfield(fig, 'numeric', 'Position', [150, 500, 100, 22]);

    lblChord = uilabel(fig, 'Text', 'Chord:', 'Position', [50, 450, 100, 22]);
    txtChord = uieditfield(fig, 'numeric', 'Position', [150, 450, 100, 22]);

    lblAngle = uilabel(fig, 'Text', 'Angle:', 'Position', [50, 400, 100, 22]);
    txtAngle = uieditfield(fig, 'numeric', 'Position', [150, 400, 100, 22]);

    % Create a "Calculate" button
    btnCalculate = uibutton(fig, 'Text', 'Calculate', 'Position', [50, 350, 100, 22], ...
        'ButtonPushedFcn', @(btnCalculate, event) calculateButtonPushed(txtWingspan.Value, txtChord.Value, txtAngle.Value));

    % Create a "Close" button
    btnClose = uibutton(fig, 'Text', 'Close', 'Position', [200, 350, 100, 22], ...
        'ButtonPushedFcn', @(btnClose, event) closeButtonPushed(fig));
end

function calculateButtonPushed(wingspan, chord, angle)
    % Call the AERONastran function with the provided inputs
    AERONastran(wingspan, chord, angle);
end

function closeButtonPushed(fig)
    % Close the figure window
    close(fig);
end

