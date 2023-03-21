global cross_connection_matrix B_enable N

%% Filter Order Setup
N = 4;

%% GUI
cross_connection_matrix = zeros(N, N);
B_enable = zeros(1, N);


if isEven(N)
    window_width = (ceil(N/2) + 1)*230;
    num_of_elements = ceil(N/2);
    
    % Create a figure window:
    GUI = uifigure('Position',[500 500 window_width 400],'Name','Filter Connection Options');
    
    uitextarea('Parent', GUI, 'Position',[60 290 100 30], 'BackgroundColor','black', 'FontColor', 'white', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Value', "Source");
    uitextarea('Parent', GUI, 'Position',[60 110 100 30], 'BackgroundColor','black', 'FontColor', 'white', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Value', "Load");

    for element = 1:1:num_of_elements
        uitextarea('Parent', GUI, 'Position',[200+230*(element-1) 300 150 6], 'BackgroundColor','black' );
        uitextarea('Parent', GUI, 'Position',[370+230*(element-1) 290 30 30], 'BackgroundColor','black', 'FontColor', 'white', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Value', string(element));
    end
    
    for element = 1:1:num_of_elements
        uitextarea('Parent', GUI, 'Position',[200+230*(element-1) 120 150 6], 'BackgroundColor','black' );
        uitextarea('Parent', GUI, 'Position',[370+230*(element-1) 110 30 30], 'BackgroundColor','black', 'FontColor', 'white', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Value', string(N - element + 1));
    end
    
    uitextarea('Parent', GUI, 'Position',[window_width-77 150 6 130], 'BackgroundColor','black' );
    
    for element = 1:1:num_of_elements
    cbx = uicheckbox(GUI,'Position',[380+230*(element-1) 330 102 15], 'ValueChangedFcn',@(cbx,event) BChanged(cbx), 'Text', "jB"+string(element));
    cbx = uicheckbox(GUI,'Position',[380+230*(element-1) 85 102 15], 'ValueChangedFcn',@(cbx,event) BChanged(cbx), 'Text', "jB"+string(N - element + 1));
    end
    
    num_of_elements = ceil(N/2) - 1;
    for element = 1:1:num_of_elements
    cbx = uicheckbox(GUI,'Position',[380+230*(element-1) 205 102 15], 'ValueChangedFcn',@(cbx,event) CouplingChanged(cbx), 'Text', "M"+string(element)+string(N - element + 1));
    cbx = uicheckbox(GUI,'Position',[500+230*(element-1) 215 102 15], 'ValueChangedFcn',@(cbx,event) CouplingChanged(cbx), 'Text', "M"+string(element)+string(N - element));
    cbx = uicheckbox(GUI,'Position',[500+230*(element-1) 190 102 15], 'ValueChangedFcn',@(cbx,event) CouplingChanged(cbx), 'Text', "M"+string(element+1)+string(N - element + 1));
    end
else
    window_width = (ceil(N/2) + 1)*230 + 40;
    num_of_elements = ceil(N/2) - 1;
    
    % Create a figure window:
    GUI = uifigure('Position',[500 500 window_width 400],'Name','Filter Connection Options');
    
    uitextarea('Parent', GUI, 'Position',[60 290 100 30], 'BackgroundColor','black', 'FontColor', 'white', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Value', "Source");
    uitextarea('Parent', GUI, 'Position',[60 110 100 30], 'BackgroundColor','black', 'FontColor', 'white', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Value', "Load");

    for element = 1:1:num_of_elements
        uitextarea('Parent', GUI, 'Position',[200+230*(element-1) 300 150 6], 'BackgroundColor','black' );
        uitextarea('Parent', GUI, 'Position',[370+230*(element-1) 290 30 30], 'BackgroundColor','black', 'FontColor', 'white', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Value', string(element));
    end
    
    uitextarea('Parent', GUI, 'Position',[370+230*(ceil(N/2)-1) 197 30 30], 'BackgroundColor','black', 'FontColor', 'white', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Value', string(ceil(N/2)));

    for element = 1:1:num_of_elements
        uitextarea('Parent', GUI, 'Position',[200+230*(element-1) 120 150 6], 'BackgroundColor','black' );
        uitextarea('Parent', GUI, 'Position',[370+230*(element-1) 110 30 30], 'BackgroundColor','black', 'FontColor', 'white', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Value', string(N - element + 1));
    end
    
    uitextarea('Parent', GUI, 'Position',[200+230*(ceil(N/2)-1) 300 189 6], 'BackgroundColor','black' );
    uitextarea('Parent', GUI, 'Position',[window_width-117 250 6 50], 'BackgroundColor','black' );
    uitextarea('Parent', GUI, 'Position',[window_width-117 120 6 50], 'BackgroundColor','black' );
    uitextarea('Parent', GUI, 'Position',[200+230*(ceil(N/2)-1) 120 189 6], 'BackgroundColor','black' );

    for element = 1:1:num_of_elements
    cbx = uicheckbox(GUI,'Position',[380+230*(element-1) 330 102 15], 'ValueChangedFcn',@(cbx,event) BChanged(cbx), 'Text', "jB"+string(element));
    cbx = uicheckbox(GUI,'Position',[380+230*(element-1) 85 102 15], 'ValueChangedFcn',@(cbx,event) BChanged(cbx), 'Text', "jB"+string(N - element + 1));
    end
    
    cbx = uicheckbox(GUI,'Position',[412+230*(ceil(N/2)-1) 205 102 15], 'ValueChangedFcn',@(cbx,event) BChanged(cbx), 'Text', "jB"+string(ceil(N/2)));

    num_of_elements = ceil(N/2) - 1;
    for element = 1:1:num_of_elements
    cbx = uicheckbox(GUI,'Position',[380+230*(element-1) 205 102 15], 'ValueChangedFcn',@(cbx,event) CouplingChanged(cbx), 'Text', "M"+string(element)+string(N - element + 1));
    if element ~= num_of_elements
        cbx = uicheckbox(GUI,'Position',[500+230*(element-1) 215 102 15], 'ValueChangedFcn',@(cbx,event) CouplingChanged(cbx), 'Text', "M"+string(element)+string(N - element));
        cbx = uicheckbox(GUI,'Position',[500+230*(element-1) 190 102 15], 'ValueChangedFcn',@(cbx,event) CouplingChanged(cbx), 'Text', "M"+string(element+1)+string(N - element + 1));
    end
    end

end

bt = uibutton(GUI,'Position',[round(window_width)/2-50 20 100 30], "Text", "OK","ButtonPushedFcn", @(bt,event) ButtonPushed(bt));



function CouplingChanged(cbx)
    global cross_connection_matrix
    val = cbx.Value;
    T = cbx.Text;
    row = round(str2double(extract(T, 2)));
    col = round(str2double(extract(T, 3)));
    if val
        cross_connection_matrix(row, col) = 1;
        cross_connection_matrix(col, row) = 1;
    else
        cross_connection_matrix(row, col) = 0;
        cross_connection_matrix(col, row) = 0;
    end
end


function BChanged(cbx)
    global B_enable
    val = cbx.Value;
    T = cbx.Text;
    col = round(str2double(extract(T, 3)));
    if val
        B_enable(col) = 1;
    else
        B_enable(col) = 0;
    end
end

function ButtonPushed(bt)
    global cross_connection_matrix B_enable N
    num_of_elements = ceil(N/2) - 1;
    valid = 1;
    for element = 1:1:num_of_elements
        if cross_connection_matrix(element, N - element) == 1 && cross_connection_matrix(element+1, N - element + 1) == 1
            msgbox('Crossing diagonal couplings cannot be enabled at the same time.','Warning', 'warn');
            valid = 0;
            break
        end
    end
    if valid
        writematrix(cross_connection_matrix, "cross_connection_matrix.csv");
        writematrix(B_enable, "B_enable.csv");
        msgbox('Cross connection matrix and FIR components are exported.','Success');
        close(bt.Parent);
        clearvars;
    end
end

function boolean1 = isEven(N)
    boolean1 = mod(N, 2) == 0;
end