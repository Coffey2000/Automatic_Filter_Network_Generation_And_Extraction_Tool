global limit_noise_suppressor N M_matrix cross_connection_matrix remaining_cross_connection_matrix B_enable TOTAL_NUM_EXTRACTION working_node ending_node
N = 4;
%epsilon = 0.8665;
epsilon = 1.1548;
syms s w

ES = [1, 2.4015-1i*0.7591, 3.6706-1i*2.1950, 2.4874-1i*3.6255, -0.1268-1i*2.0658];
FS = [1.0, -1i*0.7591, 0.7869, -1i*0.5432, 0.0208];
PS = [-1, 1i*3.1299, 2.3899];
% ES = [1, 2.4236, 4.0098, 3.8108, 2.6022];
% FS = [1, 0, 1.0729, 0, 0.1641];
% PS = [-1, 0, -2.25];

Plot_S_from_extracted_M_matrix = true;      % Enable to plot S11 and S21 from the extracted M matrix
Bandwidth = 40e6;                           % Bandwidth
center_freq = 3.6e9;                        % Center frequency
freq_start = 3.5e9;                         % Start frequency (Hz)
freq_end = 3.7e9;                           % End frequency (Hz)

steps = 1000;                               % Number of steps

Network_Extraction_Force_Ending_With_Cross_Coupling = 1;

M_matrix = zeros(N, N);
limit_noise_suppressor = 1e-3;
TOTAL_NUM_EXTRACTION = 0;
working_node = 0;
ending_node = N + 1;

valid = 0;

B_enable = readmatrix("B_enable.csv");
cross_connection_matrix = readmatrix("cross_connection_matrix.csv");
remaining_cross_connection_matrix = cross_connection_matrix;
if all(~cross_connection_matrix, 'all')
    msgbox("Connection matrix is empty.", "Warning", "warn");
    valid = 0;
elseif N~=size(B_enable, 2)
    msgbox("Size of imported FIR components doen't match with the filter order.", "Warning", "warn");
    valid = 0;
elseif N ~=size(cross_connection_matrix, 2)
    msgbox("Size of imported connection matrix doen't match with the filter order.", "Warning", "warn");
    valid = 0;
else
    num_of_elements = ceil(N/2) - 1;
    for element = 1:1:num_of_elements
        if cross_connection_matrix(element, N - element) == 1 && cross_connection_matrix(element+1, N - element + 1) == 1
            msgbox('Crossing diagonal couplings cannot be enabled at the same time.','Warning', 'warn');
            valid = 0;
            break
        else
            valid = 1;
        end
    end
end


if valid
    WB2 = waitbar(0,'Extracting network components ....');

    num_coupling = sum(cross_connection_matrix, "all")/2;
    num_FIR = sum(B_enable, "all");
    TOTAL_NUM_EXTRACTION = 2*N + 1 + num_coupling + num_FIR;

    A = sym('A',[2 TOTAL_NUM_EXTRACTION + 1]);
    B = sym('B',[2 TOTAL_NUM_EXTRACTION + 1]);
    C = sym('C',[2 TOTAL_NUM_EXTRACTION + 1]);
    D = sym('D',[2 TOTAL_NUM_EXTRACTION + 1]);
    P = sym('P',[2 TOTAL_NUM_EXTRACTION + 1]);
    Extracted_C = zeros(1, N);
    Extracted_B = zeros(1, N);
    

    [A, B, C, D, P] = EF2ABCD(A, B, C, D, P, ES, FS, PS, epsilon);
    working_node = 0;
    ending_node = N + 1;
    [A, B, C, D, P] = series_unit_INV_extraction(A, B, C, D, P);
    working_node = next_working_node(working_node, Extracted_C);
    [A, B, C, D, P] = reverse(A, B, C, D, P);
    [working_node, ending_node] = swap(working_node, ending_node);
    [A, B, C, D, P] = series_unit_INV_extraction(A, B, C, D, P);
    working_node = next_working_node(working_node, Extracted_C);
    [A, B, C, D, P] = reverse(A, B, C, D, P);
    [working_node, ending_node] = swap(working_node, ending_node);

    for CURRENT_NUM_EXTRACTION = 3:1:TOTAL_NUM_EXTRACTION
        waitbar(CURRENT_NUM_EXTRACTION/TOTAL_NUM_EXTRACTION, WB2,'Extracting network components ....');
        if CURRENT_NUM_EXTRACTION == TOTAL_NUM_EXTRACTION
            [A, B, C, D, P] = parallel_INV_extraction(A, B, C, D, P);
        elseif Extracted_C(working_node) == 0
            [Extracted_C, A, B, C, D, P] = C_extraction(Extracted_C, A, B, C, D, P);
        elseif  B_enable(working_node) == 1 && Extracted_B(working_node) == 0
            [Extracted_B, A, B, C, D, P] = B_extraction(Extracted_B, A, B, C, D, P);
        else
            [A, B, C, D, P] = reverse(A, B, C, D, P);
            [working_node, ending_node] = swap(working_node, ending_node);
            if Extracted_C(working_node) == 0
                [Extracted_C, A, B, C, D, P] = C_extraction(Extracted_C, A, B, C, D, P);
            elseif  B_enable(working_node) == 1 && Extracted_B(working_node) == 0
                [Extracted_B, A, B, C, D, P] = B_extraction(Extracted_B, A, B, C, D, P);
            elseif cross_connection_matrix(working_node, ending_node) == 1 && M_matrix(working_node, ending_node) == 0
                [A, B, C, D, P] = parallel_INV_extraction(A, B, C, D, P);
            else
                working_connection = sum(remaining_cross_connection_matrix(working_node, :), "all");
                ending_connection = sum(remaining_cross_connection_matrix(ending_node, :), "all");
                if working_connection>ending_connection
                    [A, B, C, D, P] = reverse(A, B, C, D, P);
                    [working_node, ending_node] = swap(working_node, ending_node);
                    [A, B, C, D, P] = series_unit_INV_extraction(A, B, C, D, P);
                    working_node = next_working_node(working_node, Extracted_C);
                elseif working_connection<ending_connection
                    [A, B, C, D, P] = series_unit_INV_extraction(A, B, C, D, P);
                    working_node = next_working_node(working_node, Extracted_C);
                else
                    if (abs(next_working_node(working_node, Extracted_C) - (N+1-ending_node)) >= 2) || (Network_Extraction_Force_Ending_With_Cross_Coupling && ((working_node>ceil(N/2) && next_working_node(working_node, Extracted_C)<=ceil(N/2)) || (working_node<=ceil(N/2) && next_working_node(working_node, Extracted_C)>ceil(N/2)))) 
                        [A, B, C, D, P] = reverse(A, B, C, D, P);
                        [working_node, ending_node] = swap(working_node, ending_node);
                        [A, B, C, D, P] = series_unit_INV_extraction(A, B, C, D, P);
                        working_node = next_working_node(working_node, Extracted_C);
                    else
                        [A, B, C, D, P] = series_unit_INV_extraction(A, B, C, D, P);
                        working_node = next_working_node(working_node, Extracted_C);
                    end
                end
            end
        end
    end
    close(WB2)

    
    for i = 1:1:N
        if Extracted_C(i) == 0
            failed = 1;
            break;
        end
        failed = 0;
    end

    if failed
        msgbox("Network extraction failed! Please try another coupling connection.", "Warning", "warn");
    end

    scaled_M_matrix = zeros(N,N);
    scaled_Extracted_C = ones(1, N);
    scaled_Extracted_B = zeros(1, N);


    for i = 1:1:N
        for j = 1:1:N
            if M_matrix(i,j) ~= 0
                scaled_M_matrix(i,j) = M_matrix(i,j)/sqrt(Extracted_C(i)*Extracted_C(j));
            end
        end
    end

    for i = 1:1:N
        scaled_Extracted_B(i) = Extracted_B(i)/Extracted_C(i);
    end

    SL_scaled_M_matrix = zeros(N+2, N+2);
    SL_scaled_M_matrix(2:end-1, 2:end-1) = scaled_M_matrix;
    SL_scaled_M_matrix(1, 2) = 1/sqrt(Extracted_C(1));
    SL_scaled_M_matrix(2, 1) = 1/sqrt(Extracted_C(1));
    SL_scaled_M_matrix(end-1, end) = 1/sqrt(Extracted_C(end));
    SL_scaled_M_matrix(end, end-1) = 1/sqrt(Extracted_C(end));


    RS = SL_scaled_M_matrix(1, 2)^2;
    RL = SL_scaled_M_matrix(end-1, end)^2;
    

%% Plot S Parameters from Network Extraction

    if Plot_S_from_extracted_M_matrix && ~failed
        step_size = (freq_end-freq_start)/steps;
        
        R = zeros(N, N);
        R(1,1) = RS;
        R(end, end) = RL;
        
        S11_M_matrix = zeros(1, steps + 1);
        S21_M_matrix = zeros(1, steps + 1);
        
        WB3 = waitbar(0,'Calculating S11 and S21 from extracted M Matrix ....');
        
        for f = freq_start : step_size : freq_end
        waitbar((f-freq_start)/(freq_end-freq_start), WB3,'Calculating S11 and S21 from extracted M Matrix ....');
    
        lambda = center_freq/Bandwidth*(f/center_freq-center_freq/f);
        
        A_matrix = lambda*eye(N) - 1i*R + scaled_M_matrix;
        A_matrix_inv = A_matrix^(-1);
        
        S11_M_matrix((f - freq_start)/step_size + 1) = 1 + 2*1i*RS*A_matrix_inv(1,1);
        S21_M_matrix(round((f - freq_start)/step_size + 1)) = -2*1i*sqrt(RS*RL)*A_matrix_inv(N,1);
        end
    
        close(WB3)
        
        freq = linspace(freq_start, freq_end, steps + 1);
        
        figure;
        ref = plot(freq, 20*log10(abs(S11_M_matrix)));
        hold on
        trans = plot(freq, 20*log10(abs(S21_M_matrix)));
        hold off
        
        legend([ref, trans], "S11", "S21")
        xlabel("Frequency (Hz)")
        ylabel("dB")
        title("Extracted M Matrix S11 S21 vs Frequency")
    end


else
end




%% Functions


function array_done = pad2N_inf(array, N)
    array_done = zeros(1,N+1);
    length = size(array,2);
    array_done(N-length+2:N+1) = array;
    array_done(1:N-length+1) = inf;
end


function array_done = pad2N(array)
    global N
    array_done = zeros(1,N+1);
    length = size(array,2);
    if length>=N+1
        array_done = array(end - N:end);
    else
        array_done(N-length+2:N+1) = array;
        array_done(1:N-length+1) = 0;
    end
end


function [A, B, C, D, P] = EF2ABCD(Ain, Bin, Cin, Din, Pin, ES, FS, PS, epsilon)
    global N
    syms s

    ES = pad2N(ES);
    FS = pad2N(FS);
    PS = pad2N(PS);
    
    
    AS = 1i*imag(ES + FS);
    for i = 2:2:N+1
        AS(i) = real(ES(i) + FS(i));
    end
    
    BS = 1i*imag(ES + FS);
    for i = 1:2:N+1
        BS(i) = real(ES(i) + FS(i));
    end
    
    CS = 1i*imag(ES - FS);
    for i = 1:2:N+1
        CS(i) = real(ES(i) - FS(i));
    end
    
    DS = 1i*imag(ES - FS);
    for i = 2:2:N+1
        DS(i) = real(ES(i) - FS(i));
    end
    
    Prem = PS/epsilon;

    A_current_sym = poly2sym(AS, s);
    B_current_sym = poly2sym(BS, s);
    C_current_sym = poly2sym(CS, s);
    D_current_sym = poly2sym(DS, s);
    P_current_sym = poly2sym(Prem, s);


    [A, B, C, D, P] = store_ABCD(Ain, Bin, Cin, Din, Pin, A_current_sym, B_current_sym, C_current_sym, D_current_sym, P_current_sym, 0);

end

function [A, B, C, D, P] = series_unit_INV_extraction(Ain, Bin, Cin, Din, Pin)
    global working_node M_matrix N
    [~, index_last_element, row] = find_index_and_row(Ain);

    A_current_sym = -1i*Cin(row, index_last_element);
    B_current_sym = -1i*Din(row, index_last_element);
    C_current_sym = -1i*Ain(row, index_last_element);
    D_current_sym = -1i*Bin(row, index_last_element);
    P_current_sym = Pin(row, index_last_element);

%     sym2poly(A_current_sym)
%     sym2poly(B_current_sym)
%     sym2poly(C_current_sym)
%     sym2poly(D_current_sym)

    [A, B, C, D, P] = store_ABCD(Ain, Bin, Cin, Din, Pin, A_current_sym, B_current_sym, C_current_sym, D_current_sym, P_current_sym, 0);

    if working_node == 0 || working_node == N+1
    1 + 1;
    elseif working_node<ceil(N/2)
        if M_matrix(working_node, working_node + 1) == 0
            M_matrix(working_node, working_node + 1) = 1;
            M_matrix(working_node + 1, working_node) = 1;
        else
            M_matrix(working_node, working_node - 1) = 1;
            M_matrix(working_node - 1, working_node) = 1;
        end
    else
        if M_matrix(working_node, working_node - 1) == 0
            M_matrix(working_node, working_node - 1) = 1;
            M_matrix(working_node - 1, working_node) = 1;
        else
            M_matrix(working_node, working_node + 1) = 1;
            M_matrix(working_node + 1, working_node) = 1;
        end
    end
end


function [Extracted_C, A, B, C, D, P] = C_extraction(Extracted_Cin, Ain, Bin, Cin, Din, Pin)
    global working_node
    [~, index_last_element, row] = find_index_and_row(Ain);

%     disp("input of C extraction")
%     sym2poly(Din(row, index_last_element))
%     sym2poly(Bin(row, index_last_element))

    syms s
    C_current = double(limit(noise_suppress(Din(row, index_last_element))/(s*noise_suppress(Bin(row, index_last_element))), s, inf));

%     sym2poly(noise_suppress(Din(row, index_last_element)))
%     sym2poly(noise_suppress(Bin(row, index_last_element)))
    
    A_current_sym = Ain(row, index_last_element);
    B_current_sym = Bin(row, index_last_element);
    C_current_sym = Cin(row, index_last_element) - s*C_current*Ain(row, index_last_element);
    D_current_sym = Din(row, index_last_element) - s*C_current*Bin(row, index_last_element);
    P_current_sym = Pin(row, index_last_element);

    [A, B, C, D, P] = store_ABCD(Ain, Bin, Cin, Din, Pin, A_current_sym, B_current_sym, C_current_sym, D_current_sym, P_current_sym, 0);
    Extracted_C = Extracted_Cin;
    Extracted_C(working_node) = C_current;
end



function [Extracted_B, A, B, C, D, P] = B_extraction(Extracted_Bin, Ain, Bin, Cin, Din, Pin)
    global working_node
    [~, index_last_element, row] = find_index_and_row(Ain);
    
    syms s
    B_current = imag(double(limit(noise_suppress(Din(row, index_last_element))/noise_suppress(Bin(row, index_last_element)), s, inf)));

%     sym2poly(noise_suppress(Din(row, index_last_element)))
%     sym2poly(noise_suppress(Bin(row, index_last_element)))

    A_current_sym = Ain(row, index_last_element);
    B_current_sym = Bin(row, index_last_element);
    C_current_sym = Cin(row, index_last_element) - 1i*B_current*Ain(row, index_last_element);
    D_current_sym = Din(row, index_last_element) - 1i*B_current*Bin(row, index_last_element);
    P_current_sym = Pin(row, index_last_element);

%     sym2poly(A_current_sym)
%     sym2poly(B_current_sym)
%     sym2poly(C_current_sym)
%     sym2poly(D_current_sym)

    [A, B, C, D, P] = store_ABCD(Ain, Bin, Cin, Din, Pin, A_current_sym, B_current_sym, C_current_sym, D_current_sym, P_current_sym, 0);
    Extracted_B = Extracted_Bin;
    Extracted_B(working_node) = B_current;
end

function [A, B, C, D, P] = reverse(Ain, Bin, Cin, Din, Pin)
    [~, index_last_element, row] = find_index_and_row(Ain);

%     sym2poly(Ain(row, index_last_element))
%     sym2poly(Bin(row, index_last_element))
%     sym2poly(Cin(row, index_last_element))
%     sym2poly(Din(row, index_last_element))

    A_current_sym = Din(row, index_last_element);
    B_current_sym = Bin(row, index_last_element);
    C_current_sym = Cin(row, index_last_element);
    D_current_sym = Ain(row, index_last_element);
    P_current_sym = Pin(row, index_last_element);

%     disp("D output of reverse")
%     sym2poly(D_current_sym)

    [A, B, C, D, P] = store_ABCD(Ain, Bin, Cin, Din, Pin, A_current_sym, B_current_sym, C_current_sym, D_current_sym, P_current_sym, 1);
end

function [A, B, C, D, P] = parallel_INV_extraction(Ain, Bin, Cin, Din, Pin)
    global working_node ending_node M_matrix remaining_cross_connection_matrix
    [~, index_last_element, row] = find_index_and_row(Ain);

    syms s
    M_parallel_current = double(limit(-1*noise_suppress(Pin(row, index_last_element))/noise_suppress(Bin(row, index_last_element)), s, inf));

    A_current_sym = Ain(row, index_last_element);
    B_current_sym = Bin(row, index_last_element);
    C_current_sym = Cin(row, index_last_element) + 2*M_parallel_current*Pin(row, index_last_element)+ M_parallel_current^2*Bin(row, index_last_element);
    D_current_sym = Din(row, index_last_element);
    P_current_sym = Pin(row, index_last_element) + M_parallel_current*Bin(row, index_last_element);

    [A, B, C, D, P] = store_ABCD(Ain, Bin, Cin, Din, Pin, A_current_sym, B_current_sym, C_current_sym, D_current_sym, P_current_sym, 0);
    M_matrix(working_node, ending_node) = M_parallel_current;
    M_matrix(ending_node, working_node) = M_parallel_current;
    remaining_cross_connection_matrix(working_node, ending_node) = 0;
    remaining_cross_connection_matrix(ending_node, working_node) = 0;
end


function [A, B, C, D, P] = store_ABCD(Ain, Bin, Cin, Din, Pin, A_current_sym, B_current_sym, C_current_sym, D_current_sym, P_current_sym, reverse_flag)
    [index_to_write, ~, row] = find_index_and_row(Ain);

    A = Ain;
    B = Bin;
    C = Cin;
    D = Din;
    P = Pin;
    
%     sym2poly(C_current_sym)

    deleted = 0;

    if reverse_flag
        index_to_write = index_to_write - 1;
        if row == 1
           row = 2;
           if is_A_Occupied(Ain, 2, index_to_write)
               [A, B, C, D, P] = delete_ABCD(A, B, C, D, P, 1, index_to_write);
               deleted = 1;
           end

        else
           row = 1;
           if is_A_Occupied(Ain, 1, index_to_write)
               [A, B, C, D, P] = delete_ABCD(A, B, C, D, P, 2, index_to_write);
               deleted = 1;
           end
        end
    end

    if ~deleted
        A(row, index_to_write) = A_current_sym;
        B(row, index_to_write) = B_current_sym;
        C(row, index_to_write) = C_current_sym;
        D(row, index_to_write) = D_current_sym;
        P(row, index_to_write) = P_current_sym;
    end
end



function [index_to_write, index_last_element, row] = find_index_and_row(Ain)
    global TOTAL_NUM_EXTRACTION
    index_to_write = 1;
    
    for i = 1:1:TOTAL_NUM_EXTRACTION + 1
        if ~is_A_Occupied(Ain, 1, i) && ~is_A_Occupied(Ain, 2, i)
            index_to_write = i;
            break
        end
    end

    num_sequential_full_col = 0;
    last_vacant_col = 1;
    if index_to_write ~= 1
        for j = 1:1:index_to_write - 1
            if is_A_Occupied(Ain, 1, index_to_write - j) && is_A_Occupied(Ain, 2, index_to_write - j)
                 num_sequential_full_col =  num_sequential_full_col + 1;
            else
                last_vacant_col = index_to_write - j;
                break
            end
        end
    end

    if is_A_Occupied(Ain, 1, last_vacant_col)
        last_vacant_row = 2;
    else
        last_vacant_row = 1;
    end


    if index_to_write == 1
        row = 1;
    elseif mod(num_sequential_full_col, 2) ~= 0
        row = last_vacant_row;
    else
        if last_vacant_row == 1
            row = 2;
        else
            row = 1;
        end
    end

    index_last_element = index_to_write - 1;
end



function noiseless_sym=noise_suppress(noisy_sym)
    global limit_noise_suppressor
    syms s
%     sym2poly(noisy_sym)
    noiseless = pad2N(sym2poly(noisy_sym));
    noiseless(abs(noiseless)<limit_noise_suppressor) = 0;
    noiseless_sym = poly2sym(noiseless, s);
end

function [b, a] = swap(a, b)
end

function working_node = next_working_node(working_node_current, Extracted_C)
    global N

    if working_node_current == 0
        working_node = 1;
    elseif working_node_current == N + 1
        working_node = N;
    elseif working_node_current>ceil(N/2)
        if Extracted_C(working_node_current - 1) == 0
            working_node = working_node_current - 1;
        else
            working_node = working_node_current + 1; 
        end
    else
        if Extracted_C(working_node_current + 1) == 0
            working_node = working_node_current + 1;
        else
            working_node = working_node_current - 1;
        end
    end
end

function occupied = is_A_Occupied(Ain, row, col)
    occupied = ~strcmp(string(Ain(row,col)),"A" + string(row) + "_" + string(col));
end

function [A, B, C, D, P] = delete_ABCD(Ain, Bin, Cin, Din, Pin, row, col)
    A = Ain;
    B = Bin;
    C = Cin;
    D = Din;
    P = Pin;
    
    A(row, col) = sym("A" + string(row) + "_" + string(col));
    B(row, col) = sym("B" + string(row) + "_" + string(col));
    C(row, col) = sym("C" + string(row) + "_" + string(col));
    D(row, col) = sym("D" + string(row) + "_" + string(col));
    P(row, col) = sym("P" + string(row) + "_" + string(col));
end