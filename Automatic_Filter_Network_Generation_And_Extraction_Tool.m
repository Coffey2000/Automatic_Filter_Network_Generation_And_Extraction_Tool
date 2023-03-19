global limit_noise_suppressor N M_matrix cross_connection_matrix remaining_cross_connection_matrix B_enable TOTAL_NUM_EXTRACTION working_node ending_node

%% Filter Setup
N = 3;              % Filter order
RL = 23;            % Filter Return Loss (dB)
TZ = [inf];    % Array of frequencies of transmittion zeros (rad/s)
                    % For transmission zeros at infinity, type "inf"

%% Simulation Setup
Polynomial_Solver = "recursive";                                    % Choose between recursive solver and numerical solver for polynomial generation
Enable_Network_Extraction = true;                                   % Enable Network Extraction and generate M matrix
Network_Extraction_Force_Ending_With_Cross_Coupling = true;         % Force the ending coupling to be extracted as cross coupling

%% Debugging Tool                                        
                                                         % The two rows in the ABCDP matrixes represents the direction of the expression. 
                                                         % Row 1 contains expressions for forward direction while row 2 constains expressions for reverse direction. 
                                                         % Each expression in the ABCDP matrixes represent the remaining ABCDP expressions after (column - 1) extractions 
                                                         % or each of the ABCDP expression is used for (column)th extraction.

Enable_ABCDP_simplification = true;                      % Enable simplification of the ABCD and P expressions for easier understanding
round_to_decimal_places = 3;                             % Number of decimal places to round to in simplification

%% Plotting Setup
Plot_S_from_polynomials = true;             % Enable to plot S11 and S21 from the polynomials
normalized_freq_start = -3;                 % Start frequency (rad/s)
normalized_freq_end = 3;                    % End frequency (rad/s)


Plot_S_from_extracted_M_matrix = true;      % Enable to plot S11 and S21 from the extracted M matrix
Bandwidth = 40e6;                           % Bandwidth
center_freq = 3.6e9;                        % Center frequency
freq_start = 3.5e9;                         % Start frequency (Hz)
freq_end = 3.7e9;                           % End frequency (Hz)


steps = 1000;                               % Number of steps

%% Solve for P(w) and P(S)
syms s w 

TZ = pad2N_inf(TZ, N-1);
finite_TZ = TZ(TZ ~= inf);
num_of_finite_TZ = size(finite_TZ, 2);

PW_roots = TZ;
PW = poly(PW_roots);
PW = pad2N(PW);
PW_sym = poly2sym(PW, w);
PS_sym = subs(PW_sym, w, -1i*s);

% if mod(N-num_of_finite_TZ, 2) == 0
%     PS_sym = PS_sym * 1;
% end

PS = sym2poly(PS_sym);
PS = pad2N(PS);

%% Solve for F(w) and F(S)
if Polynomial_Solver == "recursive"
    WB0 = waitbar(0,'Generating filter polynomials ....');
    U = sym('U',[1 N]);
    V = sym('V',[1 N]);
    
    for n = 1:1:N
        if n == 1
            U(n) = w - 1/TZ(n);
            V(n) = sqrt(w^2 - 1)*(1-1/(TZ(n)^2))^0.5;
        else
            U(n) = w * U(n-1) - U(n-1)/TZ(n) + sqrt(w^2 - 1) * sqrt(1 - 1/(TZ(n)^2)) * V(n-1);
            V(n) = w * V(n-1) - V(n-1)/TZ(n) + sqrt(w^2 - 1) * sqrt(1 - 1/(TZ(n)^2)) * U(n-1);
        end
    end
    
    FW_roots = roots(sym2poly(U(N)));
    FW = poly(FW_roots);
    FW(abs(FW)<1e-10) = 0;
    FW_sym = poly2sym(FW, w);
    FS_sym = subs(FW_sym, w, -1i*s);
    FS = sym2poly(FS_sym);
    FS = pad2N(FS);

elseif Polynomial_Solver == "numerical"
    WB0 = waitbar(0,'Generating filter polynomials ....');
    innerTerm = 0;
    for i = 1 : 1 : num_of_finite_TZ
        innerTerm = innerTerm + acosh((w-1/finite_TZ(i))/(1-w/finite_TZ(i)));
    end
    innerTerm = innerTerm + (N - num_of_finite_TZ)*acosh(w);
    
    CW = @(w) cosh(innerTerm);
    
    PW = poly(TZ);
    PW = pad2N(PW);
    PW_sym = @(w) poly2sym(PW, w);
    FW_sym = @(w) CW(w) * PW_sym(w);
    FW_sym = matlabFunction(FW_sym, 'vars', w);
    
    
    checkPoints = linspace(-1,1,N);
    FW_roots = zeros(1,N);
    for i=1:1:N
        FW_roots(i) = fzero(FW_sym, checkPoints(i) + 0.1);
    end
    FW_roots(abs(FW_roots)<1e-10) = 0;
    
    fplot(FW_sym, [-1,1])
    xlabel("Frequency (rad)")
    title("F(W)")
    
    PW_sym = poly2sym(PW, w);
    PS_sym = subs(PW_sym, w, -1i*s);
%     if mod(N-num_of_finite_TZ, 2) == 0
%         PS_sym = PS_sym * 1i;
%     end
    
    PS = sym2poly(PS_sym);
    PS = pad2N(PS);
    
    FW = poly(FW_roots);
    FW(abs(FW)<1e-10) = 0;
    FW_sym = poly2sym(FW, w);
    FS_sym = subs(FW_sym, w, -1i*s);
    FS = sym2poly(FS_sym);
    FS = pad2N(FS);
else
    warndlg('Please choose a valid solver.','Warning');
end

waitbar(0.5, WB0,'Generating filter polynomials ....');


%% Solve for epsilon, E(w) and E(S)
epsilon = abs(double(1/sqrt(10^(RL/10)-1)*subs(PW_sym, w, 1)/subs(FW_sym, w, 1)));
reference_epsilon = double(1/sqrt(10^(RL/10)-1));

%ES2 = conv(PS, (PS))/epsilon^2 + conv(FS, (FS));
% INTER1 = conv(PS, conj(PS))/epsilon^2;
% INTER2 = conv(FS, conj(FS));
ES2 = conv(PS, conj(PS))/epsilon^2 + conv(FS, conj(FS));
ES2_sym = poly2sym(ES2, s);
ES2_roots = roots(ES2);
ES_roots = ES2_roots(real(ES2_roots)<0);
%ES_roots = ES2_roots(1:2);

waitbar(0.75, WB0,'Generating filter polynomials ....');

reference_ES_roots = zeros(1, N);
if num_of_finite_TZ == 0
    for k = 1:1:N
        reference_ES_roots(k) = -sinh(1/N*asinh(1/reference_epsilon))*sin((pi/2)*(2*k-1)/N) + 1i * cosh(1/N*asinh(1/reference_epsilon))*cos((pi/2)*(2*k-1)/N);
    end
    ES = poly(reference_ES_roots);
else
    ES = poly(ES_roots);
end

%reference_ES2 = conv(ES, conj(ES));
% reference_ES2_roots = roots(reference_ES2);
ES_sym = poly2sym(ES, s);
EW_sym = subs(ES_sym, s, 1i*w);
EW = sym2poly(EW_sym);

waitbar(1, WB0,'Generating filter polynomials ....');
close(WB0)


%% Plot S parameters from the polynomials


if Plot_S_from_polynomials
    step_size = (normalized_freq_end-normalized_freq_start)/steps;
    
    S11_polynomial = zeros(1, steps + 1);
    S21_polynomial = zeros(1, steps + 1);
    
    
    WB1 = waitbar(0,'Calculating S11 and S21 from polynomials ....');
    
    for w_current = normalized_freq_start : step_size : normalized_freq_end
        waitbar((w_current-normalized_freq_start)/(normalized_freq_end - normalized_freq_start), WB1,'Calculating S11 and S21 from polynomials ....');
        S21_polynomial(round((w_current - normalized_freq_start)/step_size + 1)) = double(abs(subs(PW_sym, w, w_current)))/(epsilon * double(abs(subs(EW_sym, w, w_current))));
        S11_polynomial(round((w_current - normalized_freq_start)/step_size + 1)) = double(abs(subs(FW_sym, w, w_current)))/double(abs(subs(EW_sym, w, w_current)));
    end
    
    close(WB1)
    
    freq = linspace(normalized_freq_start, normalized_freq_end, steps + 1);
    %freq_shifted = center_freq/Bandwidth*(freq./center_freq - center_freq./freq);
    
    figure;
    ref = plot(freq, 20*log10(abs(S11_polynomial)));
    hold on
    trans = plot(freq, 20*log10(abs(S21_polynomial)));
    hold off
    
    legend([ref, trans], "S11", "S21")
    xlabel("Normalized Frequency (rad)")
    ylabel("dB")
    title("Polynomials S11 S21 vs Frequency")
end


%% Network Extraction


if Enable_Network_Extraction
    M_matrix = zeros(N, N);
    limit_noise_suppressor = 1e-3;
    TOTAL_NUM_EXTRACTION = 0;
    working_node = 0;
    ending_node = N + 1;
    
    valid = 0;
    
    B_enable = readmatrix("B_enable.csv");
    cross_connection_matrix = readmatrix("cross_connection_matrix.csv");
    remaining_cross_connection_matrix = cross_connection_matrix;
    if N~=size(B_enable, 2)
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
                if isEven(N)
                    [A, B, C, D, P] = parallel_INV_extraction(A, B, C, D, P);
                else
                    [A, B, C, D, P] = series_unit_INV_extraction(A, B, C, D, P);
                end
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
    

        if Enable_ABCDP_simplification
            for i = 1:1:2
                for j = 1:1:TOTAL_NUM_EXTRACTION + 1
                    A_simplified(i, j) = vpa(poly2sym(round(sym2poly(A(i,j)), round_to_decimal_places), s));
                    B_simplified(i, j) = vpa(poly2sym(round(sym2poly(B(i,j)), round_to_decimal_places), s));
                    C_simplified(i, j) = vpa(poly2sym(round(sym2poly(C(i,j)), round_to_decimal_places), s));
                    D_simplified(i, j) = vpa(poly2sym(round(sym2poly(D(i,j)), round_to_decimal_places), s));
                    P_simplified(i, j) = vpa(poly2sym(round(sym2poly(P(i,j)), round_to_decimal_places), s));
                end
            end
        end
        
        

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

function boolean1 = isEven(N)
    boolean1 = mod(N, 2) == 0;
end