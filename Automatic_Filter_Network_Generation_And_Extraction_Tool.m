global limit_noise_suppressor N M_matrix cross_connection_matrix remaining_cross_connection_matrix B_enable TOTAL_NUM_EXTRACTION working_node ending_node

%% Filter Setup
N = 5;                          % Filter order
IL = 0;                         % Filter ripple level (dB) Set either IL or RL and set the other value to 0!
RL = 25;                        % Filter Return Loss (dB)
TZ = [inf];                % Array of frequencies of transmittion zeros (rad/s)
                                % For transmission zeros at infinity, type "inf"
Q = [inf inf inf inf inf];          % Unloaded quality factor of the resonators

%% Simulation Setup
Polynomial_Solver = "recursive";                                    % Choose between recursive solver and numerical solver for polynomial generation
Enable_Network_Extraction = true;                                   % Enable Network Extraction and generate M matrix
Network_Extraction_Force_Ending_With_Cross_Coupling = false;        % Force the ending coupling to be extracted as cross coupling

Find_Extraction_Solution_Ultra_Fast = true;                        % Enable to automatically find coupling connection solutions (skip mirror cross coupling connections)
Find_Extraction_Solution_Fast = false;                              % Enable to automatically find coupling connection solutions (try with FIRs that are all off or on)
Find_Extraction_Solution_All = false;                               % Enable to automatically find coupling connection solutions (try all combinations of FIRs)
target_num_solution = 1;                                            % Targeting number of couping connection solutions to find
limit_noise_suppressor = 1e-3;                                      % When doing limit, coefficients smaller than this number are cleared.

%% Debugging Tool                                        
                                                         % The two rows in the ABCDP matrixes represents the direction of the expression. 
                                                         % Row 1 contains expressions for forward direction while row 2 constains expressions for reverse direction. 
                                                         % Each expression in the ABCDP matrixes represent the remaining ABCDP expressions after (column - 1) extractions 
                                                         % or each of the ABCDP expression is used for (column)th extraction.

Enable_ABCDP_simplification = true;                      % Enable simplification of the ABCD and P expressions for easier understanding
Enable_UV_simplification = true;                         % Enable simplification of the U and V expressions for easier understanding
round_to_decimal_places = 3;                             % Number of decimal places to round to in simplification

%% Plotting Setup
Plot_S_from_polynomials = true;             % Enable to plot S11 and S21 from the polynomials
normalized_freq_start = -3;                 % Start frequency (rad/s)
normalized_freq_end = 3;                    % End frequency (rad/s)


Plot_S_from_extracted_M_matrix = true;      % Enable to plot S11 and S21 from the extracted M matrix
Bandwidth = 38.4e6;                          % Bandwidth
center_freq = 1.6e9;                          % Center frequency
freq_start = 1.54e9;                         % Start frequency (Hz)
freq_end = 1.66e9;                           % End frequency (Hz)


steps = 1000;                               % Number of steps

%% Prepare IL, RL and Q
if IL ~= 0 && RL ~= 0
    msgbox("IL and RL cannot both be set! Set either IL or RL and set the other value to 0. Using the specified RL in this simulation ....", "Warning", "warn");
elseif IL == 0 && RL == 0
    msgbox("IL and RL cannot both be 0! Set either IL or RL and set the other value to 0. Using the specified RL in this simulation ....", "Warning", "warn");
elseif IL ~= 0
    RL = -10*log10(1-10^(-IL/10));
elseif N ~= size(Q,2)
    msgbox("Q does not have the correct size! Adjusting Q for this simulation ....", "Warning", "warn");
    Q = pad2N_for_Q(Q);
end

Return_Loss = RL;


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
        
        if Enable_UV_simplification
            try
                U_simplified(n) = vpa(poly2sym(round(sym2poly(U(n)), round_to_decimal_places), w));
                V_simplified(n) = vpa(poly2sym(round(sym2poly(V(n)), round_to_decimal_places), w));
            catch
            end
        end
    end


    FW_roots = roots(sym2poly(U(N)));
    FW = poly(FW_roots);
    FW(abs(FW)<1e-10) = 0;
    FW_sym = poly2sym(FW, w);
    FS_sym = subs(FW_sym, w, -1i*s);
    FS = sym2poly(FS_sym);
    
    if FS(1) ~= 1
        FS = FS/FS(1);
        FS_sym = poly2sym(FS, s);
    end

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

    if FS(1) ~= 1
        FS = FS/FS(1);
        FS_sym = poly2sym(FS, s);
    end

    FS = pad2N(FS);
else
    warndlg('Please choose a valid solver.','Warning');
end

waitbar(0.5, WB0,'Generating filter polynomials ....');


%% Solve for epsilon, E(w) and E(S)
epsilon = abs(double(1/sqrt(10^(RL/10)-1)*subs(PW_sym, w, 1)/subs(FW_sym, w, 1)));
%pure_Chebyshev_epsilon = double(1/sqrt(10^(RL/10)-1));

if isEven(N)
    ES2 = conv(PS, PS)/epsilon^2 + conv(FS, FS);
else
    ES2 = -1 * conv(PS, PS)/epsilon^2 + conv(FS, FS);
end

% INTER1 = conv(PS, conj(PS))/epsilon^2;
% INTER2 = conv(FS, conj(FS));
%ES2 = conv(PS, conj(PS))/epsilon^2 + conv(FS, conj(FS));
ES2_sym = poly2sym(ES2, s);
ES2_roots = roots(ES2);
ES_roots = ES2_roots(real(ES2_roots)<-1e-8);
%ES_roots = ES2_roots(1:2);

waitbar(0.75, WB0,'Generating filter polynomials ....');

% pure_Chebyshev_ES_roots = zeros(1, N);
% if num_of_finite_TZ == 0
%     for k = 1:1:N
%         pure_Chebyshev_ES_roots(k) = -sinh(1/N*asinh(1/pure_Chebyshev_epsilon))*sin((pi/2)*(2*k-1)/N) + 1i * cosh(1/N*asinh(1/pure_Chebyshev_epsilon))*cos((pi/2)*(2*k-1)/N);
%     end
%     ES = poly(pure_Chebyshev_ES_roots);
% else
%     ES = poly(ES_roots);
% end

% for k = 1:1:N
%     pure_Chebyshev_ES_roots(k) = -sinh(1/N*asinh(1/pure_Chebyshev_epsilon))*sin((pi/2)*(2*k-1)/N) + 1i * cosh(1/N*asinh(1/pure_Chebyshev_epsilon))*cos((pi/2)*(2*k-1)/N);
% end
% ES = poly(pure_Chebyshev_ES_roots);

ES = poly(ES_roots);
%reference_ES2 = conv(ES, conj(ES));
% reference_ES2_roots = roots(reference_ES2);
ES_sym = poly2sym(ES, s);
EW_sym = subs(ES_sym, s, 1i*w);
EW = sym2poly(EW_sym);
ES = pad2N(ES);

waitbar(1, WB0,'Generating filter polynomials ....');
close(WB0)


%% Calculate and Plot S parameters from the polynomials

step_size = (normalized_freq_end-normalized_freq_start)/steps;

S11_polynomial = zeros(1, steps + 1);
S21_polynomial = zeros(1, steps + 1);


WB1 = waitbar(0,'Calculating S11 and S21 from polynomials ....');

for normalized_freq_current = normalized_freq_start : step_size : normalized_freq_end
    waitbar((normalized_freq_current-normalized_freq_start)/(normalized_freq_end - normalized_freq_start), WB1,'Calculating S11 and S21 from polynomials ....');
    S21_polynomial(round((normalized_freq_current - normalized_freq_start)/step_size + 1)) = abs(subs(PW_sym, w, normalized_freq_current)/(epsilon * subs(EW_sym, w, normalized_freq_current)));
    S11_polynomial(round((normalized_freq_current - normalized_freq_start)/step_size + 1)) = abs(subs(FW_sym, w, normalized_freq_current)/subs(EW_sym, w, normalized_freq_current));
end

close(WB1)
    
if Plot_S_from_polynomials
    freq = linspace(normalized_freq_start, normalized_freq_end, steps + 1);
    
    figure;
    ref = plot(freq, 20*log10(S11_polynomial));
    hold on
    trans = plot(freq, 20*log10(S21_polynomial));
    hold off
    
    legend([ref, trans], "S11", "S21")
    xlabel("Normalized Frequency (rad)")
    ylabel("dB")
    title("Polynomials S11 S21 vs Frequency")
end


%% Network Extraction


if Enable_Network_Extraction
    M_matrix = zeros(N, N);
    TOTAL_NUM_EXTRACTION = 0;
    working_node = 0;
    ending_node = N + 1;
    
    valid = 1;
    
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
            end
        end
    end
    
    
    if valid
        try
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
                    if abs(M_matrix(working_node, ending_node) + 1) < 1e-9
                        M_matrix(working_node, ending_node) = 1;
                        M_matrix(ending_node, working_node) = 1;
                    end
                elseif Extracted_C(working_node) == 0
                    [Extracted_C, A, B, C, D, P] = C_extraction(Extracted_C, A, B, C, D, P);
                elseif B_enable(working_node) == 1 && Extracted_B(working_node) == 0
                    [Extracted_B, A, B, C, D, P] = B_extraction(Extracted_B, A, B, C, D, P);
                elseif Extracted_C(ending_node) == 0
                    [A, B, C, D, P] = reverse(A, B, C, D, P);
                    [working_node, ending_node] = swap(working_node, ending_node);
                    [Extracted_C, A, B, C, D, P] = C_extraction(Extracted_C, A, B, C, D, P);
                elseif B_enable(ending_node) == 1 && Extracted_B(ending_node) == 0
                    [A, B, C, D, P] = reverse(A, B, C, D, P);
                    [working_node, ending_node] = swap(working_node, ending_node);
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

            close(WB2)


            if Enable_ABCDP_simplification
                    A_simplified = sym('A',[2 TOTAL_NUM_EXTRACTION + 1]);
                    B_simplified = sym('B',[2 TOTAL_NUM_EXTRACTION + 1]);
                    C_simplified = sym('C',[2 TOTAL_NUM_EXTRACTION + 1]);
                    D_simplified = sym('D',[2 TOTAL_NUM_EXTRACTION + 1]);
                    P_simplified = sym('P',[2 TOTAL_NUM_EXTRACTION + 1]);
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
            
            
            current_trial_successful = 0;
    
            if is_C_full(Extracted_C)
                if is_ABCD_finished(A, B, C, D)
                    if is_response_identical(Extracted_B, Extracted_C, normalized_freq_start, normalized_freq_end, steps, S11_polynomial, S21_polynomial)
                        current_trial_successful = 1;
                    end
                end
            end

    
            if current_trial_successful
                scaled_M_matrix = zeros(N,N);
                scaled_Extracted_C = ones(1, N);
                scaled_Extracted_B = zeros(1, N);
            
                
                for i = 1:1:N
                    scaled_Extracted_B(i) = Extracted_B(i)/Extracted_C(i);
                end
               
        
                for i = 1:1:N
                    for j = 1:1:N
                        if j == i
                            scaled_M_matrix(i,j) = scaled_Extracted_B(j);
                        elseif M_matrix(i,j) ~= 0
                            scaled_M_matrix(i,j) = M_matrix(i,j)/sqrt(Extracted_C(i)*Extracted_C(j));
                        end
                    end
                end
        
            
                SL_scaled_M_matrix = zeros(N+2, N+2);
                SL_scaled_M_matrix(2:end-1, 2:end-1) = scaled_M_matrix;
                SL_scaled_M_matrix(1, 2) = 1/sqrt(Extracted_C(1));
                SL_scaled_M_matrix(2, 1) = 1/sqrt(Extracted_C(1));
                SL_scaled_M_matrix(end-1, end) = 1/sqrt(Extracted_C(end));
                SL_scaled_M_matrix(end, end-1) = 1/sqrt(Extracted_C(end));
            
            
                RS = SL_scaled_M_matrix(1, 2)^2;
                RL = SL_scaled_M_matrix(end-1, end)^2;
            else
                msgbox("Network extraction failed! Please try another coupling connection.", "Warning", "warn");
            end

    
        


%% Plot S Parameters from Network Extraction

        if Plot_S_from_extracted_M_matrix && current_trial_successful
            step_size = (freq_end-freq_start)/steps;
            
            R = zeros(N, N);
            R(1,1) = RS;
            R(end, end) = RL;
            
            S11_M_matrix = zeros(1, steps + 1);
            S21_M_matrix = zeros(1, steps + 1);
            
            WB3 = waitbar(0,'Calculating S11 and S21 from extracted M Matrix ....');

            delta = center_freq./(Bandwidth.*Q);

            for f = freq_start : step_size : freq_end
                waitbar((f-freq_start)/(freq_end-freq_start), WB3,'Calculating S11 and S21 from extracted M Matrix ....');

                lambda = center_freq/Bandwidth*(f/center_freq-center_freq/f);

                A_matrix = diag(lambda - 1i*delta)*eye(N) - 1i*R + scaled_M_matrix;
                A_matrix_inv = A_matrix^(-1);
                
                S11_M_matrix(round((f - freq_start)/step_size + 1)) = abs(1 + 2*1i*RS*A_matrix_inv(1,1));
                S21_M_matrix(round((f - freq_start)/step_size + 1)) = abs(-2*1i*sqrt(RS*RL)*A_matrix_inv(N,1));
            end
        
            close(WB3)
            
            freq = linspace(freq_start, freq_end, steps + 1);
            
            figure;
            ref = plot(freq, 20*log10(S11_M_matrix));
            hold on
            trans = plot(freq, 20*log10(S21_M_matrix));
            hold off
            
            legend([ref, trans], "S11", "S21")
            xlabel("Frequency (Hz)")
            ylabel("dB")
            title("Extracted M Matrix S11 S21 vs Frequency")
        end
        catch
            %close(WB2)
            msgbox("Network extraction failed! Please try another coupling connection.", "Warning", "warn");
        end
    
    else
    end
end



%% Automatic Find Solution

if Find_Extraction_Solution_All || Find_Extraction_Solution_Fast || Find_Extraction_Solution_Ultra_Fast

    current_trial_successful = 0;
    num_solution_found = 0;
    found_cross_connection_matrix = zeros(N, N, target_num_solution);
    found_B_enable = zeros(target_num_solution, N);
    found_M_matrix = zeros(N, N, target_num_solution);
    found_Extracted_C = zeros(target_num_solution, N);
    found_Extracted_B = zeros(target_num_solution, N);

    All_possible_B_enable = reshape((dec2bin(0:2^N-1)-'0'),[], N);

    if isEven(N)
        num_choices_cross_connections = (ceil(N/2) - 1)*3;
    else
        num_choices_cross_connections = (ceil(N/2) - 1)*3 - 2;
    end
    
    All_possible_cross_connections = reshape((dec2bin(0:2^num_choices_cross_connections-1)-'0'),[], num_choices_cross_connections);

    TOTAL_NUM_TRIAL_CROSS_CONNECTION = size(All_possible_cross_connections, 1);

    sum_TZ = 0;
    for i = 1:1:N
        if TZ(i) ~= inf
            sum_TZ = sum_TZ + TZ(i);
        end
    end

    if abs(sum_TZ) < 1e-9
        symmetric_TZ = true;
    else
        symmetric_TZ = false;
    end

    if Find_Extraction_Solution_Fast || Find_Extraction_Solution_Ultra_Fast
        TOTAL_NUM_TRIAL_B = 1;
    else
        TOTAL_NUM_TRIAL_B = size(All_possible_B_enable, 1);
    end

    WB4 = waitbar(0,'Finding network extraction solution ....');

    for CURRENT_NUM_TRIAL_B = 1:1:TOTAL_NUM_TRIAL_B

        if Find_Extraction_Solution_Fast || Find_Extraction_Solution_Ultra_Fast
            if symmetric_TZ
                B_enable = zeros(1, N);
            else
                B_enable = ones(1, N);
            end
        else
            B_enable = All_possible_B_enable(CURRENT_NUM_TRIAL_B, :);
        end
        
        if Find_Extraction_Solution_Ultra_Fast
            cross_connection_matrix_history = zeros(N, N, TOTAL_NUM_TRIAL_CROSS_CONNECTION);
        end

        for CURRENT_NUM_TRIAL_CROSS_CONNECTION = 1:1:TOTAL_NUM_TRIAL_CROSS_CONNECTION

            waitbar((TOTAL_NUM_TRIAL_CROSS_CONNECTION*(CURRENT_NUM_TRIAL_B - 1) + CURRENT_NUM_TRIAL_CROSS_CONNECTION)/(TOTAL_NUM_TRIAL_CROSS_CONNECTION*TOTAL_NUM_TRIAL_B), WB4, 'Finding network extraction solution ....');

            M_matrix = zeros(N, N);
            
            TOTAL_NUM_EXTRACTION = 0;
            working_node = 0;
            ending_node = N + 1;
           
            cross_connection_matrix = generate_cross_connection_matrix(All_possible_cross_connections(CURRENT_NUM_TRIAL_CROSS_CONNECTION, :));
            
            repeated = 0;
            if Find_Extraction_Solution_Ultra_Fast && target_num_solution ~= 1
                repeated = check_repetition(cross_connection_matrix_history);
                if ~repeated
                    cross_connection_matrix_history = populate(cross_connection_matrix_history);
                end
            end

            if ~repeated
                remaining_cross_connection_matrix = cross_connection_matrix;
    
                current_sweep_valid = 1;
                num_of_elements = ceil(N/2) - 1;
                for element = 1:1:num_of_elements
                    if cross_connection_matrix(element, N - element) == 1 && cross_connection_matrix(element+1, N - element + 1) == 1
                        current_sweep_valid = 0;
                        break
                    end
                end
        
    
                if current_sweep_valid
                    try
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
                            if CURRENT_NUM_EXTRACTION == TOTAL_NUM_EXTRACTION
                                [A, B, C, D, P] = parallel_INV_extraction(A, B, C, D, P);
                                if abs(M_matrix(working_node, ending_node) + 1) < 1e-9
                                    M_matrix(working_node, ending_node) = 1;
                                    M_matrix(ending_node, working_node) = 1;
                                end
                            elseif Extracted_C(working_node) == 0
                                [Extracted_C, A, B, C, D, P] = C_extraction(Extracted_C, A, B, C, D, P);
                            elseif B_enable(working_node) == 1 && Extracted_B(working_node) == 0
                                [Extracted_B, A, B, C, D, P] = B_extraction(Extracted_B, A, B, C, D, P);
                            elseif Extracted_C(ending_node) == 0
                                [A, B, C, D, P] = reverse(A, B, C, D, P);
                                [working_node, ending_node] = swap(working_node, ending_node);
                                [Extracted_C, A, B, C, D, P] = C_extraction(Extracted_C, A, B, C, D, P);
                            elseif B_enable(ending_node) == 1 && Extracted_B(ending_node) == 0
                                [A, B, C, D, P] = reverse(A, B, C, D, P);
                                [working_node, ending_node] = swap(working_node, ending_node);
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
                    
    
                        current_trial_successful = 0;
    
                        if is_C_full(Extracted_C)
                            if is_ABCD_finished(A, B, C, D)
                                if is_response_identical(Extracted_B, Extracted_C, normalized_freq_start, normalized_freq_end, steps, S11_polynomial, S21_polynomial)
                                    current_trial_successful = 1;
                                end
                            end
                        end
                        
    
                        if current_trial_successful
                            num_solution_found = num_solution_found + 1;
                            found_cross_connection_matrix(:, :, num_solution_found) = cross_connection_matrix;
                            found_B_enable(num_solution_found, :) = B_enable;
                            found_M_matrix(:, :, num_solution_found) = M_matrix;
                            found_Extracted_C(num_solution_found, :) = Extracted_C;
                            found_Extracted_B(num_solution_found, :) = Extracted_B;
                        end
    
                    catch
                    end
                end
            end
            if num_solution_found == target_num_solution
                break;
            end
        end

        if num_solution_found == target_num_solution
            break;
        end
    end
    close(WB4)


    found_cross_connection_matrix(:, :, num_solution_found + 1:end) = [];
    found_B_enable(num_solution_found + 1:end, :) = [];
    found_M_matrix(:, :, num_solution_found + 1:end) = [];
    found_Extracted_C(num_solution_found + 1:end, :) = [];
    found_Extracted_B(num_solution_found + 1:end, :) = [];


    if num_solution_found ~= 0
        if num_solution_found == 1
            msgbox("Successfully found " + num_solution_found + " network extraction solution!", "Successful");
        else
            msgbox("Successfully found " + num_solution_found + " network extraction solutions!", "Successful");
        end
        
        for current_solution = 1:1:num_solution_found
            scaled_M_matrix = zeros(N,N);
            scaled_Extracted_C = ones(1, N);
            scaled_Extracted_B = zeros(1, N);
        
        
            for i = 1:1:N
                scaled_Extracted_B(i) = found_Extracted_B(current_solution, i)/found_Extracted_C(current_solution, i);
            end
           
    
            for i = 1:1:N
                for j = 1:1:N
                    if j == i
                        scaled_M_matrix(i,j) = scaled_Extracted_B(j);
                    elseif found_M_matrix(i,j, current_solution) ~= 0
                        scaled_M_matrix(i,j) = found_M_matrix(i,j, current_solution)/sqrt(found_Extracted_C(current_solution, i)*found_Extracted_C(current_solution, j));
                    end
                end
            end


        
            SL_scaled_M_matrix = zeros(N+2, N+2);
            SL_scaled_M_matrix(2:end-1, 2:end-1) = scaled_M_matrix;
            SL_scaled_M_matrix(1, 2) = 1/sqrt(found_Extracted_C(current_solution,1));
            SL_scaled_M_matrix(2, 1) = 1/sqrt(found_Extracted_C(current_solution,1));
            SL_scaled_M_matrix(end-1, end) = 1/sqrt(found_Extracted_C(current_solution,end));
            SL_scaled_M_matrix(end, end-1) = 1/sqrt(found_Extracted_C(current_solution,end));
        
        
            RS = SL_scaled_M_matrix(1, 2)^2;
            RL = SL_scaled_M_matrix(end-1, end)^2;
    
            step_size = (freq_end-freq_start)/steps;
            
            R = zeros(N, N);
            R(1,1) = RS;
            R(end, end) = RL;
            
            S11_M_matrix = zeros(1, steps + 1);
            S21_M_matrix = zeros(1, steps + 1);
            
            WB5 = waitbar(0,'Calculating S11 and S21 from extracted M Matrix ....');

            delta = center_freq./(Bandwidth.*Q);
            
            for f = freq_start : step_size : freq_end
                waitbar((f-freq_start)/(freq_end-freq_start), WB5,'Calculating S11 and S21 from extracted M Matrix ....');
            
                lambda = center_freq/Bandwidth*(f/center_freq-center_freq/f);

                A_matrix = diag(lambda - 1i*delta)*eye(N) - 1i*R + scaled_M_matrix;
                A_matrix_inv = A_matrix^(-1);
                
                S11_M_matrix(round((f - freq_start)/step_size + 1)) = abs(1 + 2*1i*RS*A_matrix_inv(1,1));
                S21_M_matrix(round((f - freq_start)/step_size + 1)) = abs(-2*1i*sqrt(RS*RL)*A_matrix_inv(N,1));
            end
        
            close(WB5)
            
            freq = linspace(freq_start, freq_end, steps + 1);
            
            figure;
            ref = plot(freq, 20*log10(S11_M_matrix));
            hold on
            trans = plot(freq, 20*log10(S21_M_matrix));
            hold off
            
            legend([ref, trans], "S11", "S21")
            xlabel("Frequency (Hz)")
            ylabel("dB")
            title("Solution " + current_solution + " Extracted M Matrix S11 S21 vs Frequency")
        end
    else
        msgbox("No network extraction solution was found!", "Failed");
    end
end
    


clearvars A_matrix A_matrix_inv All_possible_B_enable All_possible_cross_connections CURRENT_NUM_EXTRACTION CURRENT_NUM_TRIAL_B ...
    CURRENT_NUM_TRIAL_CROSS_CONNECTION current_solution current_trial_successful Enable_ABCDP_simplification Enable_Network_Extraction...
    ending_connection ending_node failed Find_Extraction_Solution_All Find_Extraction_Solution_Fast i j lambda limit_noise_suppressor...
    num_coupling num_FIR num_of_elements num_of_finite_TZ Plot_S_from_polynomials Plot_S_from_extracted_M_matrix ...
    WB0 WB1 WB2 WB3 WB4 WB5 ref trans element remaining_cross_connection_matrix R RS RL S11_M_matrix S11_polynomial ...
    S21_M_matrix S21_polynomial step_size steps TOTAL_NUM_EXTRACTION TOTAL_NUM_TRIAL_B TOTAL_NUM_TRIAL_CROSS_CONNECTION ...
    valid w s working_connection working_node target_num_solution N n freq freq_end freq_start center_freq current_sweep_valid...
    f finite_TZ Network_Extraction_Force_Ending_With_Cross_Coupling normalized_freq_current normalized_freq_end normalized_freq_start...
    num_choices_cross_connections Polynomial_Solver round_to_decimal_places Enable_UV_simplification cross_connection_matrix_history delta...
    repeated sum_TZ symmetric_TZ










%% Functions


function array_done = pad2N_inf(array, N)
    array_done = zeros(1,N+1);
    length = size(array,2);
    if length>=N+1
        array_done = array(end - N:end);
    else
        array_done(N-length+2:N+1) = array;
        array_done(1:N-length+1) = inf;
    end
end

function array_done = pad2N_for_Q(array)
    global N
    array_done = zeros(1,N);
    length = size(array,2);
    if length>=N
        array_done = array(1:N);
    else
        array_done(1:length) = array;
        array_done(length+1:end) = inf;
    end
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
    
    if isEven(N)
        Prem = PS/epsilon;
    else
        Prem = PS * -1i/epsilon;
    end

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
    % sym2poly(Pin(row, index_last_element))
    % sym2poly(Bin(row, index_last_element))
    M_parallel_current = double(limit(-1*noise_suppress(Pin(row, index_last_element))/noise_suppress(Bin(row, index_last_element)), s, inf));

    A_current_sym = Ain(row, index_last_element);
    B_current_sym = Bin(row, index_last_element);
    C_current_sym = Cin(row, index_last_element) + 2*M_parallel_current*Pin(row, index_last_element) + M_parallel_current^2*Bin(row, index_last_element);
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


function sweeping_cross_connection_matrix = generate_cross_connection_matrix(connection_choise)
    global N

    sweeping_cross_connection_matrix = zeros(N, N);

    for i = 1:1:size(connection_choise, 2)
        if connection_choise(i) == 1
            group_number = ceil(i/3);
            group_index = i - 3*(group_number - 1);
            if group_index == 1 || group_index == 2
                row = group_number;
            else
                row = group_number + 1;
            end
            if group_index == 1 || group_index == 3
                column = N - group_number + 1;
            else
                column = N - group_number;
            end
            sweeping_cross_connection_matrix(row, column) = 1;
            sweeping_cross_connection_matrix(column, row) = 1;
        end
    end
end


function finished_flag = is_ABCD_finished(Ain, Bin, Cin, Din)

    A_current_1 = pad2N(sym2poly(noise_suppress(Ain(1, end))));
    A_current_2 = pad2N(sym2poly(noise_suppress(Ain(2, end))));
    if (A_current_1(end-1) ~= 0 && A_current_2(end-1) ~= 0) || isnan(A_current_1(end)) || isnan(A_current_2(end))
        finished_flag = 0;
    else
        B_current_1 = pad2N(sym2poly(noise_suppress(Bin(1, end))));
        B_current_2 = pad2N(sym2poly(noise_suppress(Bin(2, end))));
        if (B_current_1(end-1) ~= 0 && B_current_2(end-1) ~= 0) || isnan(B_current_1(end)) || isnan(B_current_2(end))
            finished_flag = 0;
        else
            C_current_1 = pad2N(sym2poly(noise_suppress(Cin(1, end))));
            C_current_2 = pad2N(sym2poly(noise_suppress(Cin(2, end))));
            if (C_current_1(end-1) ~= 0 && C_current_2(end-1) ~= 0) || isnan(C_current_1(end)) || isnan(C_current_2(end))
                finished_flag = 0;
            else
                D_current_1 = pad2N(sym2poly(noise_suppress(Din(1, end))));
                D_current_2 = pad2N(sym2poly(noise_suppress(Din(2, end))));
                if (D_current_1(end-1) ~= 0 && D_current_2(end-1) ~= 0) || isnan(D_current_1(end)) || isnan(D_current_2(end))
                    finished_flag = 0;
                else
                    finished_flag = 1;
                end
            end
        end
    end
end


function same_S_flag = is_response_identical(Extracted_B, Extracted_C, normalized_freq_start, normalized_freq_end, steps, S11_polynomial, S21_polynomial)
    global M_matrix N

    scaled_M_matrix = zeros(N,N);
    scaled_Extracted_B = zeros(1, N);

    
    for i = 1:1:N
        scaled_Extracted_B(i) = Extracted_B(i)/Extracted_C(i);
    end
   

    for i = 1:1:N
        for j = 1:1:N
            if j == i
                scaled_M_matrix(i,j) = scaled_Extracted_B(j);
            elseif M_matrix(i,j) ~= 0
                scaled_M_matrix(i,j) = M_matrix(i,j)/sqrt(Extracted_C(i)*Extracted_C(j));
            end
        end
    end


    RS = (1/sqrt(Extracted_C(1)))^2;
    RL = (1/sqrt(Extracted_C(end)))^2;
    
    step_size = (normalized_freq_end-normalized_freq_start)/steps;
    
    R = zeros(N, N);
    R(1,1) = RS;
    R(end, end) = RL;
    
    cumulative_error = 0;

    for f = normalized_freq_start : step_size : normalized_freq_end
        lambda = f;
        
        A_matrix = lambda*eye(N) - 1i*R + scaled_M_matrix;
        A_matrix_inv = A_matrix^(-1);
        
        cumulative_error = cumulative_error + abs(S11_polynomial(round((f - normalized_freq_start)/step_size + 1)) - abs(1 + 2*1i*RS*A_matrix_inv(1,1)));
        cumulative_error = cumulative_error + abs(S21_polynomial(round((f - normalized_freq_start)/step_size + 1)) - abs(-2*1i*sqrt(RS*RL)*A_matrix_inv(N,1)));
    end

    if cumulative_error/steps > 1e-10
        same_S_flag = 0;
    else
        same_S_flag = 1;
    end
end


function C_full_flag = is_C_full(Extracted_C)
    global N

    C_full_flag = 1;
    for i = 1:1:N
        if Extracted_C(i) == 0
            C_full_flag = 0;
            break;
        end
    end
end

function repeated = check_repetition(cross_connection_matrix_history)
    global cross_connection_matrix
    
    history_length = 0;
    for i = 1:1:size(cross_connection_matrix_history, 3) 
        if sum(cross_connection_matrix_history(:,:,i), "all") ~= 0
            history_length = history_length + 1;
        else
            break
        end
    end
    
    repeated = 0;
    if history_length == 0
        repeated = 0;
    else
        for i = 1:1:history_length
            %rotated = rot90(cross_connection_matrix, 2);
            if cross_connection_matrix_history(:,:,i) == rot90(cross_connection_matrix, 2)
                repeated = 1;
            end
        end
    end
end

function cross_connection_matrix_history = populate(cross_connection_matrix_history_in)
    global cross_connection_matrix

    cross_connection_matrix_history = cross_connection_matrix_history_in;

    history_length = 0;
    for i = 1:1:size(cross_connection_matrix_history_in, 3) 
        if sum(cross_connection_matrix_history_in(:,:,i), "all") ~= 0
            history_length = history_length + 1;
        else
            break
        end
    end

    cross_connection_matrix_history(:,:,history_length + 1) = cross_connection_matrix;
end