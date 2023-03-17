%% Filter Setup
N = 3;          % Filter order
RL = 23;        % Filter Return Loss (dB)
TZ = [inf];     % Array of frequencies of transmittion zeros (rad/s)
                % For transmission zeros at infinity, type "inf"

%% Simulation Setup
Solver = "recursive";         % Choose between recursive solver and numerical solver


%% Plotting Setup
Plot_S_from_polynomials = true;  % Enable to plot S11 and S21 from the polynomials
w_start = -3;                 % Start frequency (rad/s)
w_end = 3;                    % End frequency (rad/s)
steps = 1000;                 % Number of steps


%% Solve for P(w) and P(S)
syms s w 

TZ = pad2N_inf(TZ, N-1);
finite_TZ = TZ(TZ ~= inf);
num_of_finite_TZ = size(finite_TZ, 2);

PW_roots = TZ;
PW = poly(PW_roots);
PW = pad2N(PW, N);
PW_sym = poly2sym(PW, w);
PS_sym = subs(PW_sym, w, -1i*s);

if mod(N-num_of_finite_TZ, 2) == 0
    PS_sym = PS_sym * 1;
end

PS = sym2poly(PS_sym);
PS = pad2N(PS, N);

%% Solve for F(w) and F(S)
if Solver == "recursive"
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
    FS = pad2N(FS, N);

elseif Solver == "numerical"
    innerTerm = 0;
    for i = 1 : 1 : num_of_finite_TZ
        innerTerm = innerTerm + acosh((w-1/finite_TZ(i))/(1-w/finite_TZ(i)));
    end
    innerTerm = innerTerm + (N - num_of_finite_TZ)*acosh(w);
    
    CW = @(w) cosh(innerTerm);
    
    PW = poly(TZ);
    PW = pad2N(PW, N);
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
    if mod(N-num_of_finite_TZ, 2) == 0
        PS_sym = PS_sym * 1i;
    end
    
    PS = sym2poly(PS_sym);
    PS = pad2N(PS, N);
    
    FW = poly(FW_roots);
    FW(abs(FW)<1e-10) = 0;
    FW_sym = poly2sym(FW, w);
    FS_sym = subs(FW_sym, w, -1i*s);
    FS = sym2poly(FS_sym);
    FS = pad2N(FS, N);
else
    warndlg('Please choose a valid solver.','Warning');
end

%% Solve for epsilon, E(w) and E(S)
epsilon = double(1/sqrt(10^(RL/10)-1)*subs(PW_sym, w, 1)/subs(FW_sym, w, 1));
reference_epsilon = double(1/sqrt(10^(RL/10)-1));

%ES2 = conv(PS, (PS))/epsilon^2 + conv(FS, (FS));
% INTER1 = conv(PS, conj(PS))/epsilon^2;
% INTER2 = conv(FS, conj(FS));
ES2 = conv(PS, conj(PS))/epsilon^2 + conv(FS, conj(FS));
ES2_sym = poly2sym(ES2, s);
ES2_roots = roots(ES2);
ES_roots = ES2_roots(real(ES2_roots)<0);
%ES_roots = ES2_roots(1:2);

if num_of_finite_TZ == 0
    for k = 1:1:N
        reference_ES_roots(k) = -sinh(1/N*asinh(1/reference_epsilon))*sin((pi/2)*(2*k-1)/N) + 1i * cosh(1/N*asinh(1/reference_epsilon))*cos((pi/2)*(2*k-1)/N);
    end
    ES = poly(reference_ES_roots)
else
    ES = poly(ES_roots)
end

%reference_ES2 = conv(ES, conj(ES));
% reference_ES2_roots = roots(reference_ES2);
ES_sym = poly2sym(ES, s);
EW_sym = subs(ES_sym, s, 1i*w);
EW = sym2poly(EW_sym);




%% Plot S parameters from the polynomials
if Plot_S_from_polynomials
    step_size = (w_end-w_start)/steps;
    
    S11 = zeros(1, steps + 1);
    S21 = zeros(1, steps + 1);
    
    
    B = waitbar(0,'Calculating S11 and S21 from polynomials ....');
    
    for w_current = w_start : step_size : w_end
        waitbar((w_current-w_start)/(w_end - w_start), B,'Calculating S11 and S21 ....');
        S21(round((w_current - w_start)/step_size + 1)) = double(abs(subs(PW_sym, w, w_current)))/(epsilon * double(abs(subs(EW_sym, w, w_current))));
        S11(round((w_current - w_start)/step_size + 1)) = double(abs(subs(FW_sym, w, w_current)))/double(abs(subs(EW_sym, w, w_current)));
    end
    
    close(B)
    
    freq = linspace(w_start, w_end, steps + 1);
    
    figure;
    ref = plot(freq, 20*log10(abs(S11)));
    hold on
    trans = plot(freq, 20*log10(abs(S21)));
    hold off
    
    legend([ref, trans], "S11", "S21")
    xlabel("Normalized Frequency (rad)")
    ylabel("dB")
    title("Polynomials S11 S21 vs Frequency")
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [AS, BS, CS, DS, PS] = EF2ABCD(ES, FS, PS, epsilon, N);
% 
% [AS_1, BS_1, CS_1, DS_1, PS_1] = series_unit_INV_extraction(AS, BS, CS, DS, PS);
% 
% [C1, AS_2, BS_2, CS_2, DS_2, PS_2] = C_extraction(AS_1, BS_1, CS_1, DS_1, PS_1, N);
% 
% [AS_3r, BS_3r, CS_3r, DS_3r, PS_3r] = reverse(AS_2, BS_2, CS_2, DS_2, PS_2);
% 
% [AS_4, BS_4, CS_4, DS_4, PS_4] = series_unit_INV_extraction(AS_3r, BS_3r, CS_3r, DS_3r, PS_3r);
% 
% [C4, AS_5, BS_5, CS_5, DS_5, PS_5] = C_extraction(AS_4, BS_4, CS_4, DS_4, PS_4, N);
% 
% [M14, AS_6, BS_6, CS_6, DS_6, PS_6] = parallel_INV_extraction(AS_5, BS_5, CS_5, DS_5, PS_5, N);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function array_done = pad2N(array, N)
    array_done = zeros(1,N+1);
    length = size(array,2);
    array_done(N-length+2:N+1) = array;
    array_done(1:N-length+1) = 0;
end

function array_done = pad2N_inf(array, N)
    array_done = zeros(1,N+1);
    length = size(array,2);
    array_done(N-length+2:N+1) = array;
    array_done(1:N-length+1) = inf;
end

