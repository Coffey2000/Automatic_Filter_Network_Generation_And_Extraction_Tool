N = 4;
epsilon = 1.1548;
syms s w

global limit_noise_suppressor
limit_noise_suppressor = 1e-3;

ES = [1, 2.4015-1i*0.7591, 3.6706-1i*2.1950, 2.4874-1i*3.6255, -0.1268-1i*2.0658];
FS = [1.0, -1i*0.7591, 0.7869, -1i*0.5432, 0.0208];
PS = [-1, 1i*3.1299, 2.3899];

[AS, BS, CS, DS, PS] = EF2ABCD(ES, FS, PS, epsilon, N);

[AS_1, BS_1, CS_1, DS_1, PS_1] = series_unit_INV_extraction(AS, BS, CS, DS, PS);

[C1, AS_2, BS_2, CS_2, DS_2, PS_2] = C_extraction(AS_1, BS_1, CS_1, DS_1, PS_1, N);

[B1, AS_3, BS_3, CS_3, DS_3, PS_3] = B_extraction(AS_2, BS_2, CS_2, DS_2, PS_2, N);


[AS_3r, BS_3r, CS_3r, DS_3r, PS_3r] = reverse(AS_3, BS_3, CS_3, DS_3, PS_3);

[AS_4r, BS_4r, CS_4r, DS_4r, PS_4r] = series_unit_INV_extraction(AS_3r, BS_3r, CS_3r, DS_3r, PS_3r);

[C4, AS_5r, BS_5r, CS_5r, DS_5r, PS_5r] = C_extraction(AS_4r, BS_4r, CS_4r, DS_4r, PS_4r, N);

[B4, AS_6r, BS_6r, CS_6r, DS_6r, PS_6r] = B_extraction(AS_5r, BS_5r, CS_5r, DS_5r, PS_5r, N);

[M14, AS_7r, BS_7r, CS_7r, DS_7r, PS_7r] = parallel_INV_extraction(AS_6r, BS_6r, CS_6r, DS_6r, PS_6r, N);


[AS_7, BS_7, CS_7, DS_7, PS_7] = reverse(AS_7r, BS_7r, CS_7r, DS_7r, PS_7r);

[AS_8, BS_8, CS_8, DS_8, PS_8] = series_unit_INV_extraction(AS_7, BS_7, CS_7, DS_7, PS_7);

[C2, AS_9, BS_9, CS_9, DS_9, PS_9] = C_extraction(AS_8, BS_8, CS_8, DS_8, PS_8, N);

[B2, AS_10, BS_10, CS_10, DS_10, PS_10] = B_extraction(AS_9, BS_9, CS_9, DS_9, PS_9, N);

[M24, AS_11, BS_11, CS_11, DS_11, PS_11] = parallel_INV_extraction(AS_10, BS_10, CS_10, DS_10, PS_10, N);


[AS_11r, BS_11r, CS_11r, DS_11r, PS_11r] = reverse(AS_11, BS_11, CS_11, DS_11, PS_11);

[AS_12r, BS_12r, CS_12r, DS_12r, PS_12r] = series_unit_INV_extraction(AS_11r, BS_11r, CS_11r, DS_11r, PS_11r);

[C3, AS_13r, BS_13r, CS_13r, DS_13r, PS_13r] = C_extraction(AS_12r, BS_12r, CS_12r, DS_12r, PS_12r, N);

[B3,  AS_14r, BS_14r, CS_14r, DS_14r, PS_14r] = B_extraction(AS_13r, BS_13r, CS_13r, DS_13r, PS_13r, N);

[M23, AS_15r, BS_15r, CS_15r, DS_15r, PS_15r] = parallel_INV_extraction(AS_14r, BS_14r, CS_14r, DS_14r, PS_14r, N);





function array_done = pad2N(array, N)
    array_done = zeros(1,N+1);
    length = size(array,2);
    if length>N+1
        array_done = array(end - N:end);
    else
        array_done(N-length+2:N+1) = array;
        array_done(1:N-length+1) = 0;
    end
end

function [AS, BS, CS, DS, Prem] = EF2ABCD(ES,FS,PS, epsilon, N)
    ES = pad2N(ES, N);
    FS = pad2N(FS, N);
    PS = pad2N(PS, N);
    
    
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
end

function [Arem, Brem, Crem, Drem, Prem] = series_unit_INV_extraction(AS, BS, CS, DS, PS)
    Arem = -1i*CS;
    Brem = -1i*DS;
    Crem = -1i*AS;
    Drem = -1i*BS;
    Prem = PS;
end


function [Cap, Arem, Brem, Crem, Drem, Prem] = C_extraction(A, B, C, D, P, N)
    global limit_noise_suppressor
    syms s
    D_noiseless = D;
    D_noiseless(abs(D_noiseless)<limit_noise_suppressor) = 0;
    B_noiseless = B;
    B_noiseless(abs(B_noiseless)<limit_noise_suppressor) = 0;
    D_noiseless_sym = poly2sym(D_noiseless, s);
    B_noiseless_sym = poly2sym(B_noiseless, s);
    Cap = double(limit(D_noiseless_sym/(s*B_noiseless_sym), s, inf));

    Arem = A;
    Brem = B;
    Crem = pad2N(sym2poly(poly2sym(C, s) - s*Cap*poly2sym(A, s)), N);
    Drem = pad2N(sym2poly(poly2sym(D, s) - s*Cap*poly2sym(B, s)), N);
    Prem = P;
end



function [B0, Arem, Brem, Crem, Drem, Prem] = B_extraction(A, B, C, D, P, N)
    global limit_noise_suppressor
    syms s
    D_noiseless = D;
    D_noiseless(abs(D_noiseless)<limit_noise_suppressor) = 0;
    B_noiseless = B;
    B_noiseless(abs(B_noiseless)<limit_noise_suppressor) = 0;
    D_noiseless_sym = poly2sym(D_noiseless, s);
    B_noiseless_sym = poly2sym(B_noiseless, s);
    B0 = imag(double(limit(D_noiseless_sym/B_noiseless_sym, s, inf)));

    Arem = A;
    Brem = B;
    Crem = pad2N(sym2poly(poly2sym(C, s) - 1i*B0*poly2sym(A, s)), N);
    Drem = pad2N(sym2poly(poly2sym(D, s) - 1i*B0*poly2sym(B, s)), N);
    Prem = P;
end

function [Arev, Brev, Crev, Drev, Prev] = reverse(A,B,C,D,P)
    Arev = D;
    Brev = B;
    Crev = C;
    Drev = A;
    Prev = P;
end

function [Mparallel, Arem, Brem, Crem, Drem, Prem] = parallel_INV_extraction(A, B, C, D, P, N)
    global limit_noise_suppressor
    syms s
    P_noiseless = P;
    P_noiseless(abs(P_noiseless)<limit_noise_suppressor) = 0;
    B_noiseless = B;
    B_noiseless(abs(B_noiseless)<limit_noise_suppressor) = 0;
    P_noiseless_sym = poly2sym(P_noiseless, s);
    B_noiseless_sym = poly2sym(B_noiseless, s);
    Mparallel = double(limit(-P_noiseless_sym/B_noiseless_sym, s, inf));

    Arem = A;
    Brem = B;
    Crem = pad2N(sym2poly(poly2sym(C, s) + 2*Mparallel*poly2sym(P, s)+ Mparallel^2*poly2sym(B, s)), N);
    Drem = D;
    Prem = pad2N(sym2poly(poly2sym(P, s) + Mparallel*poly2sym(B, s)), N);
end







