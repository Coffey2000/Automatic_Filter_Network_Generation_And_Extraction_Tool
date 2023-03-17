N = 6;

M = zeros(8,8);
M(1,2) = 0.81664;
M(2,1) = 0.81664;
M(2,3) = 0.586767;
M(3,2) = 0.586767;
M(3,4) = 0.547577;
M(4,3) = 0.547577;
M(4,5) = 0.536432;
M(5,4) = 0.536432;
M(5,6) = 0.503952;
M(6,5) = 0.503952;
M(5,8) = -0.209083;
M(8,5) = -0.209083;
M(6,7) = 0.723611;
M(7,6) = 0.723611;
M(7,8) = 0.789421;
M(8,7) = 0.789421;

RS = 0.987513^2;
RL = 0.987513^2;
BW = (90 + 10*N)*10^6;
f0 = 3.8e9;

R = zeros(8,8);
R(1,1) = RS;
R(8,8) = RL;

S11 = zeros(1, 600001);
S21 = zeros(1, 600001);


for f = 3.5e9 : 1000 : 4.1e9
lambda = f0/BW*(f/f0-f0/f);

A = lambda*eye(8) - 1i*R + M;
A_inv = A^(-1);

S11((f-3.5e9)/1000 + 1) = 1 + 2*1i*RS*A_inv(1,1);
S21((f-3.5e9)/1000 + 1) = -2*1i*sqrt(RS*RL)*A_inv(8,1);
end

freq = linspace(3.5e9, 4.1e9, 600001);

figure;
ref = plot(freq, 20*log10(abs(S11)));
hold on
trans = plot(freq, 20*log10(abs(S21)));
hold off

legend([ref, trans], "S11", "S21")
xlabel("Frequency (Hz)")
ylabel("dB")
title("Q4a S11 S21 vs Frequency")





M2 = M + (0.1+N/100)*eye(8);

S11_2 = zeros(1, 800001);
S21_2 = zeros(1, 800001);

for f = 3.4e9 : 1000 : 4.2e9
lambda = f0/BW*(f/f0-f0/f);

A = lambda*eye(8) - 1i*R + M2;
A_inv = A^(-1);

S11_2((f-3.4e9)/1000 + 1) = 1 + 2*1i*RS*A_inv(1,1);
S21_2((f-3.4e9)/1000 + 1) = -2*1i*sqrt(RS*RL)*A_inv(8,1);
end

freq_2 = linspace(3.4e9, 4.2e9, 800001);

figure;
ref_2 = plot(freq_2, 20*log10(abs(S11_2)));
hold on
trans_2 = plot(freq_2, 20*log10(abs(S21_2)));
hold off

legend([ref_2, trans_2], "S11", "S21")
xlabel("Frequency (Hz)")
ylabel("dB")
title("Q4b S11 S21 vs Frequency")


figure;
ref_2 = plot(freq_2, 20*log10(abs(S11_2)), "--", "color", "red");
hold on
trans_2 = plot(freq_2, 20*log10(abs(S21_2)), "--", "color", "blue");
hold on
ref = plot(freq, 20*log10(abs(S11)), "color", "red");
hold on
trans = plot(freq, 20*log10(abs(S21)), "color", "blue");
hold off

legend([ref_2, trans_2, ref, trans], "S11 shifted", "S21 shifted", "S11", "S21")
xlabel("Frequency (Hz)")
ylabel("dB")
title("Q4b S11 S21 vs Frequency")




Q = 3000;
delta = f0/(BW*Q);

S11 = zeros(1, 800001);
S21 = zeros(1, 800001);

for f = 3.4e9 : 1000 : 4.2e9
lambda = f0/BW*(f/f0-f0/f);

A = (lambda - 1i*delta)*eye(8) - 1i*R + M2;
A_inv = A^(-1);

S11((f-3.4e9)/1000 + 1) = 1 + 2*1i*RS*A_inv(1,1);
S21((f-3.4e9)/1000 + 1) = -2*1i*sqrt(RS*RL)*A_inv(8,1);
end

freq = linspace(3.4e9, 4.2e9, 800001);

figure;
ref = plot(freq, 20*log10(abs(S11)));
hold on
trans = plot(freq, 20*log10(abs(S21)));
hold off

legend([ref, trans], "S11", "S21")
xlabel("Frequency (Hz)")
ylabel("dB")
title("Q4c S11 S21 vs Frequency")
