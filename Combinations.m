global N
N = 4;
All_possible_B_enable = reshape((dec2bin(0:2^N-1)-'0'),[], N);
%sweeping_cross_connection_matrix = zeros(N, N);
if isEven(N)
    num_choices_cross_connections = (ceil(N/2) - 1)*3;
else
    num_choices_cross_connections = (ceil(N/2) - 1)*3 - 2;
end

All_possible_cross_connections = reshape((dec2bin(0:2^num_choices_cross_connections-1)-'0'),[], num_choices_cross_connections);
%sweeping_cross_connection_matrix = generate_cross_connection_matrix(All_possible_cross_connections(2,:));

sweeping_cross_connection_matrix = generate_cross_connection_matrix([1 0 0]);

function sweeping_cross_connection_matrix = generate_cross_connection_matrix(connection_choise)
    global N

    sweeping_cross_connection_matrix = zeros(N, N);

    for i = 1:1:size(connection_choise, 2)
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
        if connection_choise(i) == 1
            sweeping_cross_connection_matrix(row, column) = 1;
            sweeping_cross_connection_matrix(column, row) = 1;
        end
    end
end














function boolean1 = isEven(N)
    boolean1 = mod(N, 2) == 0;
end