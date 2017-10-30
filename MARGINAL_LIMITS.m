
load('analytical2.mat')


%C_limits_2015 = [c1_min 0.0001 c1_max 0.0001; c2_min 0.0005 c2_max 0.00005; c3_min 0.0001 c3_max 0.00005; c4_min 0.00005 c4_max 0.00005; c5_min 0.0005 c5_max 0.0001];
%these should be < abs(4pi)
%both O^8 (strong coupling) : don't have exp senitivity yet (maybe with more Lum.)

%[C1_min_2015,C2_min_2015,C3_min_2015,C4_min_2015,C5_min_2015,C1_min_tol_2015,C2_min_tol_2015,C3_min_tol_2015,C4_min_tol_2015,C5_min_tol_2015] = get_lower_limit(sigma_exp_SM,-35,10^-5,sig_i,sig_ij);
%[C1_max_2015,C2_max_2015,C3_max_2015,C4_max_2015,C5_max_2015,C1_max_tol_2015,C2_max_tol_2015,C3_max_tol_2015,C4_max_tol_2015,C5_max_tol_2015] = get_upper_limit(sigma_exp_SM,35,10^-5,sig_i,sig_ij);


%[C1_min_2016,C2_min_2016,C3_min_2016,C4_min_2016,C5_min_2016,C1_min_tol_2016,C2_min_tol_2016,C3_min_tol_2016,C4_min_tol_2016,C5_min_tol_2016]=get_lower_limit(sigma_exp_SM_2016,-35,10^-5,sig_i,sig_ij);
%[C1_max_2016,C2_max_2016,C3_max_2016,C4_max_2016,C5_max_2016,C1_max_tol_2016,C2_max_tol_2016,C3_max_tol_2016,C4_max_tol_2016,C5_max_tol_2016]=get_upper_limit(sigma_exp_SM_2016,35,10^-5,sig_i,sig_ij);



%C_limits_2015 = [C1_min_2015 C1_min_tol_2015 C1_max_2015 C1_max_tol_2015;C2_min_2015 C2_min_tol_2015 C2_max_2015 C2_max_tol_2015;C3_min_2015 C3_min_tol_2015 C3_max_2015 C3_max_tol_2015;C4_min_2015 C4_min_tol_2015 C4_max_2015 C4_max_tol_2015;C5_min_2015 C5_min_tol_2015 C5_max_2015 C5_max_tol_2015];
%C_limits_2016 = [C1_min_2016 C1_min_tol_2016 C1_max_2016 C1_max_tol_2016;C2_min_2016 C2_min_tol_2016 C2_max_2016 C2_max_tol_2016;C3_min_2016 C3_min_tol_2016 C3_max_2016 C3_max_tol_2016;C4_min_2016 C4_min_tol_2016 C4_max_2016 C4_max_tol_2016;C5_min_2016 C5_min_tol_2016 C5_max_2016 C5_max_tol_2016];

%Runtime 2016 limits = 1900s = 32min

sigma_SM_29_09_17 = 2.9*sig_SM - sig_SM; % ~ 27fb - constant from SM in equation

% RUN FOR LAMBDA = 1,2,3,4
% Marginal and independent.

[C1_min_29_09_17,C2_min_29_09_17,C3_min_29_09_17,C4_min_29_09_17,C5_min_29_09_17,C1_min_tol_29_09_17,C2_min_tol_29_09_17,C3_min_tol_29_09_17,C4_min_tol_29_09_17,C5_min_tol_29_09_17]=get_lower_limit(sigma_SM_29_09_17,-10,10^-5,sig_i,sig_ij);
%[C1_max_29_09_17,C2_max_29_09_17,C3_max_29_09_17,C4_max_29_09_17,C5_max_29_09_17,C1_max_tol_29_09_17,C2_max_tol_29_09_17,C3_max_tol_29_09_17,C4_max_tol_29_09_17,C5_max_tol_29_09_17]=get_upper_limit(sigma_SM_29_09_17,35,10^-5,sig_i,sig_ij);



C_limits_29_09_17 = [C1_min_29_09_17 C1_min_tol_29_09_17 C1_max_29_09_17 C1_max_tol_29_09_17;C2_min_29_09_17 C2_min_tol_29_09_17 C2_max_29_09_17 C2_max_tol_29_09_17;C3_min_29_09_17 C3_min_tol_29_09_17 C3_max_29_09_17 C3_max_tol_29_09_17;C4_min_29_09_17 C4_min_tol_29_09_17 C4_max_29_09_17 C4_max_tol_29_09_17;C5_min_29_09_17 C5_min_tol_29_09_17 C5_max_29_09_17 C5_max_tol_29_09_17];

C_limits_t1 = [C_limits_29_09_17(:,1) C_limits_29_09_17(:,3)];

disp(strcat("Perturbative limits: [",string(-4*pi),",",string(4*pi),"]"))

disp("Lambda = 1 TeV")
disp(C_limits_t1)
disp("Lambda = 2 TeV")
disp(C_limits_t1.*4)
disp("Lambda = 3 TeV")
disp(C_limits_t1.*9)
disp("Lambda = 4 TeV")
disp(C_limits_t1.*16)



% disp('2015:')
% disp(C_limits_2015)
% disp('2016:')
% disp(C_limits_2016)


function [C1_min,C2_min,C3_min,C4_min,C5_min,C1_tol,C2_tol,C3_tol,C4_tol,C5_tol] = get_lower_limit(cross_threshold,X0,Tolerance,sig_i,sig_ij)% C = investigated c
% scan over other parameters than the one investigated
Y0 = -4*pi;
Y1 = 4*pi;
Ystep = 0.8;
Ylen = size(Y0:Ystep:Y1,2); % length array

C1_min = 0;
C2_min = 0;
C3_min = 0;
C4_min = 0;
C5_min = 0;
C1_tol = Tolerance;
C2_tol = Tolerance;
C3_tol = Tolerance;
C4_tol = Tolerance;
C5_tol = Tolerance;

%scan over parameter that is investigated

%X0 given by user, lowest bound to investigate parameters, same used for
%all 5
X0_original = X0;
X1 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);


%start c1


cross_row = zeros(1,Xlen);
C1_first_zero = Xlen;

disp('starting lower limit search for C1:')
C1_quit = 'false';
C1_iterations = 0;
while strcmp(C1_quit,'false')
    C1_iterations = C1_iterations +1,
    %C1_matrix = ones(1,Xlen*Ylen^4); % preallocate matrix of certain dimension
    %i1 = 0;
    C1_first_zero = Xlen;
    for c1 = Y0:Ystep:Y1
        for c3 = Y0:Ystep:Y1
            for c4 = Y0:Ystep:Y1
                for c5 = Y0:Ystep:Y1
                    i1 = 0;
                    for c2 = X0:Xstep:X1
                        i1 = i1+1;
                        cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
                        cross = cross - cross_threshold;
                        
                        %% Create Binary data
                        if cross >= 0
                            cross = 1;% outside of allowed region
                        else
                            cross = 0;% inside of allowed region
                        end
                        cross_row(i1) = cross;
                    end
                    %% Check if first zero appears before previous first zero
                    
                    if C1_first_zero > 1 % index starts at 1 in matlab
                        for k = 1:Xlen
                            if cross_row(k)==0 && k < C1_first_zero
                                C1_first_zero = k;
                                break
                            end
                        end
                    end
                    
                    
                end
            end
        end
    end
    
    X_array = X0:Xstep:X1,
    
    %{
    C1_rows = C1_rows(1,2:size(C1_rows,2)); % eliminate the first element which is artificial
    C1_nr_rows = size(C1_rows,2), % count of rows that contain a 0
    
    %TEST CURRENT ITERATION
    
    X_array = X0:Xstep:X1,
    
    C1_mtag = double(any(C1_matrix==0,1));%returns logical, 0 = false, 1 = true
    
    C1_first_zero = -1;
    for k = 1:Xlen
        if C1_mtag(k)==1
            C1_first_zero = k;
            break
        end
    end
    %}
    
    C1_first_zero,
    
    %     if C1_nr_rows == 1
    %         C1_min = X_array(C1_first_zero);
    %         C1_tol = Xstep;
    %         C1_quit = 'true';
    %         disp('Only 1 row remaining reached')
    %     else
    if Xstep < Tolerance
        C1_min = X_array(C1_first_zero);
        C1_quit = 'true';
        disp('Tolerance reached')
    end
    
    % PREPARE NEXT ITERATION
    
    % reducing search area LOGIC:                                                   see TEST for first operation
    % if there is a zero in a column of C1_matrix, turn C1_mtag to 0
    % After scanning through C1_matrix, C1_mtag will contain ones where no
    % rows have a 0 and 0 where at least 1 row contains a 0 in that column.
    % The new X0 is chosen at the index before the first 0 in C1_mtag, and
    % X1 at the index after the first 0. Those new values will define a new
    % Xstep with again 10 steps between 0 and X1.
    
    if C1_first_zero < 0
        disp("Lower limit out of bounds, try another starting value X0")
        C1_quit = 'true';
    elseif C1_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    elseif C1_first_zero==Xlen
        X0 = X_array(C1_first_zero-1);
        X1 = X_array(C1_first_zero);
    else
        X0 = X_array(C1_first_zero-1);
        X1 = X_array(C1_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end
%outside while-loop C1

disp('Lower limit for C1')
disp(C1_min)
disp('Tolerance C1')
disp(C1_tol)



%end c1

%{

%start c2

X0 = X0_original;
X1 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);

disp('starting lower limit search for C2:')
C2_quit = 'false';
C2_iterations = 0;
while strcmp(C2_quit,'false')
    C2_iterations = C2_iterations +1,
    C2_matrix = ones(1,Xlen*Ylen^4); % preallocate matrix of certain dimension
    i1 = 0;
    for c3 = Y0:Ystep:Y1
        for c4 = Y0:Ystep:Y1
            for c5 = Y0:Ystep:Y1
                for c1 = Y0:Ystep:Y1
                    for c2 = X0:Xstep:X1
                        i1 = i1+1;
                    cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
                    cross = cross - cross_threshold;
                    if cross >= 0
                        cross = 1;% outside of allowed region
                    else
                        cross = 0;% inside of allowed region
                    end
                    C2_matrix(1,i1) = cross;% row vector containing cross section info in binary, to be transformed into matrix (2D data)
  
                    end
                end
            end
        end
    end
    C2_matrix = vec2mat(C2_matrix,Xlen); % make matrix from the row vector
    
    
    C2_rows = 0; % keep rows from C1_matrix which contain a zero
    for row = 1:size(C2_matrix,1)
        for column = 1:size(C2_matrix,2)
            if C2_matrix(row,column) == 0
                C2_rows = [C2_rows row]; %#ok<*AGROW>
                break
            end
        end
    end

    C2_rows = C2_rows(1,2:size(C2_rows,2)); % eliminate the first element which is artificial
    C2_nr_rows = size(C2_rows,2), % count of rows that contain a 0

    %TEST CURRENT ITERATION
    
    X_array = X0:Xstep:X1,
    
    C2_mtag = double(any(C2_matrix==0,1));%returns logical, 0 = false, 1 = true
    
    C2_first_zero = -1;
    for k = 1:Xlen
        if C2_mtag(k)==1
            C2_first_zero = k;
            break
        end
    end
    
    C2_first_zero,
    
%     if C2_nr_rows == 1
%         C2_min = X_array(C2_first_zero);
%         C2_tol = Xstep;
%         C2_quit = 'true';
%         disp('Only 1 row remaining reached')
    if Xstep < Tolerance
        C2_min = X_array(C2_first_zero);
        C2_quit = 'true';
        disp('Tolerance reached')
    end
    
    % PREPARE NEXT ITERATION
    
    % reducing search area LOGIC:                                                   see TEST for first operation
    % if there is a zero in a column of C1_matrix, turn C1_mtag to 0
    % After scanning through C1_matrix, C1_mtag will contain ones where no
    % rows have a 0 and 0 where at least 1 row contains a 0 in that column.
    % The new X0 is chosen at the index before the first 0 in C1_mtag, and
    % X1 at the index after the first 0. Those new values will define a new
    % Xstep with again 10 steps between 0 and X1.
    
    if C2_first_zero < 0
        disp("Lower limit out of bounds, try another starting value X0")
        C2_quit = 'true';
    elseif C2_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    elseif C2_first_zero==Xlen
        X0 = X_array(C2_first_zero-1);
        X1 = X_array(C2_first_zero);
    else
        X0 = X_array(C2_first_zero-1);
        X1 = X_array(C2_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end

disp('Lower limit for C2')
disp(C2_min)
disp('Tolerance C2')
disp(C2_tol)

%end c2





%start c3

X0 = X0_original
X1 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);

disp('starting lower limit search for C3:')
C3_quit = 'false';
C3_iterations =0;
while strcmp(C3_quit,'false')
    C3_iterations = C3_iterations +1,
    C3_matrix = ones(1,Xlen*Ylen^4); % preallocate matrix of certain dimension
    i1 = 0;
    for c4 = Y0:Ystep:Y1
        for c5 = Y0:Ystep:Y1
            for c1 = Y0:Ystep:Y1
                for c2 = Y0:Ystep:Y1
                    for c3 = X0:Xstep:X1
                        i1 = i1+1;
                    cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
                    cross = cross - cross_threshold;
                    if cross >= 0
                        cross = 1;% outside of allowed region
                    else
                        cross = 0;% inside of allowed region
                    end
                    C3_matrix(1,i1) = cross;% row vector containing cross section info in binary, to be transformed into matrix (2D data)
  
                    end
                end
            end
        end
    end
    C3_matrix = vec2mat(C3_matrix,Xlen); % make matrix from the row vector
    
    
    C3_rows = 0; % keep rows from C1_matrix which contain a zero
    for row = 1:size(C3_matrix,1)
        for column = 1:size(C3_matrix,2)
            if C3_matrix(row,column) == 0
                C3_rows = [C3_rows row]; %#ok<*AGROW>
                break
            end
        end
    end

    C3_rows = C3_rows(1,2:size(C3_rows,2)); % eliminate the first element which is artificial
    C3_nr_rows = size(C3_rows,2), % count of rows that contain a 0

    %TEST CURRENT ITERATION
    
    X_array = X0:Xstep:X1,
    
    C3_mtag = double(any(C3_matrix==0,1));%returns logical, 0 = false, 1 = true
    
    C3_first_zero = -1;
    for k = 1:Xlen
        if C3_mtag(k)==1
            C3_first_zero = k;
            break
        end
    end
    
    C3_first_zero,
    
%     if C3_nr_rows == 1
%         C3_min = X_array(C3_first_zero);
%         C3_tol = Xstep;
%         C3_quit = 'true';
%         disp('Only 1 row remaining reached')
    if Xstep < Tolerance
        C3_min = X_array(C3_first_zero);
        C3_quit = 'true';
        disp('Tolerance reached')
    end
    
    % PREPARE NEXT ITERATION
    
    % reducing search area LOGIC:                                                   see TEST for first operation
    % if there is a zero in a column of C1_matrix, turn C1_mtag to 0
    % After scanning through C1_matrix, C1_mtag will contain ones where no
    % rows have a 0 and 0 where at least 1 row contains a 0 in that column.
    % The new X0 is chosen at the index before the first 0 in C1_mtag, and
    % X1 at the index after the first 0. Those new values will define a new
    % Xstep with again 10 steps between 0 and X1.
    
    if C3_first_zero < 0
        disp("Lower limit out of bounds, try another starting value X0")
        C3_quit = 'true';
    elseif C3_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    elseif C3_first_zero==Xlen
        X0 = X_array(C3_first_zero-1);
        X1 = X_array(C3_first_zero);
    else
        X0 = X_array(C3_first_zero-1);
        X1 = X_array(C3_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end

disp('Lower limit for C3')
disp(C3_min)
disp('Tolerance C3')
disp(C3_tol)


%end c3




%start c4


X0 = X0_original;
X1 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);

disp('starting lower limit search for C4:')
C4_quit = 'false';
C4_iterations = 0;
while strcmp(C4_quit,'false')
    C4_iterations = C4_iterations +1,
    C4_matrix = ones(1,Xlen*Ylen^4); % preallocate matrix of certain dimension
    i1 = 0;
    for c5 = Y0:Ystep:Y1
        for c1 = Y0:Ystep:Y1
            for c2 = Y0:Ystep:Y1
                for c3 = Y0:Ystep:Y1
                    for c4 = X0:Xstep:X1
                        i1 = i1+1;
                    cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
                    cross = cross - cross_threshold;
                    if cross >= 0
                        cross = 1;% outside of allowed region
                    else
                        cross = 0;% inside of allowed region
                    end
                    C4_matrix(1,i1) = cross;% row vector containing cross section info in binary, to be transformed into matrix (2D data)
  
                    end
                end
            end
        end
    end
    C4_matrix = vec2mat(C4_matrix,Xlen); % make matrix from the row vector
    
    
    C4_rows = 0; % keep rows from C1_matrix which contain a zero
    for row = 1:size(C4_matrix,1)
        for column = 1:size(C4_matrix,2)
            if C4_matrix(row,column) == 0
                C4_rows = [C4_rows row]; %#ok<*AGROW>
                break
            end
        end
    end

    C4_rows = C4_rows(1,2:size(C4_rows,2)); % eliminate the first element which is artificial
    C4_nr_rows = size(C4_rows,2), % count of rows that contain a 0

    %TEST CURRENT ITERATION
    
    X_array = X0:Xstep:X1,
    
    C4_mtag = double(any(C4_matrix==0,1));%returns logical, 0 = false, 1 = true
    
    C4_first_zero = -1;
    for k = 1:Xlen
        if C4_mtag(k)==1
            C4_first_zero = k;
            break
        end
    end
    
    C4_first_zero,
    
%     if C4_nr_rows == 1
%         C4_min = X_array(C4_first_zero);
%         C4_tol = Xstep;
%         C4_quit = 'true';
%         disp('Only 1 row remaining reached')
    if Xstep < Tolerance
        C4_min = X_array(C4_first_zero);
        C4_quit = 'true';
        disp('Tolerance reached')
    end
    
    % PREPARE NEXT ITERATION
    
    % reducing search area LOGIC:                                                   see TEST for first operation
    % if there is a zero in a column of C1_matrix, turn C1_mtag to 0
    % After scanning through C1_matrix, C1_mtag will contain ones where no
    % rows have a 0 and 0 where at least 1 row contains a 0 in that column.
    % The new X0 is chosen at the index before the first 0 in C1_mtag, and
    % X1 at the index after the first 0. Those new values will define a new
    % Xstep with again 10 steps between 0 and X1.
    
    if C4_first_zero < 0
        disp("Lower limit out of bounds, try another starting value X0")
        C4_quit = 'true';
    elseif C4_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    elseif C4_first_zero==Xlen
        X0 = X_array(C4_first_zero-1);
        X1 = X_array(C4_first_zero);
    else
        X0 = X_array(C4_first_zero-1);
        X1 = X_array(C4_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end

disp('Lower limit for C4')
disp(C4_min)
disp('Tolerance C4')
disp(C4_tol)


%end c4




%start c5


X0 = X0_original;
X1 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);

disp('starting lower limit search for C5:')
C5_quit = 'false';
C5_iterations = 0;
while strcmp(C5_quit,'false')
    C5_iterations = C5_iterations +1,
    C5_matrix = ones(1,Xlen*Ylen^4); % preallocate matrix of certain dimension
    i1 = 0;
    for c1 = Y0:Ystep:Y1
        for c2 = Y0:Ystep:Y1
            for c3 = Y0:Ystep:Y1
                for c4 = Y0:Ystep:Y1
                    for c5 = X0:Xstep:X1
                        i1 = i1+1;
                    cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
                    cross = cross - cross_threshold;
                    if cross >= 0
                        cross = 1;% outside of allowed region
                    else
                        cross = 0;% inside of allowed region
                    end
                    C5_matrix(1,i1) = cross;% row vector containing cross section info in binary, to be transformed into matrix (2D data)
  
                    end
                end
            end
        end
    end
    C5_matrix = vec2mat(C5_matrix,Xlen); % make matrix from the row vector
    
    
    C5_rows = 0; % keep rows from C1_matrix which contain a zero
    for row = 1:size(C5_matrix,1)
        for column = 1:size(C5_matrix,2)
            if C5_matrix(row,column) == 0
                C5_rows = [C5_rows row]; %#ok<*AGROW>
                break
            end
        end
    end

    C5_rows = C5_rows(1,2:size(C5_rows,2)); % eliminate the first element which is artificial
    C5_nr_rows = size(C5_rows,2), % count of rows that contain a 0

    %TEST CURRENT ITERATION
    
    X_array = X0:Xstep:X1,
    
    C5_mtag = double(any(C5_matrix==0,1));%returns logical, 0 = false, 1 = true
    
    C5_first_zero = -1;
    for k = 1:Xlen
        if C5_mtag(k)==1
            C5_first_zero = k;
            break
        end
    end
    
    C5_first_zero,
    
%     if C5_nr_rows == 1
%         C5_min = X_array(C5_first_zero);
%         C5_tol = Xstep;
%         C5_quit = 'true';
%         disp('Only 1 row remaining reached')
    if Xstep < Tolerance
        C5_min = X_array(C5_first_zero);
        C5_quit = 'true';
        disp('Tolerance reached')
    end
    
    % PREPARE NEXT ITERATION
    
    % reducing search area LOGIC:                                                   see TEST for first operation
    % if there is a zero in a column of C1_matrix, turn C1_mtag to 0
    % After scanning through C1_matrix, C1_mtag will contain ones where no
    % rows have a 0 and 0 where at least 1 row contains a 0 in that column.
    % The new X0 is chosen at the index before the first 0 in C1_mtag, and
    % X1 at the index after the first 0. Those new values will define a new
    % Xstep with again 10 steps between 0 and X1.
    
    if C5_first_zero < 0
        disp("Lower limit out of bounds, try another starting value X0")
        C5_quit = 'true';
    elseif C5_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    elseif C5_first_zero==Xlen
        X0 = X_array(C5_first_zero-1);
        X1 = X_array(C5_first_zero);
    else
        X0 = X_array(C5_first_zero-1);
        X1 = X_array(C5_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end

disp('Lower limit for C5')
disp(C5_min)
disp('Tolerance C5')
disp(C5_tol)


%end c5

%}




end


function [C1_max,C2_max,C3_max,C4_max,C5_max,C1_tol,C2_tol,C3_tol,C4_tol,C5_tol] = get_upper_limit(cross_threshold,X1,Tolerance,sig_i,sig_ij)% C = investigated c
% scan over other parameters than the one investigated
Y0 = -4*pi;
Y1 = 4*pi;
Ystep = 0.8;
Ylen = size(Y0:Ystep:Y1,2); % length array

C1_max = 0;
C2_max = 0;
C3_max = 0;
C4_max = 0;
C5_max = 0;
C1_tol = Tolerance;
C2_tol = Tolerance;
C3_tol = Tolerance;
C4_tol = Tolerance;
C5_tol = Tolerance;

%scan over parameter that is investigated

%X0 given by user, lowest bound to investigate parameters, same used for
%all 5
X1_original = X1;
X0 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);


%start c1


disp('starting upper limit search for C1:')
C1_quit = 'false';
C1_iterations = 0;
while strcmp(C1_quit,'false')
    C1_iterations = C1_iterations +1,
    C1_matrix = ones(1,Xlen*Ylen^4); % preallocate matrix of certain dimension
    i1 = 0;
    for c2 = Y0:Ystep:Y1
        for c3 = Y0:Ystep:Y1
            for c4 = Y0:Ystep:Y1
                for c5 = Y0:Ystep:Y1
                    for c1 = X0:Xstep:X1
                        i1 = i1+1;
                        cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
                        cross = cross - cross_threshold;
                        if cross >= 0
                            cross = 1;% outside of allowed region
                        else
                            cross = 0;% inside of allowed region
                        end
                        C1_matrix(1,i1) = cross;% row vector containing cross section info in binary, to be transformed into matrix (2D data)
                        
                    end
                end
            end
        end
    end
    C1_matrix = vec2mat(C1_matrix,Xlen); % make matrix from the row vector
    
    
    C1_rows = 0; % keep rows from C1_matrix which contain a zero
    for row = 1:size(C1_matrix,1)
        for column = 1:size(C1_matrix,2)
            if C1_matrix(row,column) == 0
                C1_rows = [C1_rows row]; %#ok<*AGROW>
                break
            end
        end
    end
    
    C1_rows = C1_rows(1,2:size(C1_rows,2)); % eliminate the first element which is artificial
    C1_nr_rows = size(C1_rows,2), % count of rows that contain a 0
    
    %TEST CURRENT ITERATION
    
    X_array = X0:Xstep:X1,
    
    C1_mtag = double(not(any(C1_matrix==0,1)));%returns logical, 0 = false, 1 = true
    
    C1_first_zero = -1;
    for k = 1:Xlen
        if C1_mtag(k)==1
            C1_first_zero = k;
            break
        end
    end
    
    C1_first_zero,
    
    %     if C1_nr_rows == 1
    %         C1_max = X_array(C1_first_zero);
    %         C1_tol = Xstep;
    %         C1_quit = 'true';
    %         disp('Only 1 row remaining reached')
    if Xstep < Tolerance
        C1_max = X_array(C1_first_zero);
        C1_quit = 'true';
        disp('Tolerance reached')
    end
    
    % PREPARE NEXT ITERATION
    
    % reducing search area LOGIC:                                                   see TEST for first operation
    % if there is a zero in a column of C1_matrix, turn C1_mtag to 0
    % After scanning through C1_matrix, C1_mtag will contain ones where no
    % rows have a 0 and 0 where at least 1 row contains a 0 in that column.
    % The new X0 is chosen at the index before the first 0 in C1_mtag, and
    % X1 at the index after the first 0. Those new values will define a new
    % Xstep with again 10 steps between 0 and X1.
    
    if C1_first_zero < 0
        disp("Upper limit out of bounds, try another starting value X1")
        C1_quit = 'true';
    elseif C1_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    elseif C1_first_zero==Xlen
        X0 = X_array(C1_first_zero-1);
        X1 = X_array(C1_first_zero);
    else
        X0 = X_array(C1_first_zero-1);
        X1 = X_array(C1_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end
%outside while-loop C1

disp('Upper limit for C1')
disp(C1_max)
disp('Tolerance C1')
disp(C1_tol)



%end c1


%start c2

X1 = X1_original;
X0 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);

disp('starting upper limit search for C2:')
C2_quit = 'false';
C2_iterations = 0;
while strcmp(C2_quit,'false')
    C2_iterations = C2_iterations +1,
    C2_matrix = ones(1,Xlen*Ylen^4); % preallocate matrix of certain dimension
    i1 = 0;
    for c3 = Y0:Ystep:Y1
        for c4 = Y0:Ystep:Y1
            for c5 = Y0:Ystep:Y1
                for c1 = Y0:Ystep:Y1
                    for c2 = X0:Xstep:X1
                        i1 = i1+1;
                        cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
                        cross = cross - cross_threshold;
                        if cross >= 0
                            cross = 1;% outside of allowed region
                        else
                            cross = 0;% inside of allowed region
                        end
                        C2_matrix(1,i1) = cross;% row vector containing cross section info in binary, to be transformed into matrix (2D data)
                        
                    end
                end
            end
        end
    end
    C2_matrix = vec2mat(C2_matrix,Xlen); % make matrix from the row vector
    
    
    C2_rows = 0; % keep rows from C1_matrix which contain a zero
    for row = 1:size(C2_matrix,1)
        for column = 1:size(C2_matrix,2)
            if C2_matrix(row,column) == 0
                C2_rows = [C2_rows row]; %#ok<*AGROW>
                break
            end
        end
    end
    
    C2_rows = C2_rows(1,2:size(C2_rows,2)); % eliminate the first element which is artificial
    C2_nr_rows = size(C2_rows,2), % count of rows that contain a 0
    
    %TEST CURRENT ITERATION
    
    X_array = X0:Xstep:X1,
    
    C2_mtag = double(not(any(C2_matrix==0,1)));%returns logical, 0 = false, 1 = true
    
    C2_first_zero = -1;
    for k = 1:Xlen
        if C2_mtag(k)==1
            C2_first_zero = k;
            break
        end
    end
    
    C2_first_zero,
    
    %     if C2_nr_rows == 1
    %         C2_max = X_array(C2_first_zero);
    %         C2_tol = Xstep;
    %         C2_quit = 'true';
    %         disp('Only 1 row remaining reached')
    if Xstep < Tolerance
        C2_max = X_array(C2_first_zero);
        C2_quit = 'true';
        disp('Tolerance reached')
    end
    
    % PREPARE NEXT ITERATION
    
    % reducing search area LOGIC:                                                   see TEST for first operation
    % if there is a zero in a column of C1_matrix, turn C1_mtag to 0
    % After scanning through C1_matrix, C1_mtag will contain ones where no
    % rows have a 0 and 0 where at least 1 row contains a 0 in that column.
    % The new X0 is chosen at the index before the first 0 in C1_mtag, and
    % X1 at the index after the first 0. Those new values will define a new
    % Xstep with again 10 steps between 0 and X1.
    
    if C2_first_zero < 0
        disp("Lower limit out of bounds, try another starting value X0")
        C2_quit = 'true';
    elseif C2_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    elseif C2_first_zero==Xlen
        X0 = X_array(C2_first_zero-1);
        X1 = X_array(C2_first_zero);
    else
        X0 = X_array(C2_first_zero-1);
        X1 = X_array(C2_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end

disp('Upper limit for C2')
disp(C2_max)
disp('Tolerance C2')
disp(C2_tol)

%end c2





%start c3

X1 = X1_original
X0 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);

disp('starting upper limit search for C3:')
C3_quit = 'false';
C3_iterations =0;
while strcmp(C3_quit,'false')
    C3_iterations = C3_iterations +1,
    C3_matrix = ones(1,Xlen*Ylen^4); % preallocate matrix of certain dimension
    i1 = 0;
    for c4 = Y0:Ystep:Y1
        for c5 = Y0:Ystep:Y1
            for c1 = Y0:Ystep:Y1
                for c2 = Y0:Ystep:Y1
                    for c3 = X0:Xstep:X1
                        i1 = i1+1;
                        cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
                        cross = cross - cross_threshold;
                        if cross >= 0
                            cross = 1;% outside of allowed region
                        else
                            cross = 0;% inside of allowed region
                        end
                        C3_matrix(1,i1) = cross;% row vector containing cross section info in binary, to be transformed into matrix (2D data)
                        
                    end
                end
            end
        end
    end
    C3_matrix = vec2mat(C3_matrix,Xlen); % make matrix from the row vector
    
    
    C3_rows = 0; % keep rows from C1_matrix which contain a zero
    for row = 1:size(C3_matrix,1)
        for column = 1:size(C3_matrix,2)
            if C3_matrix(row,column) == 0
                C3_rows = [C3_rows row]; %#ok<*AGROW>
                break
            end
        end
    end
    
    C3_rows = C3_rows(1,2:size(C3_rows,2)); % eliminate the first element which is artificial
    C3_nr_rows = size(C3_rows,2), % count of rows that contain a 0
    
    %TEST CURRENT ITERATION
    
    X_array = X0:Xstep:X1,
    
    C3_mtag = double(not(any(C3_matrix==0,1)));%returns logical, 0 = false, 1 = true
    
    C3_first_zero = -1;
    for k = 1:Xlen
        if C3_mtag(k)==1
            C3_first_zero = k;
            break
        end
    end
    
    C3_first_zero,
    
    %     if C3_nr_rows == 1
    %         C3_max = X_array(C3_first_zero);
    %         C3_tol = Xstep;
    %         C3_quit = 'true';
    %         disp('Only 1 row remaining reached')
    if Xstep < Tolerance
        C3_max = X_array(C3_first_zero);
        C3_quit = 'true';
        disp('Tolerance reached')
    end
    
    % PREPARE NEXT ITERATION
    
    % reducing search area LOGIC:                                                   see TEST for first operation
    % if there is a zero in a column of C1_matrix, turn C1_mtag to 0
    % After scanning through C1_matrix, C1_mtag will contain ones where no
    % rows have a 0 and 0 where at least 1 row contains a 0 in that column.
    % The new X0 is chosen at the index before the first 0 in C1_mtag, and
    % X1 at the index after the first 0. Those new values will define a new
    % Xstep with again 10 steps between 0 and X1.
    
    if C3_first_zero < 0
        disp("Upper limit out of bounds, try another starting value X0")
        C3_quit = 'true';
    elseif C3_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    elseif C3_first_zero==Xlen
        X0 = X_array(C3_first_zero-1);
        X1 = X_array(C3_first_zero);
    else
        X0 = X_array(C3_first_zero-1);
        X1 = X_array(C3_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end

disp('Upper limit for C3')
disp(C3_max)
disp('Tolerance C3')
disp(C3_tol)


%end c3




%start c4


X1 = X1_original;
X0 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);

disp('starting upper limit search for C4:')
C4_quit = 'false';
C4_iterations = 0;
while strcmp(C4_quit,'false')
    C4_iterations = C4_iterations +1,
    C4_matrix = ones(1,Xlen*Ylen^4); % preallocate matrix of certain dimension
    i1 = 0;
    for c5 = Y0:Ystep:Y1
        for c1 = Y0:Ystep:Y1
            for c2 = Y0:Ystep:Y1
                for c3 = Y0:Ystep:Y1
                    for c4 = X0:Xstep:X1
                        i1 = i1+1;
                        cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
                        cross = cross - cross_threshold;
                        if cross >= 0
                            cross = 1;% outside of allowed region
                        else
                            cross = 0;% inside of allowed region
                        end
                        C4_matrix(1,i1) = cross;% row vector containing cross section info in binary, to be transformed into matrix (2D data)
                        
                    end
                end
            end
        end
    end
    C4_matrix = vec2mat(C4_matrix,Xlen); % make matrix from the row vector
    
    
    C4_rows = 0; % keep rows from C1_matrix which contain a zero
    for row = 1:size(C4_matrix,1)
        for column = 1:size(C4_matrix,2)
            if C4_matrix(row,column) == 0
                C4_rows = [C4_rows row]; %#ok<*AGROW>
                break
            end
        end
    end
    
    C4_rows = C4_rows(1,2:size(C4_rows,2)); % eliminate the first element which is artificial
    C4_nr_rows = size(C4_rows,2), % count of rows that contain a 0
    
    %TEST CURRENT ITERATION
    
    X_array = X0:Xstep:X1,
    
    C4_mtag = double(not(any(C4_matrix==0,1)));%returns logical, 0 = false, 1 = true
    
    C4_first_zero = -1;
    for k = 1:Xlen
        if C4_mtag(k)==1
            C4_first_zero = k;
            break
        end
    end
    
    C4_first_zero,
    
    %     if C4_nr_rows == 1
    %         C4_max = X_array(C4_first_zero);
    %         C4_tol = Xstep;
    %         C4_quit = 'true';
    %         disp('Only 1 row remaining reached')
    if Xstep < Tolerance
        C4_max = X_array(C4_first_zero);
        C4_quit = 'true';
        disp('Tolerance reached')
    end
    
    % PREPARE NEXT ITERATION
    
    % reducing search area LOGIC:                                                   see TEST for first operation
    % if there is a zero in a column of C1_matrix, turn C1_mtag to 0
    % After scanning through C1_matrix, C1_mtag will contain ones where no
    % rows have a 0 and 0 where at least 1 row contains a 0 in that column.
    % The new X0 is chosen at the index before the first 0 in C1_mtag, and
    % X1 at the index after the first 0. Those new values will define a new
    % Xstep with again 10 steps between 0 and X1.
    
    if C4_first_zero < 0
        disp("Upper limit out of bounds, try another starting value X0")
        C4_quit = 'true';
    elseif C4_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    elseif C4_first_zero==Xlen
        X0 = X_array(C4_first_zero-1);
        X1 = X_array(C4_first_zero);
    else
        X0 = X_array(C4_first_zero-1);
        X1 = X_array(C4_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end

disp('Upper limit for C4')
disp(C4_max)
disp('Tolerance C4')
disp(C4_tol)


%end c4




%start c5


X1 = X1_original;
X0 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);

disp('starting upper limit search for C5:')
C5_quit = 'false';
C5_iterations = 0;
while strcmp(C5_quit,'false')
    C5_iterations = C5_iterations +1,
    C5_matrix = ones(1,Xlen*Ylen^4); % preallocate matrix of certain dimension
    i1 = 0;
    for c1 = Y0:Ystep:Y1
        for c2 = Y0:Ystep:Y1
            for c3 = Y0:Ystep:Y1
                for c4 = Y0:Ystep:Y1
                    for c5 = X0:Xstep:X1
                        i1 = i1+1;
                        cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
                        cross = cross - cross_threshold;
                        if cross >= 0
                            cross = 1;% outside of allowed region
                        else
                            cross = 0;% inside of allowed region
                        end
                        C5_matrix(1,i1) = cross;% row vector containing cross section info in binary, to be transformed into matrix (2D data)
                        
                    end
                end
            end
        end
    end
    C5_matrix = vec2mat(C5_matrix,Xlen); % make matrix from the row vector
    
    
    C5_rows = 0; % keep rows from C1_matrix which contain a zero
    for row = 1:size(C5_matrix,1)
        for column = 1:size(C5_matrix,2)
            if C5_matrix(row,column) == 0
                C5_rows = [C5_rows row]; %#ok<*AGROW>
                break
            end
        end
    end
    
    C5_rows = C5_rows(1,2:size(C5_rows,2)); % eliminate the first element which is artificial
    C5_nr_rows = size(C5_rows,2), % count of rows that contain a 0
    
    %TEST CURRENT ITERATION
    
    X_array = X0:Xstep:X1,
    
    C5_mtag = double(not(any(C5_matrix==0,1)));%returns logical, 0 = false, 1 = true
    
    C5_first_zero = -1;
    for k = 1:Xlen
        if C5_mtag(k)==1
            C5_first_zero = k;
            break
        end
    end
    
    C5_first_zero,
    
    %     if C5_nr_rows == 1
    %         C5_max = X_array(C5_first_zero);
    %         C5_tol = Xstep;
    %         C5_quit = 'true';
    %         disp('Only 1 row remaining reached')
    if Xstep < Tolerance
        C5_max = X_array(C5_first_zero);
        C5_quit = 'true';
        disp('Tolerance reached')
    end
    
    % PREPARE NEXT ITERATION
    
    % reducing search area LOGIC:                                                   see TEST for first operation
    % if there is a zero in a column of C1_matrix, turn C1_mtag to 0
    % After scanning through C1_matrix, C1_mtag will contain ones where no
    % rows have a 0 and 0 where at least 1 row contains a 0 in that column.
    % The new X0 is chosen at the index before the first 0 in C1_mtag, and
    % X1 at the index after the first 0. Those new values will define a new
    % Xstep with again 10 steps between 0 and X1.
    
    if C5_first_zero < 0
        disp("Upper limit out of bounds, try another starting value X0")
        C5_quit = 'true';
    elseif C5_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    elseif C5_first_zero==Xlen
        X0 = X_array(C5_first_zero-1);
        X1 = X_array(C5_first_zero);
    else
        X0 = X_array(C5_first_zero-1);
        X1 = X_array(C5_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end

disp('Upper limit for C5')
disp(C5_max)
disp('Tolerance C5')
disp(C5_tol)


%end c5


end
