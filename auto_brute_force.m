
load('analytical2.mat')
% 
% i = 0;
% 
% cmin = -4*pi;
% cstep = 1;
% cmax = 4*pi;
% cstep2 = 0.8;
% 
% 
% s = size(cmin:cstep:cmax,2);
% s2 = size(-4*pi:cstep2:4*pi,2);
% 
% matrix1 = ones(1,s*s2^4);
% c_matrix1 = ones(5,s*s2^4);
% 
% 
% for c4 = -4*pi:cstep2:4*pi
%     for c3 = -4*pi:cstep2:4*pi
%         for c2 = -4*pi:cstep2:4*pi
%             for c5 = -4*pi:cstep2:4*pi
%                 for c1 = cmin:cstep:cmax
%                     i = i+1;
%                     cross = sig_i * [c1 ; c2 ; c3 ; c4 ; c5] + [c1 c2 c3 c4 c5]*sig_ij*[c1 ; c2 ; c3 ; c4 ; c5];
%                     cross = cross - sigma_exp_SM;
%                     if cross >= 0
%                         cross = 1;
%                     else
%                         cross = 0;
%                     end
%                     matrix1(1,i) = cross;
%                     c_matrix1(1,i) = c1;
%                     c_matrix1(2,i) = c2;
%                     c_matrix1(3,i) = c3;
%                     c_matrix1(4,i) = c4;
%                     c_matrix1(5,i) = c5;
%                    
%                 end
%             end
%         end
%     end
% end
% 
% matrix1 = vec2mat(matrix1,s);
% 
% rows1 = 0;
% for row = 1:size(matrix1,1)
%     for column = 1:size(matrix1,2)
%         if matrix1(row,column) == 0
%             rows1 = [rows1 row]; %#ok<*AGROW>
%             break
%         end
%     end
% end
% 
% rows1 = rows1(1,2:size(rows1,2));
% nr_rows1 = size(rows1,2),
% 
% 
% scr_info = get(groot,'ScreenSize');
% figure(1)
% fig1 = figure(1);
% set(fig1, 'Position', [10 (scr_info(4)/2-100) scr_info(3)/2 scr_info(4)/2]);
% 
% 
% clr = [1 1 1; 0 0 0];
% 
% imagesc(matrix1(rows1,:))
% colormap(clr)


c1_min = -3.3918;
% cmin = -1.08*pi; %%%% nr 12
% cstep = 0.0001;
% cmax = -1.079*pi;
% cstep2 = 0.8;

c1_max = 3.3036;
% cmin = 1.0508*pi; %%%%nr 25
% cstep = 0.0001;
% cmax = 1.0517*pi;
% cstep2 = 0.8;

c2_min = -7.2678;
% cmin = -2.315*pi; %%%% nr 11
% cstep = 0.0005;
% cmax = -2.312*pi;
% cstep2 = 0.8;

c2_max = 7.2838;
% cmin = 2.31820*pi; %%%%% nr 21
% cstep = 0.00005;
% cmax = 2.31875*pi;
% cstep2 = 0.8;

c3_min = -45.3223;
% cmin = -14.4268*pi; %%%%% nr 9
% cstep = 0.0001;
% cmax = -14.426*pi;
% cstep2 = 0.8;

c3_max = 45.1946;
% cmin = 14.3858*pi; %%%% nr 7
% cstep = 0.00005;
% cmax = 14.3860*pi;
% cstep2 = 0.8;

c4_min = -5.9609;
% cmin = -1.89750*pi; %%%% nr 6
% cstep = 0.00005;
% cmax = -1.89735*pi;
% cstep2 = 0.8;


c4_max = 5.8129;
% cmin = 1.8502*pi; %%%%% nr 7
% cstep = 0.00005;
% cmax = 1.8504*pi;
% cstep2 = 0.8;

c5_min = -12.8190;
% cmin = -4.082*pi;
% cstep = 0.0005;
% cmax = -4.078*pi;
% cstep2 = 0.8;

c5_max = 11.6189;
% cmin = 3.6980*pi;
% cstep = 0.0001;
% cmax = 3.6990*pi;
% cstep2 = 0.8;


c_boundaries = [c1_min 0.0001 c1_max 0.0001; c2_min 0.0005 c2_max 0.00005; c3_min 0.0001 c3_max 0.00005; c4_min 0.00005 c4_max 0.00005; c5_min 0.0005 c5_max 0.0001]
%these should be < abs(4pi)
%both O^8 (strong coupling) : don't have exp senitivity yet (maybe with more Lum.)


[t1,~,~,~,~,t2,~,~,~,~]=get_lower_limit(sigma_exp_SM_2016,-5,10^-4,sig_i,sig_ij)


function [C1_min,C2_min,C3_min,C4_min,C5_min,C1_tol,C2_tol,C3_tol,C4_tol,C5_tol] = get_lower_limit(cross_threshold,X0,Tolerance,sig_i,sig_ij)% C = investigated c
% scan over other parameters than the one investigated
Y0 = -4*pi;
Y1 = 4*pi;
Ystep = 0.8;
Ylen = size(Y0:Ystep:Y1,2); % length array

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
X1 = 0;%only lower bound
Xstep = abs(X1-X0)/10;
Xlen = size(X0:Xstep:X1,2);


%start c1

C1_quit = 'false';

while strcmp(C1_quit,'false')

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
    
    C1_mtag = double(any(C1_matrix==0,1));%returns logical, 0 = false, 1 = true
    
    C1_first_zero = -1;
    for k = 1:Xlen
        if C1_mtag(k)==1
            C1_first_zero = k;
            break
        end
    end
    
    C1_first_zero,
    
    if C1_nr_rows == 1
        C1_min = X_array(C1_first_zero);
        C1_tol = Xstep;
        C1_quit = 'true';
    elseif Xstep < Tolerance
        C1_min = X_array(C1_first_zero);
        C1_quit = 'true';
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
        C1_err = "Lower limit out of bounds, try another starting value X0",
        C1_quit = 'true';
    elseif C1_first_zero==1  % if the the first element is zero then go 1 step back
        X0 = X0-Xstep;
        X1 = X0+Xstep;
    else
        X0 = X_array(C1_first_zero-1);
        X1 = X_array(C1_first_zero+1);
    end
    Xstep = abs(X1-X0)/10,
    
end
%outside while-loop C1
%matrix1(C1_rows,:)



%end c1


%start c2


end



