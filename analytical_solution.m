clear;clc;

% c_i denotes the values for the 5 Wilson coefficients in the same order as
% on the thesis/ or FeynRules. These values were the input to calculate the
% MG5 cross sections. These 20 sets of c_i are used in the function gen_row
% to create a row with coefficients for the constants sigma_i,ij. the
% resulting matrix is rank 20. The coefficients are obtained in this order:

%order of sigmas : ['sig1';'sig2';'sig3';'sig4';'sig5';'sig11';'sig22';'sig33';'sig44';'sig55';'sig12';'sig13';'sig14';'sig15';'sig23';'sig24';'sig25';'sig34';'sig35';'sig45'];

% MG_SM = cross MG5 - cross SM (from MG5)


c1 = [1 0 0 0 0];%done
c2 = [0 1 0 0 0];%done
c3 = [0 0 1 0 0];%done
c4 = [0 0 0 1 0];%done
c5 = [0 0 0 0 1];%done
c6 = [1 1 1 1 1];%done
c7 = [-1 -1 1 1 1];%done
c8 = [-1 -1 1 0 1];%done
c9 = [0 1 0 0 -1];%done
c10 = [0 1 1 1 0];%done
c11 = [1 0 -1 1 0];%done
c12 = [-1 0 0 1 -1];%done
c13 = [-1 0 0 -1 1];%done
c14 = [0 1 -1 1 -1];%done
c15 = [0 1 0 -1 0];%done
c16 = [0 0 -1 -1 1];%done
c17 = [1 -1 0 -1 0];%done
c18 = [1 1 0 -1 1];%done
c19 = [0 1 0 -1 1];%done
c20 = [1 -1 -1 0 1];%done

sig_SM = 0.009308;
MG_SM = [0.01557, 0.01564, 0.0102, 0.01116, 0.01022, 0.02873, 0.0203, 0.01704, 0.01527, 0.02066, 0.0168, 0.01741, 0.01567, 0.01342, 0.01839, 0.01208, 0.02113, 0.0283, 0.01983, 0.02386]- sig_SM;
MG_SM = transpose(MG_SM);

A = [gen_row(c1) ; gen_row(c2) ; gen_row(c3) ; gen_row(c4) ; gen_row(c5) ; gen_row(c6) ; gen_row(c7) ; gen_row(c8) ; gen_row(c9) ; gen_row(c10) ; gen_row(c11) ; gen_row(c12) ; gen_row(c13) ; gen_row(c14) ; gen_row(c15) ; gen_row(c16) ; gen_row(c17) ; gen_row(c18) ;gen_row(c19) ;gen_row(c20)];
rank = rank(A);

S = inv(A)*MG_SM,%Solution

%order of sigmas : ['sig1';'sig2';'sig3';'sig4';'sig5';'sig11';'sig22';'sig33';'sig44';'sig55';'sig12';'sig13';'sig14';'sig15';'sig23';'sig24';'sig25';'sig34';'sig35';'sig45'];

sig_i = [S(1) S(2) S(3) S(4) S(5)];
sig_ij = [S(6) S(11) S(12) S(13) S(14); S(11) S(7) S(15) S(16) S(17); S(12) S(15) S(8) S(18) S(19); S(13) S(16) S(18) S(9) S(20); S(14) S(17) S(19) S(20) S(10)];

sigma_exp_SM = 6.4* 0.00930792902; % This was the previous limit of ~69fb

% just to print them to full accuracy
sprintf('%.16f',sigma_exp_SM)
%LIMIT PRINT

%SIGMA_I,IJ PRINT
for i = 1:20
    sprintf('%.16f',S(i)),
end



function row = gen_row(c)
    row = [c(1) c(2) c(3) c(4) c(5) c(1)^2 c(2)^2 c(3)^2 c(4)^2 c(5)^2 2*c(1)*c(2) 2*c(1)*c(3) 2*c(1)*c(4) 2*c(1)*c(5) 2*c(2)*c(3) 2*c(2)*c(4) 2*c(2)*c(5) 2*c(3)*c(4) 2*c(3)*c(5) 2*c(4)*c(5)];
end