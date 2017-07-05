
load('analytical.mat')

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

