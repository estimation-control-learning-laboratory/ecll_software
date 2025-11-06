function [ G ] = get_TF_from_ABCD( A,B,C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

lu = size(B,2);
ly = size(C,1);


G = tf(zeros(ly,lu));
for ii = 1:lu
    for jj = 1:ly
        [NUM,DEN] = ss2tf(A,B(:,ii),C(jj,:),0);
        G(jj, ii) = tf(NUM, DEN, 1);
    end
end
        
% 
% if isfield(FLAG,'plot_MIMO_Gc')
%     if FLAG.plot_MIMO_Gc
%         counter = 0;
%         for ii = 1:lu
%             for jj = 1:ly
%                 counter = counter+1;
%                 subplot(lu, ly, counter)
%                 pzmap(G(ii,jj))
%             end
%         end
%     end
% end    


end

