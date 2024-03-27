function [uLn,vLn,wLn,gLn,uLj,vLj,wLj,gLj,gr,grj] = Big_mode_tt_kron_RK3_v2_module(uLn,vLn,wLn,gLn,uLj,vLj,wLj,gLj,umn,vmn,wmn,gmn,NT,h,Lkn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global kLyap NLyap N NX MZ dt L igpu


%[uLn,vLn,wLn,gLn]=compact_Lyap(uLn,vLn,wLn,gLn);
%[uLj,vLj,wLj,gLj]=compact_Lyap(uLj,vLj,wLj,gLj);

umean=squeeze(umn(:,1,:));
vmean=squeeze(vmn(:,1,:));
wmean=squeeze(wmn(:,1,:));
gmean=squeeze(gmn(:,1,:));


[un,vn,wn,gn] = Eigen_tt_module_kron_RK3_v2_k(uLn,vLn,wLn,gLn,umean,vmean,wmean,gmean,NT,h,Lkn(1:2:end));

[uj,vj,wj,gj] = Adjoi_tt_module_kron_RK3_v2_k(uLj,vLj,wLj,gLj,umean,vmean,wmean,gmean,NT,h,Lkn(1:2:end));

% if igpu
% [un,vn,wn,gn] = gather(un,vn,wn,gn);
% [uj,vj,wj,gj] = gather(uj,vj,wj,gj);
% [uLn,vLn,wLn,gLn] = gather(uLn,vLn,wLn,gLn);
% [uLj,vLj,wLj,gLj] = gather(uLj,vLj,wLj,gLj);
% end

[uLn,vLn,wLn,gLn,uLj,vLj,wLj,gLj,gr,grj]=biorth_adj_tt(un,vn,wn,gn,uLn,vLn,wLn,gLn,uj,vj,wj,gj,uLj,vLj,wLj,gLj,NT,h,Lkn);


% [uLn,vLn,wLn,gLn]=decompact_Lyap(uLn,vLn,wLn,gLn);
% [uLj,vLj,wLj,gLj]=decompact_Lyap(uLj,vLj,wLj,gLj);

end

