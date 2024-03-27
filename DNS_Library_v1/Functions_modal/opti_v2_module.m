function [uLn,vLn,wLn,gLn,gr] = opti_v2_module(uLn,vLn,wLn,gLn,umn,vmn,wmn,gmn,NT,h,Lkn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global kLyap NLyap N NX MZ dt L


%[uLn,vLn,wLn,gLn]=compact_Lyap(uLn,vLn,wLn,gLn);
%[uLj,vLj,wLj,gLj]=compact_Lyap(uLj,vLj,wLj,gLj);

umean=squeeze(umn(:,1,:));
vmean=squeeze(vmn(:,1,:));
wmean=squeeze(wmn(:,1,:));
gmean=squeeze(gmn(:,1,:));

[un,vn,wn,gn] = Eigen_tt_module_kron_RK3_v2_k(uLn,vLn,wLn,gLn,umean,vmean,wmean,gmean,NT,h,Lkn(1:2:end));

[uL,vL,wL,gL] = Adjoi_tt_module_kron_RK3_v2_k(un,vn,wn,gn,umean,vmean,wmean,gmean,NT,h,Lkn(1:2:end));

    % Initial norms 
    
    Eni = Ener_kn(uLn,vLn,wLn);

Lkiter=0;
for Lk=kLyap
    Lkiter=Lkiter+1;

    Lnn = find(Lkn == Lk);

    for Ln= Lnn

        if rem(Ln,2) == 1
        
%        [uLi,vLi,wLi,gLi]=pick_vec_kn(uLn,vLn,wLn,gLn,Ln);
        [uLt,vLt,wLt,gLt]=pick_vec_kn(un,vn,wn,gn,Ln);
        [uJt,vJt,wJt,gJt]=pick_vec_kn(uL,vL,wL,gL,Ln);


%        if Ln==1

            Ei=Ener_k(uJt,vJt,wJt);
            Ef=Ener_k(uLt,vLt,wLt);

%            gr(Ln)=log(Ef/Eni(Ln))/(2*h*dt);
            gr(Ln)=Ef/Eni(Ln);

            [uLn(:,:,Ln),vLn(:,:,Ln),wLn(:,:,Ln),gLn(:,:,Ln)]=Lyap_norm_kn(uJt,vJt,wJt,gJt,1,Ei);

%       else
 
            [uL0,vL0,wL0,gL0]=pick_vec_kn(uLn,vLn,wLn,gLn,Ln);
            if Ln ~= Lnn(end)
                [uLb,vLb,wLb,gLb]=pick_vec_kn(uL,vL,wL,gL,Ln+1:Lnn(end));
                %[uL0,vL0,wL0,gL0]=pick_vec(uL,vL,wL,gL,Ln,Lk)

                P1=project_adj_kn(uL0,vL0,wL0,uLb,vLb,wLb);

                [uL(:,:,Ln+1:Lnn(end)),vL(:,:,Ln+1:Lnn(end)),wL(:,:,Ln+1:Lnn(end)),gL(:,:,Ln+1:Lnn(end))] = project_out(uL0,vL0,wL0,gL0,uLb,vLb,wLb,gLb,P1);

            end
%                 Ei=Ener_k(uLi,vLi,wLi);
%                 Ef=Ener_k(uLt,vLt,wLt);
% 
%                 gr(Ln,Lk)=log(Ef/Ei)/(2*dt);
% 
%                 [uLn(:,:,Ln,Lk),vLn(:,:,Ln,Lk),wLn(:,:,Ln,Lk),gLn(:,:,Ln,Lk)]=Lyap_norm_k(uLt,vLt,wLt,gLt,1);

        end

%    end
    end

    gr(2:2:Ln) = gr(1:2:Ln); %gr = permute(gr,[3 1 2]);

    [uLn(:,:,2:2:end),vLn(:,:,2:2:end),wLn(:,:,2:2:end),gLn(:,:,2:2:end)]=Lyap_translate_k(uLn(:,:,1:2:end),vLn(:,:,1:2:end),wLn(:,:,1:2:end),gLn(:,:,1:2:end));



% [uLn,vLn,wLn,gLn]=decompact_Lyap(uLn,vLn,wLn,gLn);
% [uLj,vLj,wLj,gLj]=decompact_Lyap(uLj,vLj,wLj,gLj);

end

