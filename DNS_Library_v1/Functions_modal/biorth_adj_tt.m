function [uLn,vLn,wLn,gLn,uLj,vLj,wLj,gLj,gr,grj] = biorth_adj(un,vn,wn,gn,uLn,vLn,wLn,gLn,uj,vj,wj,gj,uLj,vLj,wLj,gLj,NT,h,Lkn)

global kLyap dt

%%%% Biorthogonalization 

    % Initial norms 
    
    Eni = Ener_kn(uLn,vLn,wLn);
    Eji = Ener_kn(uLj,vLj,wLj); 

% split fields before entering function?
Lkiter=0;
for Lk=kLyap
    Lkiter=Lkiter+1;
    
    Lnn = find(Lkn == Lk);
    
    % First mode normalization
    
    for Ln = Lnn
    
        if rem(Ln,2) == 1
    
        [uLt,vLt,wLt,gLt]=pick_vec_kn(un,vn,wn,gn,Ln);% Final basis moda
        [uJt,vJt,wJt,gJt]=pick_vec_kn(uj,vj,wj,gj,Ln); % Final basis adjoint

        Enf=Ener_k(uLt,vLt,wLt);
        Ejf=Ener_k(uJt,vJt,wJt);

        gr(Ln)=log(Enf/Eni(Ln))/(2*h*(NT-1)*dt);
        grj(Ln)=log(Ejf/Eji(Ln))/(2*h*(NT-1)*dt);
       
        
        %%% Normalize and shift adjoint to unit mode projection
        [uLn(:,:,Ln),vLn(:,:,Ln),wLn(:,:,Ln),gLn(:,:,Ln)]=Lyap_norm_kn(uLt,vLt,wLt,gLt,1,Enf);
        
        [uL0,vL0,wL0,gL0]=pick_vec_kn(uLn,vLn,wLn,gLn,Ln);

        P1 = project_adj_k(uJt,vJt,wJt,uL0,vL0,wL0);
        [uLj(:,:,Ln),vLj(:,:,Ln),wLj(:,:,Ln),gLj(:,:,Ln)]=adj_norm_k(uJt,vJt,wJt,gJt,1/P1');
        
        [uJ0,vJ0,wJ0,gJ0]=pick_vec_kn(uLj,vLj,wLj,gLj,Ln);

        
        % project out of all base vectors

        if Ln ~= Lnn(end)
            
                
        [unn,vnn,wnn,gnn]=pick_vec_kn(un,vn,wn,gn,Ln+1:Lnn(end));
        %[uL0,vL0,wL0,gL0]=pick_vec(uL,vL,wL,gL,Ln,Lk)
        [ujn,vjn,wjn,gjn]=pick_vec_kn(uj,vj,wj,gj,Ln+1:Lnn(end));
            

        P1n=project_adj_kn(uJ0,vJ0,wJ0,unn,vnn,wnn);
                
        [un(:,:,Ln+1:Lnn(end)),vn(:,:,Ln+1:Lnn(end)),wn(:,:,Ln+1:Lnn(end)),gn(:,:,Ln+1:Lnn(end))] = project_out(uL0,vL0,wL0,gL0,unn,vnn,wnn,gnn,P1n);
                               
        P2n=project_adj_kn(uL0,vL0,wL0,ujn,vjn,wjn);
              
        [uj(:,:,Ln+1:Lnn(end)),vj(:,:,Ln+1:Lnn(end)),wj(:,:,Ln+1:Lnn(end)),gj(:,:,Ln+1:Lnn(end))] = project_out(uJ0,vJ0,wJ0,gJ0,ujn,vjn,wjn,gjn,P2n);  
        
        else      
        
        end
    
        else 
            
        end
        
    end
    


end

    %%% Symmetry translation
    gr(2:2:Ln) = gr(1:2:Ln); %gr = permute(gr,[3 1 2]);
    grj(2:2:Ln) = grj(1:2:Ln); %cgr = permute(cgr,[3 1 2]);
    
    [uLn(:,:,2:2:end),vLn(:,:,2:2:end),wLn(:,:,2:2:end),gLn(:,:,2:2:end)]=Lyap_translate_k(uLn(:,:,1:2:end),vLn(:,:,1:2:end),wLn(:,:,1:2:end),gLn(:,:,1:2:end));

    [uLj(:,:,2:2:end),vLj(:,:,2:2:end),wLj(:,:,2:2:end),gLj(:,:,2:2:end)]=Lyap_translate_k(uLj(:,:,1:2:end),vLj(:,:,1:2:end),wLj(:,:,1:2:end),gLj(:,:,1:2:end));
    