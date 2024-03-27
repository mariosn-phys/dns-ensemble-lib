function [uLn,vLn,wLn,gLn] = orth_Lyap(uL,vL,wL,gL)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global kLyap NLyap

Lkiter=0;
for Lk=kLyap
    Lkiter=Lkiter+1;
    
    for Ln=1:NLyap(Lkiter)
        
         [uLt,vLt,wLt,gLt]=pick_vec(uL,vL,wL,gL,Ln,Lkiter);
        
        if Ln==1
            
            [uLn(:,:,:,Ln,Lkiter),vLn(:,:,:,Ln,Lkiter),wLn(:,:,:,Ln,Lkiter),gLn(:,:,:,Ln,Lkiter)]=Lyap_norm(uLt,vLt,wLt,gLt,1);
            
        else
            [uL0,vL0,wL0,gL0]=pick_vec(uL,vL,wL,gL,Ln,Lkiter);
            
            for itr=1:Ln-1
                
                [uLb,vLb,wLb,gLb]=pick_vec(uLn,vLn,wLn,gLn,itr,Lkiter);
                %[uL0,vL0,wL0,gL0]=pick_vec(uL,vL,wL,gL,Ln,Lk)
                
                P1=project_Lyap(uLb,vLb,wLb,uL0,vL0,wL0,Lk+1);
                
                [uLt,vLt,wLt,gLt] = project_out(uLb,vLb,wLb,gLb,uLt,vLt,wLt,gLt,P1);
                
            end
            
            [uLn(:,:,:,Ln,Lkiter),vLn(:,:,:,Ln,Lkiter),wLn(:,:,:,Ln,Lkiter),gLn(:,:,:,Ln,Lkiter)]=Lyap_norm(uLt,vLt,wLt,gLt,1);
            
        end
        
    end
end


end

