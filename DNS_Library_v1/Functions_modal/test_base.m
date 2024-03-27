function P11 = test_base(uLj,vLj,wLj,gLj,uLn,vLn,wLn,gLn)

global kLyap NLyap

for kcount1=1:kLyap
        for kcount2=kcount1
            for lcount1=1:NLyap
                for lcount2=1:NLyap
                    
                    [uL0,vL0,wL0,gL0]=pick_vec(uLj,vLj,wLj,gLj,lcount1,kcount1);
                    [uL1,vL1,wL1,gL1]=pick_vec(uLn,vLn,wLn,gLn,lcount2,kcount2);
                    P11(lcount1,lcount2,kcount1)=project_Lyap(uL0,vL0,wL0,uL1,vL1,wL1,kLyap(kcount1)+1);
                end
            end
        end
end
    
end