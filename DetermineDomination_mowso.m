function pop=DetermineDomination_mowso(pop) %pop(j).IsDominated的值为真 表示j被支配了
    N=numel(pop);

    for i=1:N
        pop(i).IsDominated=false;
    end

    for i=1:N-1
        for j=i+1:N    
            if Dominates_mowso(pop(i),pop(j)) %i支配j
               pop(j).IsDominated=true;
            end
            if Dominates_mowso(pop(j),pop(i)) %j支配i
               pop(i).IsDominated=true;
            end     
        end
    end
end