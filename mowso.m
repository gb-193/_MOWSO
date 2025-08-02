function [ps,pf]=mowso(fname,n_var,n_obj,xl,xu,Max_Gen,popsize)
%% Initialization 
nRep=popsize;
VarSize=[1 n_var];
A.Position=[];
A.Velocity=[];
A.Cost=[];
A.Best.Position=[];
A.Best.Cost=[];
A.IsDominated=[];
A.GridIndex=[];%
A.GridSubIndex=[];%
pop=repmat(A,popsize,1);
%Initialize population information
for i=1:popsize
    pop(i).Position=NumFacility_mowso(n_var,xu,xl);
    pop(i).Velocity=zeros(VarSize);
    pop(i).Cost=feval(fname,pop(i).Position);
    pop(i).Best.Position=pop(i).Position;
    pop(i).Best.Cost=pop(i).Cost;
end
% MOWSO Parameters
fmax=0.75; %  Maximum frequency of the wavy motion
fmin=0.07; %  Minimum frequency of the wavy motion
tau=4.11;
mu=2/abs(2-tau-sqrt(tau^2-4*tau));
pmin=0.5;
pmax=1.5;
a0=6.250;
a1=100;
a2=0.0005;
for i=1:popsize
    ffmin(i,:)=pop(i).Cost;
    fitness(i,:)=pop(i).Cost;
end
[fmin00,index]=min(ffmin);
fmin0=ffmin(index(1),:);
for i=1:popsize
    wbest(i,:) = pop(i).Position; % Best position initialization
end
for i=1:n_obj
    gbest(i,:)= pop(index(i)).Position;
end
% Find Pareto optimal solutions
pop=DetermineDomination_mowso(pop);
rep=pop(~[pop.IsDominated]);
%% Start the iterative process of MOWSO
for it=1:Max_Gen
    mv=1/(a0+exp((Max_Gen/2.0-it)/a1));
    s_s=abs((1-exp(-a2*it/Max_Gen))) ;
    p1=pmax+(pmax-pmin)*exp(-(4*it/Max_Gen)^2);
    p2=pmin+(pmax-pmin)*exp(-(4*it/Max_Gen)^2);
    nu=floor((popsize).*rand(1,popsize))+1;
    rmin=1; rmax=3.0;
    rr=rmin+rand()*(rmax-rmin);
    wr=abs(((2*rand()) - (1*rand()+rand()))/rr);
    f =fmin+(fmax-fmin)/(fmax+fmin);
    for i=1:popsize
        a = sign( pop(i).Position-xu)>0;
        b = sign( pop(i).Position-xl)<0;
        wo=xor(a,b);
        leader_pos=leader(rep,n_obj);
        pop(i).Velocity=mu.*pop(i).Velocity +  wr.*(leader_pos-pop(i).Position);
        if rand<mv
            pop(i).Position=  pop(i).Position.*(~wo) + (xu.*a+xl.*b); % random allocation
        else
            pop(i).Position = pop(i).Position + pop(i).Velocity.*(1/f);  % based on the wavy motion
        end
        gbest_g=gbest(obj,:);
        for j=1:n_var
            if rand<s_s
                Dist=abs(rand*( gbest_g(j) - pop(i).Position(j) ));
                if(i==1)
                    pop(i).Position(j)=gbest_g(j)+rand*Dist*sign(rand-0.5);
                else
                    WSO_Pos(i,j)= gbest_g(j)+rand*Dist*sign(rand-0.5);
                    pop(i).Position(j)=(WSO_Pos(i,j)+pop(i-1).Position(j))/2*rand;
                end
            end
        end
    end
    for i=1:popsize
        pop(i).Position = constraint_mowso(pop(i).Position,xu,xl,n_var);
        pop(i).Cost=feval(fname,pop(i).Position);
        % Evaluate the fitness
        if (all(pop(i).Cost<=fitness(i,:)) & any(pop(i).Cost<fitness(i,:)))
            wbest(i,:) = pop(i).Position; % Update the best positions
            fitness(i,:)=pop(i).Cost;   % Update the fitness
        end
        if (all(fitness(i)<=fmin0) && any(fitness(i)<fmin0))
            fmin0=fitness(i);
            if rand>0.5
                gbest(1,:) = wbest(i,:); % Update the global best positions
            else
                gbest(2,:) = wbest(i,:); % Update the global best positions
            end
        end
        if Dominates_mowso(pop(i),pop(i).Best)
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
        elseif Dominates_mowso(pop(i).Best,pop(i))
            % Do Nothing
        else
            if rand<0.5
                pop(i).Best.Position=pop(i).Position;
                pop(i).Best.Cost=pop(i).Cost;
            end
        end
    end
    % Add Non-Dominated Particles
    pop=DetermineDomination_mowso(pop);
    rep=[rep
        pop(~[pop.IsDominated])];
    rep=DetermineDomination_mowso(rep);
    rep=rep(~[rep.IsDominated]);
    [n_ps,m_pf]=size(rep);
    for i=1:n_ps
        beiyong_ps(i,:)=rep(i).Position;
        beiyong_pf(i,:)=rep(i).Cost;
    end
    max_rep=numel(rep);
    if max_rep>nRep
        Extra=max_rep-nRep;
        for e=1:Extra
            d=sort_rep_mowso(rep,n_obj);
            rep=DeleteOneRepMemebr_mowso(rep,d);
        end
        for i=1:nRep
            ps(i,:)=rep(i).Position;
            pf(i,:)=rep(i).Cost;
        end
    end
end
if max_rep<=nRep
    ps=beiyong_ps;
    pf=beiyong_pf;
end
