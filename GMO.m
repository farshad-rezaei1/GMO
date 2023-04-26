%___________________________________________________________________%
% GMO: Geometric Mean Optimizer                                     %
%                                                                   %
% Developed in MATLAB R2018b                                        %
%                                                                   %
% Inventor and programmer: Farshad Rezaei, PhD                      %
%                                                                   %
% e-Mail: farshad.rezaei@gmail.com                                  %
%         f.rezaei@alumni.iut.ac.ir                                 %
%                                                                   %
% Homepage: https://www.linkedin.com/in/farshad-rezaei-5a92559a/    %
%                                                                   %
% Main paper: Rezaei, F., Safavi, H.R., Abd Elaziz, M. et al. GMO:  %
% geometric mean optimizer for solving engineering problems.        %
% Soft Comput (2023). https://doi.org/10.1007/s00500-023-08202-z    %
%___________________________________________________________________%

% GMO: Geometric Mean Optimizer                                                                 
function [z_iter,z_final,pos_final] = GMO(np,nx,maxit,varmax,varmin,velmax,velmin,epsilon,fobj)
% disp(['Number of Iterations = ',num2str(it)]);
pp_pbest=zeros(np,nx);
pp_kbest=zeros(np,nx);
pp_guide=zeros(np,nx);
mutant=zeros(1,nx);
stdev2=zeros(nx);
index=zeros(np);
fit=zeros(maxit,np);
DFI=zeros(np);
optimal_pos=zeros(1,nx);
z_pbest=zeros(np);
z_kbest=zeros(np);
z_optimal=inf*ones(maxit);
pos_final=zeros(nx);
z_iter=zeros(maxit);
kbest_max=np;
kbest_min=2;

% Starting the optimization process
it=1;

% Initialization process of the algorithm
[pp,pv]=Initialization(np,nx,varmax,varmin,velmax,velmin);

% Objective function evaluations and determining the best-so-far solutions and objectives
for j=1:np
    z=fobj(pp(j,1:nx));
    z_pbest(j)=z;
    pp_pbest(j,1:nx)=pp(j,1:nx);
end
for i=1:nx
    stdev2(1,i)=std(pp_pbest(1:np,i));
end
max_stdev2=max(stdev2(1,1:nx));

% Calculating Mean and Std of the best-so-far solutions' objective values
ave=mean(z_pbest(1:np));
stdev=std(z_pbest(1:np));

% Determining the number of the Elite best-so-far solutions
kbest=kbest_max-(kbest_max-kbest_min)*(it/maxit);
n_best=round(kbest);

% Evaluating the Dual-Fitness Indices of the solutions- Eq.(3)
for j=1:np
    index(j) = j;
end
for j=1:np
    prod=1;
    for jj=1:np
        if jj~=j
            prod=prod*(1/(1+exp((-4)/(stdev*sqrt(exp(1)))*(z_pbest(jj)-ave)))); 
        end
    end
    fit(it,j)=prod;
end
for j=1:np-1
    for jj=j+1:np
        if fit(it,jj)>fit(it,j)
            c1=fit(it,j);
            fit(it,j)=fit(it,jj);
            fit(it,jj)=c1;
            c2=index(j);
            index(j)=index(jj);
            index(jj)=c2;
        end
    end
end

% Designating the Elite solutions
sum1=0;
for j=1:n_best
    z_kbest(j)=z_pbest(index(j));
    pp_kbest(j,1:nx)=pp_pbest(index(j),1:nx);
    sum1=sum1+fit(it,j);
end 

% Calculating the guide solutions
for j=1:np
    pp_guide(j,1:nx)=zeros(1,nx);
    for jj=1:n_best
        if index(jj)~=j
            DFI(jj)=fit(it,jj)/(sum1+epsilon);
            pp_guide(j,1:nx)=pp_guide(j,1:nx)+DFI(jj).*pp_kbest(jj,1:nx); % Eq.(5)
        end
    end
end

% Determining the best objective value and solution
for j=1:np
    if z_pbest(j)<z_optimal(it)
        z_optimal(it)=z_pbest(j);
        optimal_pos(it,1:nx)=pp_pbest(j,1:nx);
    end
end

% Saving the best-so-far objective value in the current run
z_iter(it)=z_optimal(it);

% The Main Loop
while it<maxit
    it=it+1;        
    w=1-(it/maxit); % Eq.(9)
%     disp(['Number of Iterations = ',num2str(it)]);
    for j=1:np 
        % Mutating the guide solutions- Eq.(6)
        mutant(1,1:nx)=pp_guide(j,1:nx)+w*randn(1,nx).*(max_stdev2-stdev2(1,1:nx));
        
        % Updating the velocity of the solutions- Eq.(7)
        pv(j,1:nx) = w*pv(j,1:nx)+(ones(1,nx)+(2*rand(1,nx)-ones(1,nx))*w).*(mutant(1,1:nx)-pp(j,1:nx));
        
        % Returning back the velocity of the solutions if going beyond the velocity boundaries
        flag4lbv=pv(j,:)<velmin(1,:);
        flag4ubv=pv(j,:)>velmax(1,:);
        pv(j,:)=(pv(j,:)).*(~(flag4lbv+flag4ubv))+velmin.*flag4lbv+velmax.*flag4ubv;
        
        % Updating the position of the solutions- Eq.(8)
        pp(j,:)=pp(j,:)+pv(j,:);
        
        % Returning back the position and velocity of the solutions if going beyond the position boundaries
        flag4lbp=pp(j,:)<varmin(1,:);
        flag4ubp=pp(j,:)>varmax(1,:);
        pp(j,:)=(pp(j,:)).*(~(flag4lbp+flag4ubp))+varmin.*flag4lbp+varmax.*flag4ubp; 
        pv(j,:)=(pv(j,:)).*(ones(1,nx)-2*(flag4lbp+flag4ubp)); 
    
        % Objective function evaluations and determining the best-so-far solutions and objectives
        z=fobj(pp(j,1:nx));
        if z<z_pbest(j)
            z_pbest(j)=z;
            pp_pbest(j,1:nx)=pp(j,1:nx);
        end
    end
    for i=1:nx
        stdev2(1,i)=std(pp_pbest(1:np,i));
    end
    max_stdev2=max(stdev2(1,1:nx));

    % Calculating Mean and Std of the best-so-far solutions' objective values
    ave=mean(z_pbest(1:np));
    stdev=std(z_pbest(1:np));

    % Determining the number of the Elite best-so-far solutions
    kbest=kbest_max-(kbest_max-kbest_min)*(it/maxit);
    n_best=round(kbest);

    % Evaluating the Dual-Fitness Indices of the solutions- Eq.(3)
    for j=1:np
        index(j) = j;
    end
    for j=1:np
        prod=1;
        for jj=1:np
            if jj~=j
                prod=prod*(1/(1+exp((-4)/(stdev*sqrt(exp(1)))*(z_pbest(jj)-ave))));
            end
        end
        fit(it,j)=prod;
    end
    for j=1:np-1
        for jj=j+1:np
            if fit(it,jj)>fit(it,j)
                c1=fit(it,j);
                fit(it,j)=fit(it,jj);
                fit(it,jj)=c1;
                c2=index(j);
                index(j)=index(jj);
                index(jj)=c2;
            end
        end
    end
    
    % Designating the Elite solutions
    sum1=0;
    for j=1:n_best
        z_kbest(j)=z_pbest(index(j));
        pp_kbest(j,1:nx)=pp_pbest(index(j),1:nx);
        sum1=sum1+fit(it,j);
    end
    
    % Calculating the guide solutions
    for j=1:np
        pp_guide(j,1:nx)=zeros(1,nx);
        for jj=1:n_best
            if index(jj)~=j
                DFI(jj)=fit(it,jj)/(sum1+epsilon);
                pp_guide(j,1:nx)=pp_guide(j,1:nx)+DFI(jj).*pp_kbest(jj,1:nx); % Eq.(5)
            end
        end
    end
    
    % Determining the best-so-far objective value and solution
    for j=1:np
        if z_pbest(j)<z_optimal(it)
            z_optimal(it)=z_pbest(j);
            optimal_pos(it,1:nx)=pp_pbest(j,1:nx);
        end
    end
    
    % Saving the best-so-far objective value in the current run
    z_iter(it)=z_optimal(it);
end

% Saving the final best solution and objective revealed upon the end of the optimization process
z_final=z_iter(maxit);
pos_final(1:nx)=optimal_pos(maxit,1:nx);
end