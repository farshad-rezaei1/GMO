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

% This function is to initialize the position and velocity of the solutions to start the optimization process
function [pp,pv] = Initialization(np,nx,varmax,varmin,velmax,velmin)
pp=zeros(np,nx);
pv=zeros(np,nx);

for j=1:np
    pp(j,1:nx)=(varmax-varmin).*rand(1,nx)+varmin;
    pv(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
end
end

