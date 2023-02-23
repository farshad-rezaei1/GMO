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
% Main paper: Rezaei, F., Safavi, H.R., Elaziz, M.A., Mirjalili, S. %
% Geometric Mean Optimizer for Solving Engineering Problems. Soft   %
% Computing (2023).                                                 %
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

