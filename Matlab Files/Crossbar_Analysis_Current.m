function [I_out,g,VR,T,IM]  = Crossbar_Analysis_Current(V,VT,dt,g,VR,T,IM,param)


%Inputs--------------------------------------------------------------------
%V Input Voltage, vector [rows 1] dimensions
%VT Transistor Voltage, [1 col] dimensions
%dt timestep
%g gap (state variable) of each memristor, matrix of [rows col] dimensions
%VR Voltage across Output Resistances, [1 col] dimensions
%IM current of each memristor, matrix of [rows col] dimensions
%--------------------------------------------------------------------------

%Outputs-------------------------------------------------------------------
%I_out Output Current of each column, [1 col] dimensions
%g gap (state variable) of each memristor, matrix of [rows col] dimensions
%VR Voltage across Output Resistances, [1 col] dimensions
%IM current of each memristor, matrix of [rows col] dimensions
%--------------------------------------------------------------------------

[row,col]=size(g);
R=4*ones(1,col);
RT=1E8*((5-VT)/5);
% VT=ones(1,col);
% RT(VT>=2.5)=1;
% RT(VT<2.5)=1E8;
RT=repmat(RT,row,1);

%Memristor Parameter-------------------------------------------------------
kb= param(1);
q= param(2);
a0= param(3);
Eag= param(4);
Ear= param(5);
eL= param(6);
gap_min= param(7);
gap_max= param(8);
I0= param(9);
g0= param(10);
V0= param(11);
Vel0= param(12);
gamma0= param(13);
g1= param(14);
beta= param(15);
Cth= param(16);
Tau_th= param(17);
T0=param(18);
% kb= 1.38e-23;
% q= 1.6e-19;
% a0= 0.25e-9;
% Eag= 1.501;
% Ear= 1.5;
% eL= 5e-9;
% gap_min= 0.1e-9;
% gap_max= 1.7e-9;
% I0= 6.14e-5;
% g0= 2.75e-10;
% V0= 0.43;
% Vel0= 150;
% gamma0= 16.5;
% g1= 1e-9;
% beta= 1.25;
% Cth= 3.18e-16;
% Tau_th= 2.3e-10;
% T0=298;
%--------------------------------------------------------------------------

%Memristor Voltage Drop----------------------------------------------------
VM=repmat(V',1,col)-repmat(VR,row,1);
VM(RT>0.5E7)=0.0001;
% VM=repmat(V',1,col)-repmat(VR,row,1)-IM.*RT;
%--------------------------------------------------------------------------

%Memristor Stanford Meta Model--------------------------------------------- 
    IM=I0*exp(-g/g0).*sinh(VM/V0); 
    gamma=gamma0-beta*(g/g1).*(g/g1).*(g/g1);
    dg = -Vel0 * (exp(-(q*Eag)./(kb*T)).* exp((q*a0*gamma.*VM) ./ (eL*kb*T)) - exp(-(q*Ear) ./ (kb*T)) .* exp(-(q*a0*gamma.*VM) ./ (eL*kb*T)))*dt;	
    g=g+dg;	
    T = (T + dt .* (abs(VM.*IM) ./ Cth + T0 ./ Tau_th)) ./ (1 + dt/ Tau_th);
    g(g>gap_max)=gap_max;
    g(g<gap_min)=gap_min;
%--------------------------------------------------------------------------

%Crossbar Analysis---------------------------------------------------------
        I_out= sum(IM,1);
        VR=sum(IM,1).*R;
%--------------------------------------------------------------------------    
end
