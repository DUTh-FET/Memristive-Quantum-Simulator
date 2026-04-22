%INSERT COIN (OR GATE)
function Output_states  = Crossbar_Computing(g,Qubit_State,coefficients,param)

[row,col,gate_num]=size(g);
dt=1e-12;
%Computing Pulse Generation------------------------------------------------
%     Precision
%     Pulse_Amplitude=0.1;
%     Pulse_Width=param(20);
%     Pulse_duration=Pulse_Width+0.2*Pulse_Width;
%     n=0;
%     %Find the number of digits after the decimal point in a number
%     while (floor(Pulse_duration*10^n)~=Pulse_duration*10^n)
%         n=n+1;
%     end
%     p=Pulse_duration/10^(-n);
%     %Counting the number of digits
%     p=numel(num2str(p))-1;
%     dt=10^-(n-p+precision);
%     pulse_steps=Pulse_duration/dt;
%     %col+1 because we need a refresh signal first
%     steps=pulse_steps;
%     V_pulse = Pulse_Amplitude*(sigmf((1:pulse_steps),[0.01 600])-sigmf((1:pulse_steps)-(10e-9/dt),[0.01 0]));
Pulse_Width=param(20);   
Pulse_duration=Pulse_Width+0.2*Pulse_Width;
Va=0.1;
pulse_steps=Pulse_duration/dt;
steps=pulse_steps;
V_pulse = Va*(sigmf((1:pulse_steps),[0.01 600])-sigmf((1:pulse_steps)-(10e-9/dt),[0.01 0]));
%----------------------------------------------------------------------

%Input Initial Qubit State as Voltage Vector-------------------------------
V_input = kron(Qubit_State,V_pulse);
Gain_out=V_input;
%--------------------------------------------------------------------------

Gainward=10200;
% Gainward=1.0176e+04;
for i=1:gate_num
    
    Gain=(coefficients(i))*Gainward;
    %Control Sequence of every Transistor column---------------------------
    T=ones(col,1);
    VT = 5*kron(T,V_pulse/Va);
    %----------------------------------------------------------------------
    
    VR=zeros(1,col);
    IM=zeros(row,col);
    Temp=298*ones(row,col);
    
    for j=1:steps
        [I_out,g(:,:,i),VR,Temp,IM]  = Crossbar_Analysis_Current(Gain_out(:,j)',VT(:,j)',dt,g(:,:,i),VR,Temp,IM,param);
        Differential=reshape(I_out,2,col/2)';
        Differential=Differential(:,1)-Differential(:,2);
        Gain_out(:,j)=Gain*Differential;
    end
end

Quantum_Rezult=10*Gain_out;
Output_states=zeros(1,col/2);

% hold on
for i=1:col/2   
% linetypes=['-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.'];
% plot((1:steps),Quantum_Rezult(i,:),linetypes(i),'LineWidth',2);
Output_states(i)=Quantum_Rezult(i,steps/2);
end
% hold off

end
