%INSERT COIN (OR GATE)
function [g,R,V_write_plot,VT_plot,time_plot] = Crossbar_Generation(GigaGate,param)

[row,col,gate_num]=size(GigaGate);
col=2*col;
g=1.7e-9*ones(row,col,gate_num);
R=zeros(size(g));
dt=1e-12;
% precision=4;
for i=1:size(GigaGate,3)
    %Memristors of Crossbar that need to be set (W matrix)-----------------
    W_1=GigaGate(:,:,i);
    W_2=GigaGate(:,:,i);
    W_1(W_1==-1)=0;
    W_1=kron(W_1,[1 0]);
    W_2(W_2==1)=0;
    W_2=kron(W_2,[0 -1]);
    W=W_1+W_2;
    %----------------------------------------------------------------------
    
    %Write Pulse Generation------------------------------------------------
    Pulse_Amplitude=param(19);
    Pulse_Width=param(20);
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
    %col+1 because we need a refresh signal first
%     steps=(col+1)*pulse_steps;
%     V_pulse = Pulse_Amplitude*(sigmf((1:pulse_steps),[0.01 600])-sigmf((1:pulse_steps)-(10e-9/dt),[0.01 0]));
    p1=linspace(0,Pulse_Amplitude,0.1*(Pulse_Width/dt));
    p2=linspace(Pulse_Amplitude,Pulse_Amplitude,0.8*(Pulse_Width/dt));
    p3=linspace(Pulse_Amplitude,0,0.1*(Pulse_Width/dt));
    p4=linspace(0,0,0.2*(Pulse_Width/dt));
    V_pulse=[p1 p2 p3 p4];
    resolution=100;
    V_pulse_downsampled=downsample(V_pulse,resolution);
    pulse_steps=length(V_pulse);
    steps=(col+1)*pulse_steps;
    %----------------------------------------------------------------------
    
    %Write Sequence of every row-------------------------------------------
    W=[-1*ones(row,1) W];
    V_write = kron(W,V_pulse)+0.001;
    %----------------------------------------------------------------------
    
    %Plot Write equence of every row---------------------------------------
    V_write_plot(:,:,i) = kron(W,V_pulse_downsampled)+0.001;
    time_plot=dt*resolution*(1:length(V_write_plot));
    %----------------------------------------------------------------------
    
    %Control Sequence of every Transistor column---------------------------
    T=eye(col,col);
    T=[1*ones(col,1) T];
    T_amplitude=5;
    VT = T_amplitude*kron(T,V_pulse/Pulse_Amplitude);
    VT_plot(:,:,i)=T_amplitude*kron(T,V_pulse_downsampled/Pulse_Amplitude);
    %----------------------------------------------------------------------
    
    VR=zeros(1,col);
    IM=zeros(row,col);
    Temp=298*ones(row,col);
    
    for j=1:steps
        [~,g(:,:,i),VR,Temp,IM]  = Crossbar_Analysis_Current(V_write(:,j)',VT(:,j)',dt,g(:,:,i),VR,Temp,IM,param);
    end
  
    %Resistance Read-------------------------------------------------------
    I0= param(9);
    g0= param(10);
    V0= param(11);
    VM=0.1;
    I=I0*exp(-g(:,:,i)/g0).*sinh(VM/V0);
    R(:,:,i)=VM./I;
    %----------------------------------------------------------------------
end

end
