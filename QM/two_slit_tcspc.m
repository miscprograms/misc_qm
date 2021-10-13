% 2-slit expt

%     |          |
%     |          |
%     |          |A
%     S1         |  A'
%     |          |   x = A' to A
%     | d        |M
%     |          |
%     S2         |
%     |          |B
%     |          |
%     |    L     |

% distance S1 to S2 =d; distance from A' to A=x
% distance S1 to A= d1=sqrt(L.^2+ x.^2)  
% distance S2 to A= d2=sqrt(L.^2+ (x+d).^2)  

% theta= angle between mid point of slits and A.
% 2 d sin(theta)= n* lamda

% x varies from [-2d , d]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time correlated single photon counting (TCSPC)

% Source sine wave of frequency "f" and amplitude "A", modelled as a collection of large number of photons,
% each of width "dt << T=1/f" and amplitude "dA << A"
        %we assume a DC offset in sine wave, so that Total number of photons at any time instant is always positive.
        %Rods and cones in our eyes respond to carrier wavelength, but filter out DC offset.
        
% randomly choose screen location for each particle or set of particles
% arriving at slit. uniform probability density function. Similar to Coin Toss experiment
  
%Screen/detector accumulates particles, like in a photographic film, in Tonomura's electron interference experiment   
%TCSPC :Each carrier frequency cycle divided into 50 time bins. Particles arriving at screen are accumulated at suitable time bin.

% As we keep running simulation, more and more Accumulated particles will continue to retain sinusoidal interference pattern, 
%but will saturate screen display because they exceed max value 256. 
% Same thing happens with Tonomura electron interference experiment, so they cut off video after a short time.

% This simulation is for light interference. For electron interference, use
% suitable values for particle speed and wavelength.

%Video tutorial for 2-slit quantum interference simulation using Time Correlated Single Photon Counting (TCSPC)
%https://www.ocf.berkeley.edu/~araman/files/papers/two_slit_tcspc_video.mp4

c=3e8; %light speed in vacuum

lamda=500e-9; %wavelength 500 nm
f=c/lamda; %light carrier frequency



d=lamda*4; %distnace between 2 slits
L=d*10; %distance of screen/detetcor from slit plane

x_step=d/10; 
x_max=d*10;
x_min=-x_max;

y_min=x_min; y_max=x_max;y_step=x_step;


xArray=[x_min: x_step: x_max];
yArray=[y_min: y_step: y_max];
M=length(xArray);
N=length(yArray);

%simulate photons one by one
N0=50; %20; %100;
numTimeBins=N0;

T0=1/f;t_step=T0/N0;
t=[0:t_step:T0-t_step];


% timeBinsArray Simulates screen/detector (row,col, time instant)
timeBinsArray=zeros(M,N,N0);
screenArray=zeros(M,N); %screen/detector (row,col)

numLoop=M*N*N0*1000; %1e6;
timeCurrent=0;

slideCount=0;
screenArray_save=zeros(100,M,N);

slit_1_signal=zeros(1,N0+1);
slit_2_signal=zeros(1,N0+1);

figure(1)
hold off
colormap(gray(256)); 


%Repeat loop 
for loopCount=1:100*10000, %*10, %*50, %280000, %100*10000, %numLoop,
    
   for sampleCount=1:N0,
        
    %Generate sinusoidal RV. Used only for computing time of arrival of each particle at detector.
    % light signal= sin(2*pi*f*t)= sin(2*pi*f*n/fs)  = sin(2*pi*n/N0).
    % fs/f= N0. n=[0,1,..N0-1]  N0=50 time bins within a single carrier cycle
    rv_1= rand(1,1)*2-1; % rv_1 has values in [-1, 1] range. rv_1 = sin(2*pi*n/N0)
    rv_2=asin(rv_1);     % rv_2 = 2*pi*n/N0
    
    %sampleNum= round((rv_2+pi/2)*N0/pi); % sampleNum = 2*n + N0/2.  n=[0,1,..N0-1]
    
    sampleNum= round((rv_2)*N0/(2*pi)); % sampleNum = n.  n=[0,1,..N0-1]
    
    if sampleNum==0, 
        sampleNum=sampleNum+1;
    end
    
    timeCurrent=t_step*(sampleNum-1); %sampleNum n=[0,1,..N0-1] randomly chosen as per above lines. N0=50 time bins
    %timeCurrent is time of arrival of particle at slit, chosen randomly.
    
    for sample=sampleNum, %N0,
        
        %randomly choose screen location for each particle or set of particles arriving at slit
        row=round(rand(1,1)*M);
        col=round(rand(1,1)*N);
        if row==0
          row=1;
        end
        if col==0
           col=1;
        end
        
        %find time bins and number of photons reaching detector time bin
        x=xArray(row); y=yArray(col);
        d1=sqrt(L.^2+ (x+d/2).^2+y.^2); %slit S1
        t1=timeCurrent+d1/c; 
        d2=sqrt(L.^2+ (x-d/2).^2+y.^2);  %S2
        t2=timeCurrent+d2/c;
        
        
        D1=1;
        D2=1;


        
        %Choose slit number randomly
        slitNum=rand(1,1)>.5;
        
        %we assume a DC offset in sine wave, so that Total number of photons at any time instant is always positive.
        %Rods and cones in our eyes respond to carrier wavelength, but filter out DC offset.
 
        %we assume 10 particles in each slit, to speed up simulation. Else,simulation is 10 times slower.
        %they get attenuated to one or two particles when they travel from slit to a specific location in screen
        K1=5*2;  

        
        % 1/R.^2 attenuation not used in this simulation
        if 0
        K1= K1*7e9;
        
        D1=4*pi*d1.^2; %1; %(2*pi*d1).^2;
        D2=4*pi*d2.^2; %1; %(2*pi*d2).^2;
        
        %D1=D1/(4*pi*L.^2);
        %D2=D2/(4*pi*L.^2);
        end        
        
%Screen/detector is timeBinsArray(row,col,T1)

      %Slit 1

        if slitNum==1,
             
            slit_1_signal(sample)=1*K1;slit_2_signal(sample)=0;
            N1=round(slit_1_signal(sample)/D1);
            
            %TCSPC done here. Each carrier frequency cycle divided into 50
            %time bins. Particles arriving at screen are accumulated at suitable time bin.
            T1=round(mod(t1,T0)*numTimeBins/T0); if T1==0, T1=N0; end 
        
            %accumulate particles in timeBinsArray(row,col,T1)
            value=timeBinsArray(row,col,T1) + N1;
            
            %if accumulated particles > 256, satruate at 256.
            if value>256,
               value=256;
            end
             
            timeBinsArray(row,col,T1)= value;

        end
        
      %Slit 2
        if slitNum==0,
            
            slit_1_signal(sample)=0;slit_2_signal(sample)=1*K1;
            N2=round(slit_2_signal(sample)/D2);
 
           
            %TCSPC done here. Each carrier frequency cycle divided into 50
            %time bins. Particles arriving at screen are accumulated at suitable time bin.

            T2=round(mod(t2,T0)*numTimeBins/T0); if T2==0, T2=N0; end 
        
            %accumulate particles in timeBinsArray(row,col,T1)
            value=timeBinsArray(row,col,T2) + N2;

            %if accumulated particles > 256, satruate at 256.
            if value>256,
               value=256;
            end
            
            timeBinsArray(row,col,T2)= value;
            
        end
       

        
        %Demodulation done by Rods and Cones in eye: I cos(wt)+Qsin(wt) => sqrt(I.^2 + Q.^2)
          
        %Take the maximum of timeBinsArray samples, accumulated as per TCSPC, within a single carrier frequency period T0= 1/f0
        if nnz(timeBinsArray(row,col,1:N0))==N0,
        
           sum1= round(max(timeBinsArray(row,col,1:N0)));
            
            if sum1==0, sum1=1; end
            
            screenArray(row,col)=uint8(sum1); %*5;
                   
        end %if nnz()
        

    end %for sample
    
    
  end %sampleCount

%Display accumulated particles over time. Screen/detector accumulates particles, like in a photographic film, in Tonomura's electron interference experiment        
    if mod(loopCount,10000)==0,
       loopCount

       
       slideCount=slideCount+1;
       if slideCount<=100,
       screenArray_save(slideCount,:,: )=screenArray;
       end
       

        hold off
        image(transpose(screenArray(:,:)));
        axis image;

        pause(0.1); %(1); %(0.01);
        %pause

    end
    

end %for loopCount

%save screenArray_save;




