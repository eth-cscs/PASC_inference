>> [af, sf] = filters1;     	% analysis and synthesis filters
>> J = 5;                       % five stages
>> L = 8*2^(J+1);		
>> N = L/2^J;
>> x = zeros(L,8*L);            % create zero signal
>> w = double_f2D(x,J,af);      % analysis filter bank (5 stages)
>> w{J}{1}(N/2,N/2+0*N) = 1;	% set single coefficient to 1 for each of 
>> w{J}{2}(N/2,N/2+1*N) = 1;	% the eight 2-D wavelets
>> w{J}{3}(N/2,N/2+2*N) = 1;
>> w{J}{4}(N/2,N/2+3*N) = 1;
>> w{J}{5}(N/2,N/2+4*N) = 1;
>> w{J}{6}(N/2,N/2+5*N) = 1;
>> w{J}{7}(N/2,N/2+6*N) = 1;
>> w{J}{8}(N/2,N/2+7*N) = 1;
>> y = double_i2D(w,J,sf);		% synthesis filter bank (5 stages)
>> figure(1)		
>> clf
>> imagesc(y);                  % construct an image of all eight wavelets
>> axis equal	
>> axis off
>> colormap(gray(128))          % convert to a grayscale image