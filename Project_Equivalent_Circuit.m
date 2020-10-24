
%% PART1
% Initialize workspace, load the E2 circuit model as well as the E2 dynamic data
addpath readonly
load readonly/E2model.mat; % load parameter values already created for the E2 cell -- this is a single R-C model
% load readonly/E2model2RC.mat; % this is a two R-C model of the E2 cell
load readonly/E2_DYN_P25.mat; % load raw test data for the E2 cell at 25 degC

% Resample at consistent 1Hz rate.
deltaT = 1; 
time = DYNData.script1.time - DYNData.script1.time(1);    
t = (0:deltaT:time(end));
voltage = interp1(time,DYNData.script1.voltage,t);
current = interp1(time,DYNData.script1.current,t);
time = t;

% function [vk,rck,hk,zk,sik,OCV] = simCellTemp(ik,temp,deltaT,model,z0,iR0,h0)
%
% ik - current in amperes, where (+) is discharge. Size is N x 1.
% temp  - temperature (degC). Size is N x 1.
% deltaT = sampling interval of data in seconds. Size is 1 x 1 (a scalar)
% model - standard model structure
% z0 - initial SOC. Size is 1 x 1.
% iR0 - initial resistor currents as column vector. Size is Nr x 1 where Nr is 
%       number of R-C pairs in model.
% h0 - initial hysteresis state. Size is 1 x 1.
%
% vest - predicted cell voltage. Size is N x 1.
% rck - predicted resistor currents. Size is N x Nr (first row is set to iR0')
% hk - predicted dynamic hysteresis states. Size is N x 1 (first entry is h0)
% zk - predicted cell state of charge. Size is N x 1 (first entry is z0)
% sik - sign of input current. Size is N x 1.
% OCV - predicted cell open circuit voltage. Size is N x 1.
function [vest,rck,hk,zk,sik,OCV] = simCellTemp(ik,temp,deltaT,model,z0,iR0,h0)

  % Force data to be column vector(s) in case user entered data incorrectly
  ik = ik(:); iR0 = iR0(:); temp = temp(:);
  N = length(ik); Nr = length(iR0);
  % initialize some outputs
  vest = zeros(N,1); rck = zeros(N,Nr); hk = zeros(N,1); zk = zeros(N,1); 
  sik = zeros(N,1); OCV = zeros(N,1);
  rck(1,:) = iR0'; hk(1) = h0; zk(1) = z0; sik(1) = 0;
  OCV(1) = OCVfromSOCtemp(z0,temp(1),model);


  
  for k = 2:length(ik),
  
  T = temp(k-1); % sample code uses only single temperature -- you will need to change this!

  RCfact = exp(-deltaT./abs(getParamESC('RCParam',T,model)))';
  if length(RCfact) ~= Nr,
    error('iR0 does not have the correct number of entries');
  end
  G = getParamESC('GParam',T,model);
  Q = getParamESC('QParam',T,model);
  M = getParamESC('MParam',T,model);
  M0 = getParamESC('M0Param',T,model);
  RParam = getParamESC('RParam',T,model);
  R0Param = getParamESC('R0Param',T,model);
  etaParam = getParamESC('etaParam',T,model);
  
    etaik(k-1) = ik(k-1); 
  if (etaik(k-1) < 0),
    etaik(k-1) = etaParam*ik(k-1);
  end

  % Simulate the dynamic states of the model
    rck(k,:) = rck(k-1,:)*diag(RCfact) + (1-RCfact')*etaik(k-1);
    zk(k) = zk(k-1) - etaik(k-1)*deltaT/(Q*3600);
  end
  
  if any(zk>1.1),
    warning('Current may have wrong sign as SOC > 110%');
  end
  
  % Hysteresis stuff
  
  for k=2:length(ik),
  T = temp(k-1); % sample code uses only single temperature -- you will need to change this!
  G = getParamESC('GParam',T,model);
  Q = getParamESC('QParam',T,model);
      
    fac(k-1)=exp(-abs(G*etaik(k-1)*deltaT/(3600*Q)));
    
    hk(k)=fac(k-1)*hk(k-1)+(fac(k-1)-1)*sign(ik(k-1));
    sik(k) = sign(ik(k));
    if abs(ik(k))<Q/100, sik(k) = sik(k-1); end
  end
    
  % Compute output equation
  for k = 1:length(ik),
    T = temp(k); 
    M = getParamESC('MParam',T,model);
    M0 = getParamESC('M0Param',T,model);
    RParam = getParamESC('RParam',T,model);
    R0Param = getParamESC('R0Param',T,model);
  
    OCV(k) = OCVfromSOCtemp(zk(k),T,model);
  
    vest(k) = OCV(k) - rck(k,:)*RParam' - R0Param*ik(k) + M*hk(k) + M0*sik(k);
  end
  

end


% Execute simCellTemp to determine voltage and other internal states/variables

temp = 25*ones(size(current)); 
% temp = linspace(25,45,length(current)); % uncomment to simulate temperature ramp 

[vest,rck,hk,zk,sik,OCV] = simCellTemp(current,temp,deltaT,model,1,0,0);


subplot(1,2,1)
plot(time/3600,voltage,time/3600,vest); % factor of 3600 converts seconds -> hours
xlabel('Time (hr)'); ylabel('Voltage (V)'); title('Comparing measured to simulated voltage');
legend('Measured voltage','Simulated voltage');

% Now, plot the voltage prediction error
subplot(1,2,2)
plot(time/3600,1000*(voltage(:)-vest(:)));
xlabel('Time (hr)'); ylabel('Voltage (mV)'); title('Voltage prediction error');


%% PART2


% Initialize workspace, load the E2 circuit model as well as the E2 dynamic data
addpath readonly
load readonly/E2model.mat; % load parameter values already created for the E2 cell
load readonly/E2_DYN_P25.mat; % load raw test data for the E2 cell at 25 degC

% Resample at consistent 1Hz rate.
deltaT = 1; 
time = DYNData.script1.time - DYNData.script1.time(1);    
t = (2000:deltaT:2500); 
current = interp1(time,DYNData.script1.current,t);



% function [vpack,vcell,icell,zcell,qcell,rcell] = simPCMTemp(Ns,Np,current,temp,deltaT,model)
%
% Simulate parallel-connected-module packs (cells are connected in parallel
% to make modules; these modules are connected in series to make packs).
% Note: This function must work for models having a single R-C pair. It is
% not required to work for models having multiple R-C pairs.
%
% Ns - number of modules connected in series to make a pack
% Np - number of cells connected in parallel in each module
% current - battery pack current, where (+) is discharge. Size is N x 1.
% temp  - temperature (degC). Size is N x 1 (all cells at same temperature).
% deltaT = sampling interval in data (s)
% model - standard model structure
% delta - variability vector: variable initial SOC if delta(1)==1; 
%                         variable total capacity if delta(2)==1;
%                         variable series resistance if delta(3)==1;
%
% vpack - battery pack voltage. Size is N x 1.
% vcell - individual cell voltages. Size is N x Ns x Np
% icell - individual cell currents. Size is N x Ns x Np
% zcell - individual cell states of charge. Size is N x Ns x Np
% qcell - individual cell capacities. Size is N x Ns x Np
% rcell - individual cell series resistances. Size is N x Ns x Np

function [vpack,vcell,icell,zcell,qcell,rcell] = simPCMTemp(Ns,Np,current,temp,deltaT,model,delta)

  % Force current to be column vector in case user entered data incorrectly
  current = current(:); N = length(current); temp = temp(:);
  % Initialize function outputs
  vpack = zeros(N,1); vcell = zeros(N,Ns,Np); icell = zeros(N,Ns,Np);
  zcell = zeros(N,Ns,Np); qcell = zeros(N,Ns,Np); rcell = zeros(N,Ns,Np);
  
  % Do some error checking on the function inputs
  Nr = length(getParamESC('RCParam',25,model)); % number of R-C pairs.
  if Nr ~= 1,
    error('This code does not work for models having multiple R-C pairs.');
  end
  if length(temp) ~= N,
    error('Input "temp" vector not the correct dimension.');
  end
  
  % Initialize states for ESC cell model
  if delta(1),
    z = reshape(linspace(0.3,0.7,Ns*Np),Ns,Np); % Different initial SOCs
  else
    z = 0.5*ones(Ns,Np);
  end
  irc = zeros(Ns,Np);
  h   = zeros(Ns,Np);

  for k = 1:N,
  
  T = temp(k); % sample code uses only single temperature -- you will need to change this!
  
  
  % Default initialization for cells within the pack
  q  = getParamESC('QParam',T,model)*ones(Ns,Np); 
  rc = exp(-deltaT./abs(getParamESC('RCParam',T,model)))'*ones(Ns,Np);
  r  = (getParamESC('RParam',T,model))';
  m  = getParamESC('MParam',T,model)*ones(Ns,Np);
  g  = getParamESC('GParam',T,model)*ones(Ns,Np);
  r0 = getParamESC('R0Param',T,model)*ones(Ns,Np); 
  rt = 0.000125; % 125 microOhm resistance for each tab

  
  if delta(2), 
    q = reshape(linspace(0.95,1.05,Ns*Np),Ns,Np).*q; 
  end
  
  if delta(3),
    r0 = reshape(linspace(0.95,1.05,Ns*Np),Ns,Np).*r0; 
  end
  r0 = r0 + 2*rt; % add tab resistance to cell resistance

  
  
    v = OCVfromSOCtemp(z,T,model); % get OCV for each cell: Ns * Np matrix
    v = v + m.*h - r.*irc; % add in capacitor voltages and hysteresis

    V = (sum(v./r0,2) - current(k))./sum(1./r0,2);
    ik = (v-repmat(V,1,Np))./r0;

    z = z - (1/3600)*ik./q;  % Update each cell SOC
    irc = rc.*irc + (1-rc).*ik; % Update capacitor voltages
    fac = exp(-abs(g.*ik)./(3600*q));
    h = fac.*h + (fac-1).*sign(ik); % Update hysteresis voltages

    vpack(k)     = sum(V); % Store pack voltage
    vcell(k,:,:) = v - ik.*r0; % Store cell voltages
    zcell(k,:,:) = z; % Store cell SOCs
    icell(k,:,:) = ik; % Store cell currents
    qcell(k,:,:) = q; % Store cell capacities
    rcell(k,:,:) = r0; % Store cell resistances
  end % for k

  

end

% Execute simCell to determine voltage and other internal states/variables

temp = 25*ones(size(current)); % for now, use constant 25 degC temperature.
% temp = linspace(25,45,length(current)); % uncomment to simulate temperature ramp 

Ns = 2; Np = 2; 
[vpack,vcell,icell,zcell,qcell,rcell] = simPCMTemp(Ns,Np,current,temp,deltaT,model,[1,0,0]);



% Plot the individual SOC vs. time for all cells in all 
% series PCMs. There is one subplot for each PCM.
t = (0:(length(zcell(:,:,1))-1))/60; 
xplots = round(1.0*ceil(sqrt(Ns))); yplots = ceil(Ns/xplots); 
for k = 1:Ns,
  zr=squeeze(100*zcell(:,k,:));
  subplot(yplots,xplots,k); plot(t,zr); axis([0 ceil(max(t)) 0 100]);
  title(sprintf('Cells in PCM %d',k)); 
  ylabel('SOC (%)'); xlabel('Time (min)'); 
end



%% PART3


% Initialize workspace, load the E2 circuit model as well as the E2 dynamic data
addpath readonly
load readonly/E2model.mat; % load parameter values already created for the E2 cell
load readonly/E2_DYN_P25.mat; % load raw test data for the E2 cell at 25 degC

% Resample at consistent 1Hz rate.
deltaT = 1; 
time = DYNData.script1.time - DYNData.script1.time(1);    
t = (2000:deltaT:2500); % select short segment to speed up simulation
voltage = interp1(time,DYNData.script1.voltage,t);
current = interp1(time,DYNData.script1.current,t);
time = t;


% function [vpack,vcell,icell,zcell,qcell,rcell] = simSCMTemp(Ns,Np,current,temp,deltaT,model)
%
% Simulate series-connected-module packs (cells are connected in series
% to make modules; these modules are connected in parallel to make packs).
% Note: This function must work for models having a single R-C pair. It is
% not required to work for models having multiple R-C pairs.
%
% Ns - number of cells connected in series to make a module
% Np - number of modules connected in parallel in each pack
% current - battery pack current, where (+) is discharge. Size is N x 1.
% temp  - temperature (degC). Size is N x 1 (all cells at same temperature).
% deltaT = sampling interval in data (s)
% model - standard model structure
% delta - variability vector: variable initial SOC if delta(1)==1; 
%                         variable total capacity if delta(2)==1;
%                         variable series resistance if delta(3)==1;
%
% vpack - battery pack voltage. Size is N x 1.
% vcell - individual cell voltages. Size is N x Ns x Np
% icell - individual cell currents. Size is N x Ns x Np
% zcell - individual cell states of charge. Size is N x Ns x Np
% qcell - individual cell capacities. Size is N x Ns x Np
% rcell - individual cell series resistances. Size is N x Ns x Np

function [vpack,vcell,icell,zcell,qcell,rcell] = simSCMTemp(Ns,Np,current,temp,deltaT,model,delta) %#ok<FNDEF>

  % Force current to be column vector in case user entered data incorrectly
  current = current(:); N = length(current); temp = temp(:);
  % Initialize function outputs
  vpack = zeros(N,1); vcell = zeros(N,Ns,Np); icell = zeros(N,Ns,Np);
  zcell = zeros(N,Ns,Np); qcell = zeros(N,Ns,Np); rcell = zeros(N,Ns,Np);
  
  % Do some error checking on the function inputs
  Nr = length(getParamESC('RCParam',25,model)); % number of R-C pairs.
  if Nr ~= 1,
    error('This code does not work for models having multiple R-C pairs.');
  end
  if length(temp) ~= N,
    error('Input "temp" vector not the correct dimension.');
  end
  
  % Initialize states for ESC cell model
  if delta(1),
    z = reshape(linspace(0.3,0.7,Ns*Np),Ns,Np); % Different initial SOCs
  else
    z = 0.5*ones(Ns,Np);
  end
  irc = zeros(Ns,Np);
  h   = zeros(Ns,Np);

  for k = 1:N,
  
  T = temp(k); % sample code uses only single temperature -- you will need to change this!
  
  
  
  
  % Default initialization for cells within the pack
  q  = getParamESC('QParam',T,model)*ones(Ns,Np); 
  rc = exp(-deltaT./abs(getParamESC('RCParam',T,model)))'*ones(Ns,Np);
  r  = (getParamESC('RParam',T,model))';
  m  = getParamESC('MParam',T,model)*ones(Ns,Np);
  g  = getParamESC('GParam',T,model)*ones(Ns,Np);
  r0 = getParamESC('R0Param',T,model)*ones(Ns,Np); 
  rt = 0.000125; % 125 microOhm resistance for each tab

 
  if delta(2), 
    q = reshape(linspace(0.95,1.05,Ns*Np),Ns,Np).*q; 
  end
  
  if delta(3),
    r0 = reshape(linspace(0.95,1.05,Ns*Np),Ns,Np).*r0; 
  end
  r0 = r0 + 2*rt; % add tab resistance to cell resistance

  
  
    v = OCVfromSOCtemp(z,T,model); % get OCV for each cell: Ns * Np matrix
    v = v + m.*h - r.*irc; % add in capacitor voltages and hysteresis

    V = (sum(sum(v,1)./sum(r0,1),2)-current(k))./sum(1./sum(r0,1),2); % Bus V
    ik = (sum(v,1)-repmat(V,1,Np))./sum(r0,1); % 1*Np cell currents
    ik = repmat(ik,Ns,1); % Ns*Np cell currents

    z = z - (1/3600)*ik./q;  % Update each cell SOC
    irc = rc.*irc + (1-rc).*ik; % Update capacitor voltages
    fac = exp(-abs(g.*ik)./(3600*q));
    h = fac.*h + (fac-1).*sign(ik); % Update hysteresis voltages

    vpack(k)     = V; % Store pack voltage
    vcell(k,:,:) = v - ik.*r0; % Store cell voltages
    zcell(k,:,:) = z; % Store cell SOCs
    icell(k,:,:) = ik; % Store cell currents
    qcell(k,:,:) = q; % Store cell capacities
    rcell(k,:,:) = r0; % Store cell resistances
  end % for k

 

end


% Execute simCell to determine voltage and other internal states/variables

temp = 25*ones(size(current)); % for now, use constant 25 degC temperature.
% temp = linspace(25,45,length(current)); % uncomment to simulate temperature ramp 

Ns = 2; Np = 2; 
[vpack,vcell,icell,zcell,qcell,rcell] = simSCMTemp(Ns,Np,current,temp,deltaT,model,[1,0,0]);



% Plot the individual SOC vs. time for all cells in all 
% series PCMs. There is one subplot for each PCM.
t = (0:(length(zcell(:,:,1))-1))/60; 
xplots = round(1.0*ceil(sqrt(Ns))); yplots = ceil(Ns/xplots); means = [];
for k = 1:Np,
  zr=squeeze(100*zcell(:,:,k));
  subplot(yplots,xplots,k); plot(t,zr); axis([0 ceil(max(t)) 0 100]);
  title(sprintf('Cells in SCM %d',k)); 
  ylabel('SOC (%)'); xlabel('Time (min)'); 
end