function [y]=vfpredmodel(NUMFEMALES, TAUS, START_DAY, PHLOEM_TEMPS)
%     Predict MPB phenology for an Arizona population with given 
%     temperature, population inputs, with terms accounting for
%     the accumulation of variance across temperatures and life
%     stages.
%
%     Number of attacking females given as a row vector
%     (NUMFEMALES). Length should be equal and corresponding to 
%     TAUS.
%
%     Attack days (TAUS) given as a row vector in Julian days 
%     corresponding to NUMFEMALES.
% 
%     START_DAY is first Julian day of PHLOEM_TEMPS records.
% 
%     PHLOEM_TEMPS is phloem temperatures as column vector given
%     in C degrees for each hour of a year, starting at midnight
%      on START_DAY. Vector should be 2*24*365=17520 hrs long.
% 
%     Example command use:
%         numfemales = [9 4 1 5 39 73 20 11 1 1 2];
%         taus = [216 218 222 224 226 229 231 233 236 238 240];
%         phloem = ones(1,17520).*25;
%         startday = 168;
%         time2=[startday:(startday+729)];
%         pout=vfpredmodel_fix(numfemales, taus, startday, phloem);
%         figure
%         plot(time2,pout)
%         legend('Oviposition','Egg','L1','L2','L3','L4', ...
%         'Pupae','Adult Emergence')
%  
%% Input Parameters %%

%Matrix of rates parameters for southern population
global p       %Set p to global so it can be used by ovipos()
% p = [omega   Psi        Tb       Tm       DeltaB    DeltaM]
%     except for p(8,:), teneral adult
p =[0.0913    0.0309    6.6000   30.9000    1.3613    1.9889 
    %Oviposition
    0.2045    0.0326    6.0251   31.9309    0.5410    5.5031
    %Eggs
    0.1517    0.0521    4.6029   31.7661    0.0117    5.4256
    %First Instar
    0.1374    0.0431    5.9791   31.8337    0.0413    4.4534
    %Second Instar
    0.1856    0.0170    6.0115   31.2656         0    4.3079
    %Third Instar
    0.1694    0.0545   14.9990   31.4364         0    5.2947
    %Fourth Instar
    0.1658    0.0166    6.3504   30.8041         0    3.5426
    %Fifth Instar
%   r_max       T_b       T_m
    0.0197    11.3539   27.2079   NaN          NaN      NaN];      
    %Teneral adult (Brière curve)

% Set variance parameters
sigma=[ 0.032               %Oviposition         
    0.038                   %Eggs          
    -0.2182                 %First Instar
    0.1524                  %Second Instar
    0.165                   %Third Instar
    0.1354                  %Fourth Instar
    0.0673                  %Pupae
    0.0009];                %Teneral Adult
                 
%Calculation of nus from sigma
nus= (sigma.^2)./2;

%% Defining time vectors %%
tmps=PHLOEM_TEMPS;                                      
%Define temperature as input temperature
nh=24;                                                  
%Number of hours in a day
ndays=365+365;                                          
% 2*365 = 730, 2 years of time
nt=nh*ndays;                                            
%Total number of hours in 730 days
tmin=START_DAY;                                         
%Smallest time is start day of phloem temps
tmax=ndays+tmin;                                        
%Max day is two years from start date
dt=(tmax-tmin)/nt;                                      
%Time step (1hr)
tmps=reshape(tmps,nh,ndays);                            
%Reshape temperatures so each column is a whole day of data
treal=linspace(tmin,tmax-(tmax-tmin)/ndays,ndays);      
%Real time vector that begins at tmin
tcalc=treal-tmin;                                       
%Time vector that starts at zero - convenient for setting up 
%convolution matrix
tol=1e-8;                                               
%Tolerance used for numerical calculations


%% Lifestage calculations

%Use oviposition model for first input
pout=zeros(8,ndays);                                            
%Initialize matrix for output of emergence distributions
pinput = ovipos(NUMFEMALES,TAUS,PHLOEM_TEMPS,START_DAY);        
%Call oviposition model for initial input
pout(1,:)= pinput;                                              
%First distribution is ovipositing adults (istg=1)

%Calculate rates and integrate variabilities for each stage
for istg=2:8
    
    %get rates for current stage:
        %eggs
    if (istg == 2)                             
        rates = getrates(tmps, p(2,:));
        %L1
    elseif (istg == 3)
        rates = getrates(tmps, p(3,:));
        %L2
    elseif (istg == 4)
        rates = getrates(tmps, p(4,:));
        %L3
    elseif (istg == 5)
        rates = getrates(tmps, p(5,:));
        %L4
    elseif (istg == 6)
        rates = getrates(tmps, p(6,:));
        %Pupae
    elseif (istg == 7)
        rates = getrates(tmps, p(7,:));
        %Adult emergence/teneral adult
    elseif (istg == 8)
        rates=(briere(tmps,p(8,:)));
     
    end %of calculating rates for current stage
    
    rates=sum(rates);               
    % add up developmental increment for this day by summing over 
    % hours in the day (down columns)
    crates=dt*cumtrapz(rates);      
    % cumulative rates over days for this stage
    crates1=1-crates;               
    % to make the calculation efficient
    nu=nus(istg);                   
    % variance for this stage
    
    %pre-calculate some factors which will be used over and over 
    %in loop for efficiency
    texp(2:ndays)=1./(4*nu*tcalc(2:ndays));                 
    %1/denominator of exp of Green's function
    tden(2:ndays)=1./sqrt(4*nu*pi*tcalc(2:ndays).^3);       
    %1/sqrt(stuff) of Green's function
    
   
    % the following loop  integrates the contribution of variance
    % for non-constant temperatures by summing along diagonals	    
    % in an upper triangular matrix

    for i=1:ndays-1
        
        i1=i+1:ndays;
        
        % pttau is the Green's function which weights the 
        % contributions of variances and development
        pttau=exp(-texp(i1-i+1).*(crates1(i1)+crates(i)).^2).*tden(i1-i+1);
        
        %because of singularity in pttau it is wise to normalize
        wts=trapz([0 pttau]);
        if (wts<tol)	
            % don't normalize those pttaus which are just small
            wts=1;
        end		% of normalization
        
        pout(istg,i1)=pout(istg,i1)+pttau*pinput(i)/wts;
        
    end			% of integrating variance
    
    pinput=pout(istg,:);  
    % input for next life stage is end of last stage
    
end		% of life stage calculations
y=pout; %Output

%% Functions called by vfpredmodel %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y] = ovipos(NUMFEMALES,TAUS,PHLOEM_TEMPS,START_DAY)
%     Creates MPB oviposition PDF for an Arizona population with given 
%     temperature, population inputs, with terms accounting for the
%     accumulation of variance across temperatures. This model is derived
%     from the McManis et al (2019) oviposition rate model. Parameters for
%     getrates() are the Southern MPB population oviposition parameters
%     in this thesis.
% 
%     Number of attacking females given as a row vector (NUMFEMALES).
%     Length should beequal and corresponding to TAUS.
% 
%     Attack days (TAUS) given as a row vector in Julian days corresponding
%     to NUMFEMALES.
% 
%     START_DAY is first Julian day of PHLOEM_TEMPS records.
% 
%     PHLOEM_TEMPS is phloem temperatures as column vector given in C 
%     degrees for each hour of a year, starting at midnight on START_DAY.
%     Vector should extend at least 45 days past TAUS(end).
%
%     Example command use: Using the same inputs as vfpredmodel(), the PDF
%     can be plotted using a time vector created similiarly to lines 196
%     to 201.
%
%% Define variables and Initialize Matrices
   global p                       
% d value for calculation of t0
d = 1;
eggfree = 1;
%Sigma from McManis et al. 2019 paper
q= 0.32;        
%Creating uniform vector of f-values
fvalues = linspace(0.05,0.95,19);   
%Creating vector of taus from input TAUS
taus = TAUS; 
%Vector for number of females created from NUMFEMALES for weighting PDFs
numfemales = NUMFEMALES;            
%Initializing matrices
all = [];                           
final = [];

%% Creating time/temp vector that begins at START_DAY
%defining dt
dt=(1/24); 
%First day of time is start date for phloem temps
tmin= START_DAY; 
%Last day is day of last attack plus 45 days
tmax=TAUS(end)+45; 
%nt is number of blocks needed
nt=(tmax-tmin)/dt+1;   
%Time separated in units of days but in increments of hours
t=linspace(tmin,tmax,nt); 
%Eliminate last block that is extra
t(end)= [];              
 %Create vector of phloem temps for tmin-tmax above
T = PHLOEM_TEMPS(1:length(t));     

%% Calculate rates, cumulative rates, t0 (delay time)
%Calculate oviposition rates
rates = getrates(T,p(1,:)); 
%Get cumulative sum of rates
cumrates = dt.*cumtrapz(rates); 
%Taking the integral that will be later used for delay times
int = cumtrapz(getrates_t0(T)).*dt;                  
%Define taus for the tau loop
usetau = TAUS;                     

%Begin loop which creates PDFs for each tau
for m = 1:length(usetau)
    tau = usetau(1,m);   
    %Find index of tau in t
    tau_index=find(t>=tau,1);         
    %Find the index value for end of egg-free distance excavation
    endeggfreeindex = find((int - int(tau_index))>=eggfree,1);  
    %Use that index to find the day t0 ends, subtract tau to get t0                                                             
    t0 = t(endeggfreeindex)-tau;                                

    %Calculate cumulative rates R with delay time removed
    if isempty(endeggfreeindex)== 1                                           
        R = 0;                                                                                             
        t0 =0;
    else
         R=max(0,cumrates-cumrates(endeggfreeindex));
    end %End of R calculation                                                         
    
 % Loop calculating PDFs for uniform f for current tau
    %Loop to calculate PDFs for uniform f vector
    for k = 1:length(fvalues)
        %Use specific f for loop
        f = fvalues(k);
        %Create tt to avoid division by zero below
        tt=max(dt,(t-tau-t0));    
        %Epsilon used in PDF calculation   
        eps = (R+log(f))./tt; 
        %Derivative of epsilon used in PDF calcuation
        deps = (t>t0+tau).*abs(-(R+log(f))./(tt.^2)+(rates./tt) );
        %PDF calculation
        PDF =(1/sqrt(2.*pi.*q.*q)).*exp(((-1./(2.*q.*q)).*(eps.^2))).*deps;         
        %Normalization of PDF
        nPDF = PDF./trapz(PDF*dt); 
        %Placing PDF for specific f in a matrix
        all(k,:)=nPDF;                                                              
    end %of calculation for uniform f vector
    
    %Sum down columns to get resultant P(some egg)
    p_egg=sum(all);   
    %Normalize PDF
    p_egg=p_egg./(trapz(p_egg*dt)); 
    %Place PDF for this tau in final matrix
    final(size(final,1)+1,:) = p_egg;                                               
       
end %of tau loop 
%% Weighting and normalization of PDF for input tree
%Weighting by number of starting NUMFEMALES
tree = diag(numfemales(1,:))*final; 
%Summing PDF for all taus
tree = sum(tree);
%Normalizing PDF
tree = tree./(trapz(tree)*dt);                                                      

%% Reshape PDF from hours to days and correctly size vector for vfpredmodel
%Reshape vector into matrix where each column is a day
treematrix = reshape(tree,24,tmax-tmin);
%Sum to get PDF in terms of days
treeday = dt*sum(treematrix);  
%Create vector of length required for vfpred
treeinput = zeros(1,730);
%Place final PDF in vector
treeinput(1:(length(treeday)))=treeday;         

[y] = treeinput;                                %Define output


function [rates] = getrates(T,p) 
% Rates function as seen in Régnière et al. (2012)
% and McManis et al. (2018).
% T is temperature vector
% p is vector of parameters from McManis et al. (2018)
% p(1) = omega
% p(2) = Psi 
% p(3) = Tb
% p(4) = Tm
% p(5) = DeltaB
% p(6) = DeltaM
[rates] = max(0, p(2).*((exp(p(1).*(T-p(3))))-(((p(4)-T) ...
./(p(4)-p(3))).*exp((-1*p(1)*(T-p(3)))./p(5)))- ...
(((T-p(3))./(p(4)-p(3))).*exp((p(1).*(p(4)-p(3)))-...
((p(4)-T)./p(6))))));

function [y] = briere(tmps,p)
%Brière curve for teneral adult rates
% p = p(8,:)
% p = (rmax, T_B, T_M)
%Example code: 
% % p = [0.0197 11.3539 27.2079];
% % T = linspace(0,40,100);
% % output = briere(T,p);
% % plot(T,output)
m=2;  %for sqrt Brière
Tb=p(2); a=p(1); Tm=p(3);
% find the curve's max using Brière formula:
Topt=(2*m*Tm+(m+1)*Tb+ sqrt( 4*m^2*Tm^2+(m+1)^2*Tb^2- ...
4*m^2*Tb*Tm))/(4*m+2);
% value at max, to normalize the shape function to peak of 1
norm=Topt.*( Topt -p(2) ).*(( abs(1-(Topt./p(3))) ).^0.5);
% now the max rate is exactly the parameter a=p(1);
r = a.*tmps/norm.*( tmps -Tb ).*(( abs(1-(tmps./Tm)) ).^0.5);
output = r.*(tmps>=p(2)).*(tmps<=p(3)); 
% make sure output is zero outside developmental range
[y] = output;

function [rates] = getrates_t0(tmps) 
% Rates function as seen in Régnière et al (2012) and McManis et al (2009).
% T is temberature vector
% b is vector of barameters from McManis thesis bg 62 (2018)
% b(1) = omega
% b(2) = Psi 
% b(3) = Tb
% b(4) = Tm
% b(5) = DeltaB
% b(6) = DeltaM
b = [0.0632    0.1773    5.8992   29.6069  2.5514   2.7269];
[rates] = max(0, b(2).*((exp(b(1).*(tmps-b(3))))-(((b(4)-tmps)...
./(b(4)-b(3))).*exp((-1*b(1)*(tmps-b(3)))./b(5)))-(((tmps-b(3))...
./(b(4)-b(3))).*exp((b(1).*(b(4)-b(3)))-((b(4)-tmps)./b(6))))));