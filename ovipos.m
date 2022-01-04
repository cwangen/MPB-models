
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

%Oviposition parmaters
p =[0.0913    0.0309    6.6000   30.9000    1.3613    1.9889]                            
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
% Rates function as seen in Régnière et al (2012) and McManis et al (2009).
% T is temperature vector
% p is vector of parameters from McManis thesis pg 62 (2018)
% p(1) = omega
% p(2) = Psi 
% p(3) = Tb
% p(4) = Tm
% p(5) = DeltaB
% p(6) = DeltaM
[rates] = max(0, p(2).*((exp(p(1).*(T-p(3))))-(((p(4)-T)./(p(4)-p(3)))... 
.*exp((-1*p(1)*(T-p(3)))./p(5)))-(((T-p(3))./(p(4)-p(3)))...
.*exp((p(1).*(p(4)-p(3)))-((p(4)-T)./p(6))))));

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
