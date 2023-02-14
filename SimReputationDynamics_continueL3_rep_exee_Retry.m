function[TIME,MEndAveM]=SimReputationDynamics_continueL3_rep_exee_Retry(N,p,q,nIt,l)
% Simulates the reputation dynamics as described in Section 1.2.
% X, Matrix containing the estimated pairwise cooperation frequencies
% M, Estimated average image matrix
% MEnd, Image Matrix in the end of the simulation run.
% AssRule, Matrix that contains all assessment rules present in population
% ActRule, Matrix that contains all action rules present in the population
% PopComp, Vector that specifies how many players use each strategy
% ep, Constant error probability to commit a perception error
% q, Constant probability to observe a third-party interaction
% nIt, number of iterations of the indirect reciprocity game.
START=datetime;
%% Setting up the objects
w=0.015;
for D=1:2
if D==1
    DD=2;
else
    DD=10;
end
d=0.01*DD;
for r=1:l
% N=sum(PopComp); % N, Size of populatio
% nS=length(PopComp); % nS, Number of different strategies
MEnd=zeros(1,N); % Initializing the output
MC=ones(N,N); % Current image matrix; initially everyone is good
% xP=ones(1,PopComp(1)); for j=2:nS, xP=[xP, j*ones(1,PopComp(j))]; end
% N-dim vector. The i-th entry is the strategy index of player i.
%% Simulating the interactions
for t=1:nIt
%% Choosing a donor and a recipient and letting them interact
Do=randi(N); % Selecting a donor
Re=Do; while Re==Do, Re=randi(N); end % Selecting a different receiver
stD=MC(Do,Do); stR=MC(Do,Re); % Defining the players¡¯ standings
if Do<=N*p  % If donor is mutant
%ActRule=(1-u)*(w*stD+stR-(d+w)*stD*stR)+u*(w*stD+stR-(d+w)*stD*stR);
ActRule=w*stD+stR-w*stD*stR;
%ActRule=1+(-1+w)*stD+(1-w)*stD*stR;
else
%ActRule=(1-u)*(w*stD*(1-stR)+stR)+u*(w*stD*(1-stR)+stR);
ActRule=w*stD+stR-w*stD*stR;
%ActRule=1+(-1+w)*stD+(1-w)*stD*stR;
end
cp=ActRule; % cp=1 if donor cooperates
%% Updating the donor¡¯s reputation
for Obs=1:N % Going through all all individuals as potential observers
if Obs==Do || Obs==Re || rand(1)<q % If individual observes interaction
stD=MC(Obs,Do); stR=MC(Obs,Re); % Retrieving the players¡¯ standings
if Obs<=N*p
%AssRule=cp+w*stD*stR-(d+w)*stD*cp*stR; %IS, Ayz=0
%AssRule=-stR+w*stD*stR+cp*stR+(-d-w)*stD*cp*stR+1;%L3 Ayz=1
%AssRule=stD+cp-stD*cp+(w-1)*stD*stR+(1-d-w)*stD*cp*stR; %L1
%AssRule=cp*(stD-1)-stR+w*stD*stR+2*cp*stR+(-1-d-w)*stD*cp*stR+1; %L4
%AssRule=stD+(-1+w)*stD*stR+cp*stR+(-d-w)*stD*cp*stR; %L7
%AssRule=stD+cp-2*stD*cp+(w-1)*stD*stR+(2-w-d)*stD*cp*stR; %L2
%AssRule=-stD*cp-stR+w*stD*stR+cp*stR+(1-d-w)*stD*cp*stR+1; %L5
%AssRule=-cp-stR+w*stD*stR+2*cp*stR+(-d-w)*stD*cp*stR+1; %L6
AssRule=stD-stD*cp+(-1+w)*stD*stR+cp*stR+(1-d-w)*stD*cp*stR; %L8
%AssRule=(1-u)*(-stR+w*stD*stR+cp*stR-(d+w)*stD*cp*stR+1)+u*(cp+w*stD*stR-(d+w)*stD*cp*stR); %(1-u)*SS+u*IS
else
%AssRule=cp+w*stD*stR*(1-cp); %Ayz=0
%AssRule=-stR+w*stD*stR+cp*stR-w*stD*cp*stR+1; %Ayz=1
%AssRule=stD+cp-stD*cp+(w-1)*stD*stR+(1-w)*stD*cp*stR; %L1
%AssRule=cp*(stD-1)-stR+w*stD*stR+2*cp*stR+(-1-w)*stD*cp*stR+1; %L4
%AssRule=stD+(-1+w)*stD*stR+cp*stR-w*stD*cp*stR; %L7
%AssRule=stD+cp-2*stD*cp+(w-1)*stD*stR+(2-w)*stD*cp*stR; %L2
%AssRule=-stD*cp-stR+w*stD*stR+cp*stR+(1-w)*stD*cp*stR+1; %L5
%AssRule=-cp-stR+w*stD*stR+2*cp*stR-w*stD*cp*stR+1; %L6
AssRule=stD-stD*cp+(-1+w)*stD*stR+cp*stR+(1-w)*stD*cp*stR; %L8
%AssRule=(1-u)*(-stR+w*stD*stR+cp*stR-w*stD*cp*stR+1)+u*(cp+w*stD*stR*(1-cp)); %(1-u)*SS+u*IS
end
if AssRule > 1
AssRule = 1; 
end
MC(Obs,Do)=AssRule;
end
end
end
%% Calculating output variables by averaging over all players with same strategy
MEnd=MC;% Image matrix in the end is the last current image matrgvix
%MEnd00=sum(sum(MEnd(1:N/2,1:N/2)))/(N*N/4);
%MEnd01=sum(sum(MEnd(1:N/2,N/2+1:N)))/(N*N/4);
%MEnd10=sum(sum(MEnd(N/2+1:N,1:N/2)))/(N*N/4);
%MEnd11=sum(sum(MEnd(N/2+1:N,N/2+1:N)))/(N*N/4);
MEnd00=sum(sum(MEnd(1:N*p,1:N*p)))/(N^2*p^2);
MEnd01=sum(sum(MEnd(1:N*p,N*p+1:N)))/(N^2*p*(1-p));
MEnd10=sum(sum(MEnd(N*p+1:N,1:N*p)))/(N^2*p*(1-p));
MEnd11=sum(sum(MEnd(N*p+1:N,N*p+1:N)))/(N^2*(1-p)^2);
MEndAve(r,:)=[MEnd00 MEnd01 MEnd10 MEnd11];
end
%MEndAveM(10*(W-1)+D,:)=[w,d,mean(MEndAve)];
MEndAveM(D,:)=[d,mean(MEndAve),std(MEndAve)/sqrt(l)];
end
FIN=datetime;
TIME=FIN-START;