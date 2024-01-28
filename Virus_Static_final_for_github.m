%---------------------SVS: Base Balance in DNA Sequence Design----------------------
%Random process, each result will be different

clear all;
clc;
tic
%----------------------initialization------------------------%
Nc=20;     %The number of times a virus evolves within an individual
numDNA=20; %Virus
TotalDNA=numDNA;
NowDNA=numDNA;
VirusnumSave=7;
Tran_Time=10;%Epidemic Time
numIndividual=200; %Host
NowIndividual=numIndividual;
M_max=40; %Lm
N_max=40; %Ln
R0=[];

%Distribution location
m(1:numIndividual/4)=unifrnd(0.1,floor(M_max/2),[1 numIndividual/4]);
n(1:numIndividual/4)=unifrnd(floor(N_max/2),N_max,[1 numIndividual/4]);
m(numIndividual/4+1:numIndividual/2)=unifrnd(floor(M_max/2),M_max,[1 numIndividual/4]);
n(numIndividual/4+1:numIndividual/2)=unifrnd(0.1,floor(N_max/2),[1 numIndividual/4]);
m(numIndividual/2+1:numIndividual-10)=5*rand(1,numIndividual/2-10)+randi([floor(M_max/2) M_max-5],1,numIndividual/2-10);
n(numIndividual/2+1:numIndividual-10)=5*rand(1,numIndividual/2-10)+randi([floor(N_max/2) N_max-5],1,numIndividual/2-10);
m(numIndividual-9:numIndividual)=randperm(floor(M_max/2),10);
n(numIndividual-9:numIndividual)=randperm(floor(M_max/2),10);

Distance_Table=zeros(numIndividual,numIndividual);%Comparison table of distances between all individuals
Trans_Distance=3;%Transmission events only occur when the distance between two individuals is less than 3
W1=0.7;%The weight of DNA score in determining whether an infection has occurred
W2=0.3;%The weight of the probability of transmission events in determining whether an infection has occurred (determined by distance)
I_probability=0.7;%The probability of type I infection after infection
Susceptibility=0.7;%The probability of transitioning to susceptible type after treatment
TD=40;%The threshold for mortality after infection, exceeding which is considered death. The lower the threshold, the higher the mortality rate
immunobarrier=0.8;%Immune barrier establishment index, when the number of immunized individuals reaches this proportion, it is considered to establish an immune barrier
Die_Time_max=Tran_Time;
Cure_Time_max=Tran_Time;
Virus_List=[];%List of Existing Viruses
Individual_List=[];%List of Existing Hosts


diary('E:\DNA_design_Virus_Static\Record.txt');
fprintf('\n')
fprintf('--------------------------------------------------------------------------\n')
fprintf('\n')
fprintf('Begin Time£º\n')
Time=clock;
disp([num2str(Time(1)),'-',num2str(Time(2)),'-',num2str(Time(3)),'   ',num2str(Time(4)),':',num2str(Time(5))])



global DNAInfo   
DNAInfo=struct('Sequence',[],'Generation',[],'Similarity',[],'H_measure',[],'Continuity',[],'GCcontent',[],'Hairpin',[],'Tm',[],'Goal',[],'Transmission',[],'Virulence',[],'Compare_Similarity',[],'Compare_H_measure',[]);
DNAInfo=repmat(DNAInfo,1,numDNA);


%Creating annotations for individual information structures (which can be considered as individuals) with a large amount of information
%Virus, if Virus=0, then the individual is not infected and has no virus in their body. If Virus=x, then the infected person carries virus number x
%Type itself, one of the six categories;
%Trans_ Num and Trans_ The number and distance of individuals within the distance propagation range are determined by their two-dimensional spatial position after initialization;
%Protect uses protective equipment to reduce the risk of infection;
%Location_ X and Location_ The position of Y itself, coordinates (m, n), is a fixed value in the static model, and is yet to be determined in the dynamic model, so it will not be considered temporarily;
%Die_ Time and Die_ Probability of death from illness, time and probability of death from infection, random number;
%Die_ Time_ Click the countdown to death, equal to 0 when death occurs
%Die_ Will it die? 1 means it will, 0 means it won't
%T2toT1_ Time and T2toT1_ If the probability is a type II infection, the probability and time of conversion to type I are random numbers;
%Cure healing time, random number. Simultaneously used as a timer, when 0, cured
%Infection_ The number of times time has been infected
%Explanation of six types:
%0=susceptible, not yet infected, but susceptible to virus infection, is the target of virus search;
%1=Infectious infection (Type I), already infected, with the ability to spread, can and can only infect susceptible individuals;
%2=Non communicable infection (type II), already infected and not contagious;
%3=immune, unable to be infected;
%4=cured, the virus is killed, and then randomly transformed into immune individuals or susceptible individuals;
%5=Death due to illness, the host is killed, and the virus carried by the individual apoptosis, making it impossible for the deceased individual to infect again.
global IndividualInfo  
IndividualInfo=struct('Virus',[],'Type',[],'Trans_who',[],'Trans_num',[],'Trans_distance',[],'Protect',[],'Location_X',[],'Location_Y',[],'Die_time',[],'Die_probability',[],'Die_time_click',[],'Die_will',[],'T2toT1_time',[],'T2toT1_probability',[],'Cure',[],'OldCureTime',[],'Infection_time',[],'R0',[],'Infstart',[]);
IndividualInfo=repmat(IndividualInfo,1,numIndividual);


initial=[0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3];
for i=1:numDNA
changesort=randperm(size(initial,2));
Dna(i,:)=initial(changesort);
end 

DNA=char(numDNA,20); 
for i=1:numDNA 
    for ii=1:20
        switch Dna(i,ii)
            case 0
                DNA(i,ii)='A';
                num2str(DNA(i,ii));
            case 1
                DNA(i,ii)='G';
                num2str(DNA(i,ii));
            case 2
                DNA(i,ii)='T';
                num2str(DNA(i,ii));
            case 3
                DNA(i,ii)='C';
                num2str(DNA(i,ii));
        end  
    end
end
fprintf("The First generated DNA sequences:\n")
for i=1:numDNA
    if i<10
        fprintf(' %d:::     5''-%s-3''\n',i,DNA(i,:))
    else
        fprintf('%d:::     5''-%s-3''\n',i,DNA(i,:))
    end
end

for i=1:numDNA
DNAInfo(i).Sequence=DNA(i,:);
DNAInfo(i).Generation=1;
DNAInfo(i).Similarity=[];
DNAInfo(i).H_measure=[];
DNAInfo(i).Continuity=[];
DNAInfo(i).GCcontent=[];
DNAInfo(i).Hairpin=[];
DNAInfo(i).Tm=[];
DNAInfo(i).Goal=[];
DNAInfo(i).Transmission=[];
DNAInfo(i).Virulence=Virulence_Determination(i);
DNAInfo(i).Compare_Similarity=[];
DNAInfo(i).Compare_H_measure=[];
end

for i=1:numIndividual
IndividualInfo(i).Virus=0; 
IndividualInfo(i).Type=Init_individual_type();
IndividualInfo(i).Trans_who=[];
IndividualInfo(i).Trans_num=[];
IndividualInfo(i).Trans_distance=[];
IndividualInfo(i).Protect=0;
IndividualInfo(i).Location_X=m(i);
IndividualInfo(i).Location_Y=n(i);
IndividualInfo(i).Die_time=randi([1,floor(Die_Time_max)]);
IndividualInfo(i).Die_probability=rand(1,1);
IndividualInfo(i).Die_time_click=IndividualInfo(i).Die_time;
IndividualInfo(i).Die_will=0;
IndividualInfo(i).T2toT1_time=randi([1,Tran_Time]);
IndividualInfo(i).T2toT1_probability=rand(1,1);
IndividualInfo(i).Cure=randi([1,floor(Cure_Time_max)]);
IndividualInfo(i).OldCureTime=IndividualInfo(i).Cure;
IndividualInfo(i).Infection_time=0;
IndividualInfo(i).R0=0;
IndividualInfo(i).Infstart=1;
end

Ind_Loc_T0=struct('X',[],'Y',[]);
Ind_Loc_T1=struct('X',[],'Y',[]);
Ind_Loc_T2=struct('X',[],'Y',[]);
Ind_Loc_T3=struct('X',[],'Y',[]);
Ind_Loc_T5=struct('X',[],'Y',[]);

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-%
%**********************main*************************%
e=0;
w=waitbar(0,'Start'); 
for i=1:numDNA 
    Virus_List(i)=i;
end
for i=1:numIndividual
    Individual_List(i)=i;
end
Virus_and_Individual(numDNA,numIndividual);
Distance_Table=Transmission_Who_Num_Distance(NowIndividual,Trans_Distance);
Evaluate(DNA,Virus_List);
initype3=Initial_type3(numIndividual);
[Type0x,Type0y,Type1x,Type1y,Type2x,Type2y,Type3x,Type3y,Type5x,Type5y]=Individual_Location_Obatin(numIndividual);
Ind_Loc_T0(1).X=Type0x;
Ind_Loc_T0(1).Y=Type0y;
Ind_Loc_T1(1).X=Type1x;
Ind_Loc_T1(1).Y=Type1y;
Ind_Loc_T2(1).X=Type2x;
Ind_Loc_T2(1).Y=Type2y;
Ind_Loc_T3(1).X=Type3x;
Ind_Loc_T3(1).Y=Type3y;
Ind_Loc_T5(1).X=Type5x;
Ind_Loc_T5(1).Y=Type5y;

Real_Tran_Time=0; 
for i=1:Tran_Time   
    Real_Tran_Time=Real_Tran_Time+1; 
    fprintf('Epidemic Time: %d\n',i)
    if i<Tran_Time
    process=i/(Tran_Time+1);    
    waitbar(i/(Tran_Time+1),w,['Doing ' num2str(100*process) '%']);
    elseif i==Tran_Time
    waitbar(i/(Tran_Time+1),w,'Almost 100%'); 
    end
    diary on;
    for j=Individual_List
        fprintf(".") 
        switch IndividualInfo(j).Type
            case 0 
                
            case 1 
                Die_Determination(IndividualInfo(j).Virus,j,TD); 
                if IndividualInfo(j).Die_will==1 && IndividualInfo(j).Die_time_click>=1 
                    IndividualInfo(j).Die_time_click=IndividualInfo(j).Die_time_click-1; 
                    Evo_New_DNA=Viral_Evolution(DNA,Virus_List,Nc,IndividualInfo(j).Virus);
                    DNA=Evo_New_DNA;
                    
                    if IndividualInfo(j).Infstart==1
                        for k=IndividualInfo(j).Trans_who 
                            f=1; 
                            tran_Pro=Transmission_Probability(IndividualInfo(j).Trans_distance(f)); 
                            [Inf,~,formula]=Infection_Determination(W1,W2,DNAInfo((IndividualInfo(j).Virus)).Transmission,tran_Pro);
                            if Inf==1 && IndividualInfo(k).Type==0 
                                TotalDNA=TotalDNA+1;
                                IndividualInfo(k).Type=I_or_II_Determination(I_probability);
                                IndividualInfo(k).Infstart=0;
                                IndividualInfo(k).Virus=TotalDNA; 
                                IndividualInfo(k).Infection_time=IndividualInfo(k).Infection_time+1;
                                DNAInfo(TotalDNA)=DNAInfo(IndividualInfo(j).Virus);
                                Virus_List(size(Virus_List,2)+1)=IndividualInfo(k).Virus; 
                                if k<j 
                                    DNAInfo(IndividualInfo(k).Virus).Generation=DNAInfo(IndividualInfo(j).Virus).Generation+1;
                                end
                                DNAInfo(IndividualInfo(k).Virus).Virulence=Virulence_Determination(IndividualInfo(k).Virus); 
                                DNA(IndividualInfo(k).Virus,:)=DNAInfo(IndividualInfo(j).Virus).Sequence; 
                            end
                            f=f+1;
                        end
                    end
                elseif IndividualInfo(j).Die_will==1 && IndividualInfo(j).Die_time_click==0 
                    [now_dna,now_individual,individual_List,virus_new_list]=After_Die(j,Individual_List);
                    NowDNA=now_dna;                
                    NowIndividual=now_individual;
                    Virus_List=virus_new_list;
                    Individual_List=individual_List;
                elseif IndividualInfo(j).Die_will==0 && IndividualInfo(j).Cure>=1
                    IndividualInfo(j).Cure=IndividualInfo(j).Cure-1; 
                    Evo_New_DNA=Viral_Evolution(DNA,Virus_List,Nc,IndividualInfo(j).Virus);
                    DNA=Evo_New_DNA;
                    
                    if IndividualInfo(j).Infstart==1
                        for k=IndividualInfo(j).Trans_who 
                            f=1; 
                            tran_Pro=Transmission_Probability(IndividualInfo(j).Trans_distance(f)); 
                            [Inf,~,formula]=Infection_Determination(W1,W2,DNAInfo((IndividualInfo(j).Virus)).Transmission,tran_Pro);
                            if Inf==1 && IndividualInfo(k).Type==0 
                                TotalDNA=TotalDNA+1;
                                IndividualInfo(k).Type=I_or_II_Determination(I_probability);
                                IndividualInfo(k).Infstart=0;
                                IndividualInfo(k).Virus=TotalDNA; 
                                IndividualInfo(k).Infection_time=IndividualInfo(k).Infection_time+1;
                                DNAInfo(TotalDNA)=DNAInfo(IndividualInfo(j).Virus);
                                Virus_List(size(Virus_List,2)+1)=IndividualInfo(k).Virus; 
                                if k<j 
                                    DNAInfo(IndividualInfo(k).Virus).Generation=DNAInfo(IndividualInfo(j).Virus).Generation+1;
                                end
                                DNAInfo(IndividualInfo(k).Virus).Virulence=Virulence_Determination(IndividualInfo(k).Virus); 
                                DNA(IndividualInfo(k).Virus,:)=DNAInfo(IndividualInfo(j).Virus).Sequence; 
                            end
                            f=f+1;
                        end
                    end
                elseif IndividualInfo(j).Die_will==0 && IndividualInfo(j).Cure==0 
                    [now_dna,virus_new_list]=After_Cure(j,Individual_List,Susceptibility); 
                    NowDNA=now_dna; 
                    Virus_List=virus_new_list;
                end
            case 2 
                if IndividualInfo(j).T2toT1_probability<0.5 
                    if IndividualInfo(j).Cure>=1
                        IndividualInfo(j).Cure=IndividualInfo(j).Cure-1;
                        Evo_New_DNA=Viral_Evolution(DNA,Virus_List,Nc,IndividualInfo(j).Virus);
                        DNA=Evo_New_DNA;
                    elseif IndividualInfo(j).Cure==0
                        [now_dna,virus_new_list]=After_Cure(j,Individual_List,Susceptibility); 
                        NowDNA=now_dna; 
                        Virus_List=virus_new_list;
                    end
                elseif IndividualInfo(j).T2toT1_probability>=0.5 
                    if IndividualInfo(j).T2toT1_time>=1 
                        IndividualInfo(j).T2toT1_time=IndividualInfo(j).T2toT1_time-1;
                        Evo_New_DNA=Viral_Evolution(DNA,Virus_List,Nc,IndividualInfo(j).Virus);
                        DNA=Evo_New_DNA;
                    elseif IndividualInfo(j).T2toT1_time==0
                        IndividualInfo(j).Type=1; 
                        IndividualInfo(j).T2toT1_time=randi([1,Tran_Time]);
                    end
                end
            case 3 
                
            case 4 
               
            case 5 
                
        end
    end
    for s=1:numIndividual
        IndividualInfo(s).Infstart=1;
    end
    
    [Type0x,Type0y,Type1x,Type1y,Type2x,Type2y,Type3x,Type3y,Type5x,Type5y]=Individual_Location_Obatin(numIndividual);
    Ind_Loc_T0(i+1).X=Type0x;
    Ind_Loc_T0(i+1).Y=Type0y;
    Ind_Loc_T1(i+1).X=Type1x;
    Ind_Loc_T1(i+1).Y=Type1y;
    Ind_Loc_T2(i+1).X=Type2x;
    Ind_Loc_T2(i+1).Y=Type2y;
    Ind_Loc_T3(i+1).X=Type3x;
    Ind_Loc_T3(i+1).Y=Type3y;
    Ind_Loc_T5(i+1).X=Type5x;
    Ind_Loc_T5(i+1).Y=Type5y;
    
    R0(i+1)=R0_calculate(numIndividual,W1,W2);
    R0(1)=0;
    PVL(i+1)=size(Virus_List,2);
    PVL(1)=numDNA;
    PILN(i+1)=0;
    PILN(1)=numDNA;
    for pilnfor=1:numIndividual
        PILN(i+1)=PILN(i+1)+IndividualInfo(pilnfor).Infection_time;
    end
    fprintf('\n')
    [Imm,Immunity]=Herd_Immunity(NowIndividual,immunobarrier);
    Evaluate(DNA,Virus_List) 
    if Imm==1 
        fprintf('Immunologic barrier£¬Epidemic Time: %d\n',i)
        break 
    end  
end
waitbar(1,w,'Sorting'); 

for j=Virus_List
    Evo_New_DNA=Viral_Evolution(DNA,Virus_List,Nc,j);
    DNA=Evo_New_DNA;
end
EvaluateNewSimHm(DNA,Virus_List) 

[numDie,numI,numII]=Number_Die_I_II(numIndividual);
Main_Parameter_Print(Nc,numDNA,Tran_Time,numIndividual,M_max,N_max,PILN,PVL,R0,Immunity,numDie,initype3,numI,numII,W1,W2,TD,I_probability,Susceptibility,Die_Time_max,Cure_Time_max,immunobarrier,formula);
Fianl_Print(NowDNA,numIndividual) 
fprintf('\n------------------Best result-----------------\n')
Final_nonrepeat_List=DNAInfo_Final(DNA,Virus_List); 
index=Virus_Sort_From_High_to_Low(DNA,Final_nonrepeat_List); 
Save_List=Final_Save(VirusnumSave,index,DNA);
Compare_Final_Print(Save_List,DNA);%As a reference for the paper, a set of sequences was ultimately selected based on the minimum sum of Sim and Hm
Plot_Virus_Infection_R0_Time(PVL,PILN,R0)
Save_List_Table(Save_List)
Infection_Animation(Ind_Loc_T0,Ind_Loc_T1,Ind_Loc_T2,Ind_Loc_T3,Ind_Loc_T5,M_max,N_max,Real_Tran_Time,Time);
waitbar(1,w,'Finished£¡'); 
fprintf('***********************************************************************************\n\n\n\n\n\n\n\n\n')
diary off;
toc
%**********************Finish*************************%
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-%


function Virus_and_Individual(numDNA,numIndividual) 
global IndividualInfo
Ind=randperm(numIndividual);
for i=1:numIndividual
    if IndividualInfo(Ind(i)).Type==0 && numDNA~=0
        IndividualInfo(Ind(i)).Virus=numDNA;
        IndividualInfo(Ind(i)).Type=1;
        IndividualInfo(Ind(i)).Infection_time=1;
        numDNA=numDNA-1;
    end
end
end


function init_individual_type=Init_individual_type()
i=randi([1,10],1,1);
if i<2
   init_individual_type=3; 
else 
   init_individual_type=0;   
end
end


function distable=Transmission_Who_Num_Distance(nowIndividual,Trans_Distance)
global IndividualInfo
distable=zeros(length(nowIndividual),length(nowIndividual));
for i=1:nowIndividual   
    for ii=i:nowIndividual
        nn=sqrt(((IndividualInfo(i).Location_X-IndividualInfo(ii).Location_X)^2+(IndividualInfo(i).Location_Y-IndividualInfo(ii).Location_Y)^2));
        distable(i,ii)=nn;
        distable(ii,i)=distable(i,ii);  
    end
end 
for i=1:nowIndividual 
    w=1; 
    n=1;
    d=1; 
    for ii=1:nowIndividual
        if distable(i,ii)<Trans_Distance && i~=ii 
            IndividualInfo(i).Trans_who(w)=ii;
            IndividualInfo(i).Trans_num=n;
            IndividualInfo(i).Trans_distance(d)=distable(i,ii);
            w=w+1;
            n=n+1;
            d=d+1;
        end
    end
end
end



function Tran_Pro=Transmission_Probability(Real_Distance)
if Real_Distance>=0.5 && Real_Distance<=3
    Tran_Pro=1.2-0.4*Real_Distance;
elseif Real_Distance<0.5
    Tran_Pro=1;
elseif Real_Distance>3
    Tran_Pro=0;
end
end


function [Inf_Det,Inf_prob,formula]=Infection_Determination(W1,W2,DNA_score,T)
Inf_prob=(W1*(tanh(DNA_score/60))+W2*T)/3;
formula='(W1*(tanh(DNA_score/60))+W2*T)/3';
if rand(1,1)<Inf_prob
    Inf_Det=1;
else
    Inf_Det=0;
end
end


function I_or_II=I_or_II_Determination(I_probability)
if rand(1,1)<I_probability
    I_or_II=1;
else 
    I_or_II=2;
end
end


function virulence=Virulence_Determination(i) 
global DNAInfo
if DNAInfo(i).Generation==1
    virulence=100;
elseif DNAInfo(i).Generation>1
    virulence=(100/DNAInfo(i).Generation)+40;
end
end


function Die_Determination(i,j,TD)
global DNAInfo
global IndividualInfo
if DNAInfo(i).Virulence*IndividualInfo(j).Die_probability>=TD
    IndividualInfo(i).Die_will=1;
else
    IndividualInfo(i).Die_will=0;
end
end


function [now_dna,now_individual,individual_List,virus_new_list]=After_Die(i,Individual_List)
global IndividualInfo
    IndividualInfo(i).Type=5;
    IndividualInfo(i).Virus=0;
    Individual_List(find(Individual_List==i))=[];
    list=1;
    Virus_List=[];
    for L=Individual_List
        if IndividualInfo(L).Type==1||IndividualInfo(L).Type==2
            Virus_List(list)=IndividualInfo(L).Virus;
            list=list+1;
        end
    end   

    individual_List=Individual_List;
    virus_new_list=Virus_List;
    now_dna=size(Virus_List,2);
    now_individual=size(Individual_List,2);
end


function [now_dna,virus_new_list]=After_Cure(i,Individual_List,Susceptibility)
global IndividualInfo
if rand(1,1)<Susceptibility
   IndividualInfo(i).Type=0;
   IndividualInfo(i).Cure=IndividualInfo(i).OldCureTime;
else
   IndividualInfo(i).Type=3;
end
IndividualInfo(i).Virus=0;
list=1;
Virus_List=[];
for L=Individual_List
    if IndividualInfo(L).Type==1||IndividualInfo(L).Type==2
        Virus_List(list)=IndividualInfo(L).Virus;
        list=list+1;
    end
end
virus_new_list=Virus_List;
now_dna=size(Virus_List,2);
end


function Final_nonrepeat_List=DNAInfo_Final(dna,Virus_List)
for i=1:size(Virus_List,2)
    List(i)=Virus_List(i);
end
global DNAInfo
F=1;
Final_nonrepeat_List=[];
for i=1:size(List,2)-1
    for j=i+1:size(List,2)
         if string(DNAInfo(List(i)).Sequence)==string(DNAInfo(List(j)).Sequence) && DNAInfo(List(i)).Generation<=DNAInfo(List(j)).Generation && List(i)~=0 && List(j)~=0 
             %List(find(List==List(i)))=0;
             Virus_List(i)=0;
         elseif string(DNAInfo(List(i)).Sequence)==string(DNAInfo(List(j)).Sequence) && DNAInfo(List(i)).Generation>DNAInfo(List(j)).Generation && List(i)~=0 && List(j)~=0 
             %List(find(List==List(j)))=0;
             Virus_List(j)=0;
         end
    end
end 

for i=1:size(Virus_List,2)
    if Virus_List(i)>0
        Final_nonrepeat_List(F)=Virus_List(i);
        F=F+1;
    end
end
Evaluate(dna,Final_nonrepeat_List)
end


function Evo_New_DNA=Viral_Evolution(DNA,Virus_List,Nc,i) 

global DNAInfo
while 1   
    Pal=Palindrome(DNA(i,:));
    if Pal>=7
        changesort=randperm(size(DNA(i,:),2)); 
        DNA(i,:)=DNA(i,changesort);
    else
        break
    end
end
DNA_OriPal=DNA; 
DNA_Original=DNA;
for k=1:Nc
    
    
    hair=Hairpin(DNA(i,:));
    if hair~=0   
        changesort=randperm(size(DNA(i,1:20),2)); 
        DNA(i,1:20)=DNA(i,changesort);
    end
    for bigchange=1:5 
        choose=ceil(rand(1,1)*length(DNA(i,:)));
        if  DNA(i,choose)=='A'
            CH=['T','C']; 
            change=CH(randperm(length(CH),1));
            DNA(i,choose)=change;
        elseif DNA(i,choose)=='T'
            CH=['A','G','C'];
            change=CH(randperm(length(CH),1));
            DNA(i,choose)=change;
        elseif DNA(i,choose)=='G'
            CH=['A','T','C'];
            change=CH(randperm(length(CH),1));
            DNA(i,choose)=change;
        elseif DNA(i,choose)=='C'
            CH=['A','T','G'];
            change=CH(randperm(length(CH),1));
            DNA(i,choose)=change;
        end
    end


    
    DNA_Original=DNA;
    con=continuity(DNA(i,:));
    if con>0
        for j=1:18
            if DNA(i,j)==DNA(i,j+1) && DNA(i,j+1)==DNA(i,j+2) && DNA(i,j)=='A'
                CH=['T','G','C'];
                change=CH(randperm(length(CH),1));
                DNA(i,j+1)=change;
            elseif DNA(i,j)==DNA(i,j+1) && DNA(i,j+1)==DNA(i,j+2) && DNA(i,j)=='T'
                CH=['A','G','C'];
                change=CH(randperm(length(CH),1));
                DNA(i,j+1)=change;
            elseif DNA(i,j)==DNA(i,j+1) && DNA(i,j+1)==DNA(i,j+2) && DNA(i,j)=='G'
                CH=['A','T','C'];
                change=CH(randperm(length(CH),1));
                DNA(i,j+1)=change;
            elseif DNA(i,j)==DNA(i,j+1) && DNA(i,j+1)==DNA(i,j+2) && DNA(i,j)=='C'
                CH=['A','T','G'];
                change=CH(randperm(length(CH),1));
                DNA(i,j+1)=change;
            end
            if continuity(DNA(i,:))>continuity(DNA_Original(i,:))
                DNA=DNA_Original;   
            end
        end
    end
    
    DNA_Original=DNA;
    gc=gccontent(DNA(i,:));
    if gc<0.5
        r=0;
        remark=[];
        remark=repmat(remark,1,20);
        for g=1:20 
            if DNA(i,g)=='A' || DNA(i,g)=='T'
                r=r+1;
                remark(r)=g;
            end
        end
        n=length(remark);
        choose=ceil(rand(1,1)*n);
        CH=['C','G'];
        change=CH(randperm(length(CH),1));
        DNA(i,choose)=change;
    end
    if gc>0.5
        r=0;
        remark=[];
        remark=repmat(remark,1,20);
        for g=1:20 
            if DNA(i,g)=='G' || DNA(i,g)=='C'
                r=r+1;
                remark(r)=g;
            end
        end
        n=length(remark);
        choose=ceil(rand(1,1)*n);
        CH=['A','T'];
        change=CH(randperm(length(CH),1));
        DNA(i,choose)=change;
    end
    if gc==0.5
    end
    Pal=Palindrome(DNA(i,:));
    if Pal>=5
        DNA(i,:)=DNA_OriPal(i,:);
    end
    GTO=goaltotal(DNA_Original,i,Virus_List);
    GT=goaltotal(DNA,i,Virus_List);
    if GTO>GT
        DNA=DNA_Original;
    else
        [sim,~]=similarity(DNA(i,:),DNA(Virus_List,:),Virus_List,i);
        [hme,~]=h_measure(DNA(i,:),DNA(Virus_List,:),Virus_List,i);
        con=continuity(DNA(i,:));
        gc=gccontent(DNA(i,:));
        hair=Hairpin(DNA(i,:));
        [tm,~]=TmBioBox(DNA(i,:));
        DNAInfo(i).Sequence=DNA(i,:);
        DNAInfo(i).Virulence=Virulence_Determination(i); 
        DNAInfo(i).Similarity=sim;
        DNAInfo(i).H_measure=hme;
        DNAInfo(i).Continuity=con;
        DNAInfo(i).GCcontent=gc;
        DNAInfo(i).Hairpin=hair;
        DNAInfo(i).Tm=tm;
        DNAInfo(i).Goal=GT;
        DNAInfo(i).Transmission=GT;
    end
end
Evo_New_DNA=DNA;
end


function [total_dna,now_dna]=Virus_Replication(totalDNA,i,nowDNA,j,I_probability)
global DNAInfo
global IndividualInfo
total_dna=totalDNA+1;
DNAInfo(total_dna)=DNAInfo(i);
IndividualInfo(j).Virus=i;
IndividualInfo(j).Type=I_or_II_Determination(I_probability);
now_dna=nowDNA+1;
end


function [Imm,Immunity]=Herd_Immunity(NowIndividual,immunobarrier)
global IndividualInfo
Immunity=0;
for i=1:NowIndividual
    if IndividualInfo(i).Type==3
        Immunity=Immunity+1;
    end
end
if Immunity/NowIndividual>=immunobarrier
    Imm=1;
else
    Imm=0;
end
end

function [numDie,numI,numII]=Number_Die_I_II(numIndividual)
global IndividualInfo
numDie=0;
numI=0;
numII=0;
for i=1:numIndividual
    if IndividualInfo(i).Type==5
        numDie=numDie+1;
    elseif IndividualInfo(i).Type==1
        numI=numI+1;
    elseif IndividualInfo(i).Type==2
        numII=numII+1;
    end
end
end


function type3=Initial_type3(numIndividual)
type3=0;
global IndividualInfo
for i=1:numIndividual
    if IndividualInfo(i).Type==3
        type3=type3+1;
    end
end
end


function [Type0x,Type0y,Type1x,Type1y,Type2x,Type2y,Type3x,Type3y,Type5x,Type5y]=Individual_Location_Obatin(numIndividual)
global IndividualInfo
Type0x=[];
Type1x=[];
Type2x=[];
Type3x=[];
Type5x=[];
Type0y=[];
Type1y=[];
Type2y=[];
Type3y=[];
Type5y=[];
for i=1:numIndividual
   switch  IndividualInfo(i).Type
       case 0
           Type0x(i)=IndividualInfo(i).Location_X;
           Type0y(i)=IndividualInfo(i).Location_Y;
       case 1
           Type1x(i)=IndividualInfo(i).Location_X;
           Type1y(i)=IndividualInfo(i).Location_Y;
       case 2
           Type2x(i)=IndividualInfo(i).Location_X;
           Type2y(i)=IndividualInfo(i).Location_Y;
       case 3
           Type3x(i)=IndividualInfo(i).Location_X;
           Type3y(i)=IndividualInfo(i).Location_Y;
       case 5
           Type5x(i)=IndividualInfo(i).Location_X;
           Type5y(i)=IndividualInfo(i).Location_Y;
   end
end
Type0x=Type0x(find(Type0x~=0));
Type0y=Type0y(find(Type0y~=0));
Type1x=Type1x(find(Type1x~=0));
Type1y=Type1y(find(Type1y~=0));
Type2x=Type2x(find(Type2x~=0));
Type2y=Type2y(find(Type2y~=0));
Type3x=Type3x(find(Type3x~=0));
Type3y=Type3y(find(Type3y~=0));
Type5x=Type5x(find(Type5x~=0));
Type5y=Type5y(find(Type5y~=0));
end


%---------------------------------------%
function Main_Parameter_Print(Nc,numDNA,Tran_Time,numIndividual,M_max,N_max,PILN,PVL,R0,Immunity,numDie,initype3,numI,numII,W1,W2,TD,I_probability,Susceptibility,Die_Time_max,Cure_Time_max,immunobarrier,formula)
fprintf('\n Main parameters and result information \n')
fprintf('The number of times the virus evolves in the body: %d\n',Nc)
fprintf('Initial DNA (virus): %d\n',numDNA)
fprintf('Initial host: %d\n',numIndividual)
fprintf('Biggest epidemic Time: %d\n',Tran_Time)
fprintf('DNA score infection weight coefficient W1: %.1f\n',W1)
fprintf('Transmission event infection weight coefficient W2: %.1f\n',W2)
fprintf('Formula for infection detection: %s\n',formula)
fprintf('Death threshold TD: %d\n',TD)
fprintf('The probability of being a type I infected person after being infected_ Probability: %.1f\n',I_probability)
fprintf('After treatment, the probability of being a susceptible type: %.1f\n',Susceptibility)
fprintf('The maximum possible time of death for each infected individual (each individual is different) Die_ Time_ Max: %.1f\n',Die_Time_max)
fprintf('The maximum possible self-healing time Cure for each infected individual (each individual is different)_ Time_ Max: %.1f\n',Cure_Time_max)
fprintf('Immunobarrier index: %.1f\n',immunobarrier)

fprintf('Space length Lm: %d\n',M_max)
fprintf('Space width Ln: %d\n',N_max)
fprintf('Per capita land area: %.2f\n',(M_max*N_max)/numIndividual)
fprintf('Area of transmission for type I infected individuals: %.2f\n',3*3*3.14)
fprintf('The average number of people within the transmission range of one type I infected person: %.2f\n',(3*3*3.14*numIndividual)/(M_max*N_max))
fprintf('Number of natural immune types: %d   %.2f\n',initype3,initype3/numIndividual)
fprintf('Number of final immunization types: %d   %.2f\n',Immunity,Immunity/numIndividual)
fprintf('Final death toll: %d    %.2f\n',numDie,numDie/numIndividual)
fprintf('Remaining number of Type I infections: %d    %.2f\n',numI,numI/numIndividual)
fprintf('Remaining Type II infections: %d    %.2f\n',numII,numII/numIndividual)
fprintf('Accumulated number of infections over time:\n')
PILN
fprintf('The current number of infected individuals over time: \n')
PVL
fprintf('R0: \n')
R0
fprintf('\n')
end


%---------------------------------------%
function Fianl_Print(NowDNA,numIndividual)
global DNAInfo
global IndividualInfo
for i=1:numIndividual
    switch IndividualInfo(i).Type
        case 0
            fprintf('\n')
            fprintf('The final situation of the %d individual:\n',i)
            fprintf('Individuals are susceptible types\n')
            fprintf('Location: [%.2f,%.2f]\n',IndividualInfo(i).Location_X,IndividualInfo(i).Location_Y)
            fprintf('Probability of death and time required for death: %.2f£¬%d\n',IndividualInfo(i).Die_probability,IndividualInfo(i).Die_time)
            fprintf('Probability of Type II to Type I conversion and required time: %.2f£¬%d\n',IndividualInfo(i).T2toT1_probability,IndividualInfo(i).T2toT1_time)
            fprintf('Time required for cure: %d\n',IndividualInfo(i).Cure)
        case 1
            fprintf('\n')
            fprintf('The final situation of the %d individual:\n',i)
            fprintf('Individual is a type I carrier\n')
            fprintf('Location: [%.2f,%.2f]\n',IndividualInfo(i).Location_X,IndividualInfo(i).Location_Y)
            fprintf('Probability of death and time required for death: %.2f£¬%d\n',IndividualInfo(i).Die_probability,IndividualInfo(i).Die_time)
            fprintf('Probability of Type II to Type I conversion and required time: %.2f£¬%d\n',IndividualInfo(i).T2toT1_probability,IndividualInfo(i).T2toT1_time)
            fprintf('Time required for cure: %d\n',IndividualInfo(i).Cure)
            fprintf('Carrying viruses: %d\n',IndividualInfo(i).Virus)
            fprintf('Sequence: %s\n',DNAInfo(IndividualInfo(i).Virus).Sequence)
            fprintf('Generation=%d\n',DNAInfo(IndividualInfo(i).Virus).Generation)
            fprintf('Similarity=%.2f\n',DNAInfo(IndividualInfo(i).Virus).Similarity/(NowDNA-1))
            fprintf('H_measure=%.2f\n',DNAInfo(IndividualInfo(i).Virus).H_measure/NowDNA)
            fprintf('Continuity=%d\n',DNAInfo(IndividualInfo(i).Virus).Continuity)
            fprintf('GC_content=%.2f\n',DNAInfo(IndividualInfo(i).Virus).GCcontent)
            fprintf('Tm=%.2f\n',DNAInfo(IndividualInfo(i).Virus).Tm)
            fprintf('Hairpinnum=%d\n',DNAInfo(IndividualInfo(i).Virus).Hairpin)
            fprintf('Transmission=%.2f\n',DNAInfo(IndividualInfo(i).Virus).Transmission)
            fprintf('Virulence=%.2f\n',DNAInfo(IndividualInfo(i).Virus).Virulence)
        case 2
            fprintf('\n')
            fprintf('The final situation of the %d individual:\n',i)
            fprintf('Individual is a type II carrier\n')
            fprintf('Location: [%.2f,%.2f]\n',IndividualInfo(i).Location_X,IndividualInfo(i).Location_Y)
            fprintf('Probability of death and time required for death: %.2f£¬%d\n',IndividualInfo(i).Die_probability,IndividualInfo(i).Die_time)
            fprintf('Probability of Type II to Type I conversion and required time: %.2f£¬%d\n',IndividualInfo(i).T2toT1_probability,IndividualInfo(i).T2toT1_time)
            fprintf('Time required for cure: %d\n',IndividualInfo(i).Cure)
            fprintf('Carrying viruses: %d\n',IndividualInfo(i).Virus)
            fprintf('Sequence: %s\n',DNAInfo(IndividualInfo(i).Virus).Sequence)
            fprintf('Generation=%d\n',DNAInfo(IndividualInfo(i).Virus).Generation)
            fprintf('Similarity=%.2f\n',DNAInfo(IndividualInfo(i).Virus).Similarity/(NowDNA-1))
            fprintf('H_measure=%.2f\n',DNAInfo(IndividualInfo(i).Virus).H_measure/NowDNA)
            fprintf('Continuity=%d\n',DNAInfo(IndividualInfo(i).Virus).Continuity)
            fprintf('GC_content=%.2f\n',DNAInfo(IndividualInfo(i).Virus).GCcontent)
            fprintf('Tm=%.2f\n',DNAInfo(IndividualInfo(i).Virus).Tm)
            fprintf('Hairpinnum=%d\n',DNAInfo(IndividualInfo(i).Virus).Hairpin)
            fprintf('Transmission=%.2f\n',DNAInfo(IndividualInfo(i).Virus).Transmission)
            fprintf('Virulence=%.2f\n',DNAInfo(IndividualInfo(i).Virus).Virulence)
        case 3
            fprintf('\n')
            fprintf('The final situation of the %d individual:\n',i)
            fprintf('Individuals have immune capacity\n')
            fprintf('Location: [%.2f,%.2f]\n',IndividualInfo(i).Location_X,IndividualInfo(i).Location_Y)
        case 4
            fprintf('\n')
            fprintf('The final situation of the %d individual:\n',i)
            fprintf('The individual has just been cured\n')
            fprintf('Location: [%.2f,%.2f]\n',IndividualInfo(i).Location_X,IndividualInfo(i).Location_Y)
            fprintf('Probability of death and time required for death: %.2f£¬%d\n',IndividualInfo(i).Die_probability,IndividualInfo(i).Die_time)
            fprintf('Probability of Type II to Type I conversion and required time: %.2f£¬%d\n',IndividualInfo(i).T2toT1_probability,IndividualInfo(i).T2toT1_time)
            fprintf('Time required for cure: %d\n',IndividualInfo(i).Cure)
        case 5
            fprintf('\n')
            fprintf('The final situation of the %d individual:\n',i)
            fprintf('Individual has died \n')
            fprintf('Location: [%.2f,%.2f]\n',IndividualInfo(i).Location_X,IndividualInfo(i).Location_Y)
    end
end   
end


function index=Virus_Sort_From_High_to_Low(dna,Virus_List_nonrepeat)
global DNAInfo
Evaluate(dna,Virus_List_nonrepeat)
for i=1:size(Virus_List_nonrepeat,2)
    V(i)=DNAInfo(Virus_List_nonrepeat(i)).Transmission;
end
[~,des]=sort(V,'descend');
for i=1:size(Virus_List_nonrepeat,2)
    index(i)=Virus_List_nonrepeat(des(i));
end
fprintf('The viruses are ranked in descending order of transmission power:\n')
for j=index
    fprintf('\n')
    fprintf('Virus number: %d\n',j)
    fprintf('Sequence: %s\n',DNAInfo(j).Sequence)
    fprintf('Generation=%d\n',DNAInfo(j).Generation)
    fprintf('Similarity=%.2f\n',DNAInfo(j).Similarity/(size(Virus_List_nonrepeat,2)-1))
    fprintf('H_measure=%.2f\n',DNAInfo(j).H_measure/size(Virus_List_nonrepeat,2))
    fprintf('Continuity=%d\n',DNAInfo(j).Continuity)
    fprintf('GC_content=%.2f\n',DNAInfo(j).GCcontent)
    fprintf('Tm=%.2f\n',DNAInfo(j).Tm)
    fprintf('Hairpinnum=%d\n',DNAInfo(j).Hairpin)
    fprintf('Transmission=%.2f\n',DNAInfo(j).Transmission)
    fprintf('Virulence=%.2f\n',DNAInfo(j).Virulence)
end
end

function Save_List=Final_Save(VirusnumSave,index,dna)
global DNAInfo
if size(index,2)>=10
    shortindex=index(1:10); 
else
    shortindex=index;
end
FinalsaveV=nchoosek(shortindex,VirusnumSave);
[Vsize,~]=size(FinalsaveV);
Simmax=0;
Hmmax=0;
for i=1:Vsize
    Smax=0;
    Hmax=0;
    for j=1:VirusnumSave  
    [Hm,Sim]=HmSm(DNAInfo(FinalsaveV(i,j)).Sequence,dna(FinalsaveV(i,1:VirusnumSave),:),FinalsaveV(i,1:VirusnumSave),FinalsaveV(i,j));
    Smax=Smax+Sim;
    Hmax=Hmax+Hm;
    end
    if i==1
        Simmax=Smax;
        Hmmax=Hmax;
        Totalmax=Simmax+Hmmax;
        List=FinalsaveV(i,1:VirusnumSave);
    elseif i>1 && Totalmax>Smax+Hmax 
        Totalmax=Simmax+Hmmax;
        List=FinalsaveV(i,1:VirusnumSave);
    end
end
Save_List=List;
end

function Compare_Final_Print(Save_List,dna)
global DNAInfo
for i=1:size(Save_List,2)
    [Hm,Sim]=HmSm(DNAInfo(Save_List(i)).Sequence,dna(Save_List,:),Save_List,Save_List(i));
    DNAInfo(Save_List(i)).Compare_Similarity=Sim;
    DNAInfo(Save_List(i)).Compare_H_measure=Hm;
end
fprintf('\n A set of sequences selected based on the minimum sum of Sim and Hm:')
aveGen=0;
aveSim=0;
aveHm=0;
aveCon=0;
aveGC=0;
aveTm=0;
var=zeros(1,size(Save_List,2));
varTm=0;
i=1;
aveHair=0;
for j=Save_List
    fprintf('\n')
    fprintf('Virus number: %d\n',j)
    fprintf('Sequence: %s\n',DNAInfo(j).Sequence)
    fprintf('Generation=%d\n',DNAInfo(j).Generation)
    fprintf('Compare_Similarity=%d\n',DNAInfo(j).Compare_Similarity)
    fprintf('Compare_H_measure=%d\n',DNAInfo(j).Compare_H_measure)
    fprintf('Continuity=%d\n',DNAInfo(j).Continuity)
    fprintf('GC_content=%.2f\n',DNAInfo(j).GCcontent)
    fprintf('Tm=%.2f\n',DNAInfo(j).Tm)
    fprintf('Hairpinnum=%d\n',DNAInfo(j).Hairpin)
    fprintf('Transmission=%.2f\n',DNAInfo(j).Transmission)
    fprintf('Virulence=%.2f\n',DNAInfo(j).Virulence)
    aveGen=aveGen+DNAInfo(j).Generation;
    aveSim=aveSim+DNAInfo(j).Compare_Similarity;
    aveHm=aveHm+DNAInfo(j).Compare_H_measure;
    aveCon=aveCon+DNAInfo(j).Continuity;
    aveGC=aveGC+DNAInfo(j).GCcontent;
    aveTm=aveTm+DNAInfo(j).Tm;
    aveHair=aveHair+DNAInfo(j).Hairpin;
    var(i)=DNAInfo(j).Tm;
    i=i+1;
end
varTm=sum((var(1,:)-mean(var)).^2)/length(var);
fprintf('-----------------\n')
fprintf('AVGGen=%.2f\n',aveGen/size(Save_List,2))
fprintf('AVGSim=%.2f\n',aveSim/size(Save_List,2))
fprintf('AVGHm=%.2f\n',aveHm/size(Save_List,2))
fprintf('AVGCon=%.2f\n',aveCon/size(Save_List,2))
fprintf('AVGGC=%.2f\n',aveGC/size(Save_List,2))
fprintf('AVGTm=%.2f     VARTm=%.4f\n',aveTm/size(Save_List,2),varTm)
fprintf('AVGHair=%.2f\n',aveHair/size(Save_List,2))

end

function R0=R0_calculate(numIndividual,W1,W2)
global IndividualInfo
global DNAInfo
R0=[];
R0_total=0;
infnum=0; 
for i=1:numIndividual
    if IndividualInfo(i).Type==1 || IndividualInfo(i).Type==2 
        infnum=infnum+1;
    end
    A=0;
    if IndividualInfo(i).Type==1 && IndividualInfo(i).Die_will==1
        A=IndividualInfo(i).Die_time; 
    elseif IndividualInfo(i).Type==1 && IndividualInfo(i).Die_will==0
        A=IndividualInfo(i).Cure; 
    end
    B=0;
    C=0;
    f=1;
    if IndividualInfo(i).Type==1
    for j=IndividualInfo(i).Trans_who 
        tran_Pro=Transmission_Probability(IndividualInfo(i).Trans_distance(f)); 
        [~,Inf_prob,~]=Infection_Determination(W1,W2,DNAInfo((IndividualInfo(i).Virus)).Transmission,tran_Pro);
        if IndividualInfo(j).Type==0
            B=B+1;  
            C=C+Inf_prob;
        end
        f=f+1;
    end
    end
    R0_i=A*B*C;
    R0_total=R0_total+R0_i;
end
R0=R0_total/infnum;
end

function Plot_Virus_Infection_R0_Time(PVL,PILN,R0)
figure
subplot(1,3,1)
x1=0:1:size(PVL,2)-1;
plot(x1,PVL)
xlabel('Epidemic Time')
ylabel('Current infection number')
title('Current number of infected individuals and epidemic time curve')
ylim([min(PVL)-10 max(PVL)+5])
xlim([0 size(PVL,2)+1])
text(x1,PVL,'o');

subplot(1,3,2)
x2=0:1:size(PILN,2)-1;
plot(x2,PILN)
xlabel('Epidemic Time')
ylabel('Accumulated number of infections')
title('Accumulated number of infections and epidemic time curve')
ylim([min(PILN)-10 max(PILN)+5])
xlim([0 size(PILN,2)+1])
text(x2,PILN,'*');

subplot(1,3,3)
x3=0:1:size(R0,2)-1;
plot(x3,R0)
xlabel('Epidemic Time')
ylabel('R0')
title('R0 and time curve')
ylim([0 max(R0)+2])
xlim([0 size(R0,2)+1])
text(x3,R0,'¡÷');
set(gcf, 'unit', 'centimeters', 'position', [12 1 35 12]);
end


function Save_List_Table(Save_List)
global DNAInfo
data=zeros(8,9);
Gen=0;
Comsim=0;
Comhm=0;
Ccon=0;
Cgc=0;
Chair=0;
TTM=0;
var=zeros(1,size(Save_List,2));
Ctran=0;
Cvir=0;
for i=1:size(Save_List,2)
    data(i,:)=[DNAInfo(Save_List(i)).Generation,DNAInfo(Save_List(i)).Compare_Similarity,DNAInfo(Save_List(i)).Compare_H_measure,DNAInfo(Save_List(i)).Continuity,DNAInfo(Save_List(i)).GCcontent,DNAInfo(Save_List(i)).Hairpin,DNAInfo(Save_List(i)).Tm,DNAInfo(Save_List(i)).Transmission,DNAInfo(Save_List(i)).Virulence];
    Gen=Gen+DNAInfo(Save_List(i)).Generation;
    Comsim=Comsim+DNAInfo(Save_List(i)).Compare_Similarity;
    Comhm=Comhm+DNAInfo(Save_List(i)).Compare_H_measure;
    Ccon=Ccon+DNAInfo(Save_List(i)).Continuity;
    Cgc=Cgc+DNAInfo(Save_List(i)).GCcontent;
    Chair=Chair+DNAInfo(Save_List(i)).Hairpin;
    TTM=TTM+DNAInfo(Save_List(i)).Tm;
    var(i)=DNAInfo(Save_List(i)).Tm;
    Ctran=Ctran+DNAInfo(Save_List(i)).Transmission;
    Cvir=Cvir+DNAInfo(Save_List(i)).Virulence;
end
varTm=sum((var(1,:)-mean(var)).^2)/length(var);
data(8,:)=[Gen/7,Comsim/7,Comhm/7,Ccon/7,Cgc/7,Chair/7,TTM/7,Ctran/7,Cvir/7];
data(9,:)=[0,0,0,0,0,0,varTm,0,0];
column_name={'Generation','Compare_Similarity','Compare_H_measure','Continuity','GCcontent','Hairpin','Tm','Transmission','Virulence'};
row_name={DNAInfo(Save_List(1)).Sequence,DNAInfo(Save_List(2)).Sequence,DNAInfo(Save_List(3)).Sequence,DNAInfo(Save_List(4)).Sequence,DNAInfo(Save_List(5)).Sequence,DNAInfo(Save_List(6)).Sequence,DNAInfo(Save_List(7)).Sequence,'Average','Variance'};
f=figure;
uitable(f,'Data',data,'Position',[20 20 1100 200],'Columnname',column_name,'Rowname',row_name);
set(gcf, 'unit', 'centimeters', 'position', [5 3 30 7]);
end



function Infection_Animation(Ind_Loc_T0,Ind_Loc_T1,Ind_Loc_T2,Ind_Loc_T3,Ind_Loc_T5,M_max,N_max,Real_Tran_Time,nowTime)
figure
filename=['E:\DNA_design_Virus_Static\' num2str(nowTime(1)),'_',num2str(nowTime(2)),'_',num2str(nowTime(3)),'_',num2str(nowTime(4)),'_',num2str(nowTime(5)) '.gif'];
for i=1:Real_Tran_Time+1
    if i==1
        hold on
        scatter(Ind_Loc_T0(i).X,Ind_Loc_T0(i).Y,'filled','b');
        scatter(Ind_Loc_T1(i).X,Ind_Loc_T1(i).Y,'filled','r');
        scatter(Ind_Loc_T2(i).X,Ind_Loc_T2(i).Y,'filled','m');
        scatter(Ind_Loc_T3(i).X,Ind_Loc_T3(i).Y,'filled','g');
        scatter(Ind_Loc_T5(i).X,Ind_Loc_T5(i).Y,'filled','k');
        xlabel('The initial stage of virus attack')
        le=legend(['Susceptibility: ',num2str(size(Ind_Loc_T0(i).X,2))],['Type I infection: ',num2str(size(Ind_Loc_T1(i).X,2))],['Type II infection: ',num2str(size(Ind_Loc_T2(i).X,2))],['Immune: ',num2str(size(Ind_Loc_T3(i).X,2))],['Die: ',num2str(size(Ind_Loc_T5(i).X,2))]);
        set(le,'Location','southoutside');
        set(le,'Orientation','horizon','Box','off');
        hold off
    else
        hold on
        scatter(Ind_Loc_T0(i).X,Ind_Loc_T0(i).Y,'filled','b');
        scatter(Ind_Loc_T1(i).X,Ind_Loc_T1(i).Y,'filled','r');
        scatter(Ind_Loc_T2(i).X,Ind_Loc_T2(i).Y,'filled','m');
        scatter(Ind_Loc_T3(i).X,Ind_Loc_T3(i).Y,'filled','g');
        scatter(Ind_Loc_T5(i).X,Ind_Loc_T5(i).Y,'filled','k');
        le=legend(['Susceptibility: ',num2str(size(Ind_Loc_T0(i).X,2))],['Type I infection: ',num2str(size(Ind_Loc_T1(i).X,2))],['Type II infection: ',num2str(size(Ind_Loc_T2(i).X,2))],['Immune: ',num2str(size(Ind_Loc_T3(i).X,2))],['Die: ',num2str(size(Ind_Loc_T5(i).X,2))]);
        set(le,'Location','southoutside');
        set(le,'Orientation','horizon','Box','off');
        xlabel(['Virus epidemic time£º',num2str(i-1)])
        hold off
    end
    ylim([0 M_max])
    xlim([0 N_max])
    set(gcf, 'unit', 'centimeters', 'position', [1 1 20 20]);
    MakeGif(filename,i);
    pause(2.5);
    drawnow;
end
end

function MakeGif(filename,index)  
    f = getframe(gcf);  
    imind = frame2im(f);  
    [imind,cm] = rgb2ind(imind,256);  
    if index==1  
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',2.5);
    else  
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',2.5);
    end  
end  



function Pal=Palindrome(X)
Xp=reverse(X); 
for c=1:20 
    switch Xp(c)
        case 'A'
            Xp(c)='T';
        case 'T'
            Xp(c)='A';
        case 'G'
            Xp(c)='C';
        case 'C'
            Xp(c)='G';
    end
end
s=zeros(1,40);
for k=1:19  
    x=X;
    y=circshift(Xp,k);
    for i=k:20
        if  x(i)==y(i)
            s(k)=s(k)+1;
        end
    end
end
for kk=1:19 
    y=Xp;
    x=circshift(X,kk);
    for i=kk+1:20
        if  x(i)==y(i)
            s(19+kk)=s(kk+19)+1;
        end
    end
end
for i=1:20 
    if  X(i)==Xp(i)
        s(40)=s(40)+1;
    end
end
[ss,~]=max(s);
Pal=ss;
end


function [S,Smax]=similarity(xo,yo,Virus_List,flag)
F=find(Virus_List==flag);
s=zeros(1,40);
smax=zeros(1,size(Virus_List,2));
stotal=0;
Smax=0;
for iis=1:size(Virus_List,2)
    if F~=iis
        for k=1:19  
            x=xo;
            y=circshift(yo(iis,:),k);
            for i=k:20
                if  x(i)==y(i)
                    s(k)=s(k)+1;
                end
            end
        end
        for kk=1:19 
            y=yo;
            x=circshift(xo,kk);
            for i=kk+1:20
                if  x(i)==y(iis,i)
                    s(19+kk)=s(kk+19)+1;
                end
            end
        end
        for i=1:20
            if  xo(i)==yo(iis,i)
                s(40)=s(40)+1;
            end
        end
        [ss,~]=max(s);
        s(:,:)=0;
        stotal=stotal+ss;
        smax(iis)=ss;
    end
end
[Smax,~]=max(smax);
S=stotal;
end

function [hm,Hmax]=h_measure(xo,yo,Virus_List,flag)
F=find(Virus_List==flag);
s=zeros(1,40);
hmax=zeros(1,size(Virus_List,2));
htotal=0;

xo=reverse(xo); 
for c=1:20 
    switch xo(c)
        case 'A'
            xo(c)='T';
        case 'T'
            xo(c)='A';
        case 'G'
            xo(c)='C';
        case 'C'
            xo(c)='G';
    end
end

for ii=1:size(Virus_List,2)
    if F~=ii
        for k=1:19  
            x=xo;
            y=circshift(yo(ii,:),k);
            for i=k:20
                if  x(i)==y(i)
                    s(k)=s(k)+1;
                end
            end
        end
        for kk=1:19 
            y=yo;
            x=circshift(xo,kk);
            for i=kk+1:20
                if  x(i)==y(ii,i)
                    s(19+kk)=s(kk+19)+1;
                end
            end
        end
        for i=1:20
            if  xo(i)==yo(ii,i)
                s(40)=s(40)+1;
            end
        end
        [ss,~]=max(s);
        s(:,:)=0;
        htotal=htotal+ss;
        hmax(ii)=ss;
    end
end
[Hmax,~]=max(hmax);
hm=htotal;
end

function [Hmflag,Simflag]=HmSm(xo,yo,Virus_List,flag)
F=find(Virus_List==flag);
[m,l]=size(yo);
Gap1=round(l/4);  
DNA1=xo;
ReverseDNA1=seqreverse(DNA1);
for i=1:20
    switch ReverseDNA1(i)
        case 'A'
            ReverseDNA1(i)='T';
        case 'T'
            ReverseDNA1(i)='A';
        case 'C'
            ReverseDNA1(i)='G';
        case 'G'
            ReverseDNA1(i)='C';
    end
end
SimValue=zeros(m,1);
HmValue=zeros(m,1);
for j=1:size(Virus_List,2)

    for g=0:Gap1  
        tempIntDNAy=[yo(j,:), zeros(1,g)+5, yo(j,:)];
        for i=-l+1:l-1
            ShiftValue=shift(tempIntDNAy,i);
            currentSim=s_dis(DNA1,ShiftValue)+s_con(DNA1,ShiftValue);
            currentHm=s_dis(ReverseDNA1,ShiftValue)+s_con(ReverseDNA1,ShiftValue);
            SimValue(j)=max(SimValue(j),currentSim);
            HmValue(j)=max(HmValue(j),currentHm);
        end
    end
end
SimValue(F)=0;

Hmflag=sum(HmValue);
Simflag=sum(SimValue);
end


function IntDNAShifted= shift( IntDNA_nb,i )
if i==0
    IntDNAShifted=IntDNA_nb;
    return;
end
l=size(IntDNA_nb,2);
temp=zeros(1,abs(i))+5;
if (0<i) && (i<l)   
    IntDNAShifted=[temp, IntDNA_nb(1:l-i)];
    return;
end
if (i<0) && (i>-l) 
    IntDNAShifted=[IntDNA_nb(abs(i)+1:l),temp];
    return;
end
if abs(i)>=l 
    IntDNAShifted=zeros(1,l)+5;
    return;
end
end

function IntValue= s_con( IntDNAx,IntDNAy,S_con)

S_CON = 6;
if nargin==2
    S_con=S_CON;
end
IntValue=0;
l=size(IntDNAx,2);
for i=1:l
    IntValue=IntValue+T(ceq(IntDNAx,IntDNAy,i),S_con);
end
end

function IntValue = s_dis(IntDNAx,IntDNAy,S_dis)

S_DIS = 0.17;
if nargin==2
    S_dis=S_DIS;
end
Sigma_eq=0;
l=size(IntDNAx,2);
for i=1:l
    Sigma_eq=Sigma_eq+eqq(IntDNAx(i),IntDNAy(i));
end
temp=S_dis*length_nb(IntDNAy)/2;
IntValue=T(Sigma_eq,temp);
end

function LenNoBlank= length_nb(IntDNA_nb)
LenNoBlank=sum(IntDNA_nb~=5);
end

function IntValue= T(i,j)
if i>j 
    IntValue=i;
else
    IntValue=0;
end
end

function BoolValue= eqq(IntDNAa,IntDNAb)
if size(IntDNAa,2)==size(IntDNAb,2)
    if sum(IntDNAa==IntDNAb)==size(IntDNAa,2)
        BoolValue=1;
    else
        BoolValue=0;
    end
else
    BoolValue=0;
    error('DNA different length');    
end
end

function IntValue = ceq( IntDNAx,IntDNAy,i)
IntValue=0;
l=size(IntDNAx,2);
if i>l
    error('i OUT EDGE');
end
if i~=1
    if eqq(IntDNAx(i),IntDNAy(i))==0
        j=1;
        while j<=l-i && eqq(IntDNAx(i+j),IntDNAy(i+j))
            IntValue=IntValue+1;
            j=j+1;
        end
    else
        return;
    end
else %i==1
    if eqq(IntDNAx(i),IntDNAy(i))==0
        j=1;
        while j<=l-i && eqq(IntDNAx(i+j),IntDNAy(i+j))
            IntValue=IntValue+1;
            j=j+1;
        end        
    else
        j=0;
        while j<=l-i && eqq(IntDNAx(i+j),IntDNAy(i+j))
            IntValue=IntValue+1;
            j=j+1;
        end
    end    
end
end



function con=continuity(dna)
c=zeros(1,6);
j=1;
dna(21)='0';
c1=0;
c2=0;
for i=1:18  
    if dna(i)==dna(i+1) && dna(i+1)==dna(i+2) 
        c(j)=c(j)+1;
        if dna(i)~=dna(i+3)
            j=j+1;
        end
    end
end
if c(1)~=0
    for i=1:size(c,2)
        c1=c1+c(i);
        if c(i)~=0
            c2=c2+1;
        end
    end
    c1=c1+2*c2;
end
con=c1;
end


function gc=gccontent(x)
p=0;
for i=1:20
    if x(i)=='G' || x(i)=='C'
        p=p+1;
    end
end
gc=p/20;
end

function [Tm,GC] = TmBioBox(x)
SaltValue=1; 
PrimerconcValue=10^(-8); 
m=size(x,1);
GC=zeros(m,1);
Tm=GC;
for i=1:m
    oligoprops=oligoprop(x(i,:),'Salt', SaltValue,'Primerconc', PrimerconcValue);
    Tm(i)=round(oligoprops.Tm(5)*10000)/10000; 
end
end

function Values = Hairpin( DNAs )
R_MIN=6;
P_MIN=6;
PINLEN=3;
if nargin==1
    R_min=R_MIN;
    P_min=P_MIN;
end
[m,l]=size(DNAs);
Values=zeros(m,1);
for k=1:m
    IntDNA=DNAs(k,:);
    for p=P_min:(l-R_min)/2
        for r=R_min:(l-2*p)
            for i=0:(l-2*p-r)
                Sigma_bp=0;
                for j=0:pinlen(p,r,i,IntDNA)-1
                    Sigma_bp=Sigma_bp+bp(IntDNA(p+i-j),IntDNA(p+i+r+1+j));
                end
                Values(k)=Values(k)+T(Sigma_bp,pinlen(p,r,i,IntDNA)/ PINLEN);
            end
        end
    end
end
end

function minValue = pinlen( p,r,i,IntDNA )
l=size(IntDNA,2);
minValue=min(p+i,l-p-i-r);
end

function BoolValue = bp( IntDNAa,IntDNAb )
l=size(IntDNAa,2);
m=size(IntDNAb,2);
if l==m
    temp=3-IntDNAa;
    if sum(temp==IntDNAb)==l
        BoolValue=1;
    else
        BoolValue=0;
    end
else
    BoolValue=0;
    error('DNA different length');
end
end



function Evaluate(dna,Virus_list)
goal=[];
goal=repmat(goal,1,size(Virus_list,2));
global DNAInfo
for i=Virus_list
    goal(i)=0;
end
for i=Virus_list
    [sim,~]=similarity(dna(i,:),dna(Virus_list,:),Virus_list,i);
    [hme,~]=h_measure(dna(i,:),dna(Virus_list,:),Virus_list,i);
    con=continuity(dna(i,:));
    Tm=TmBioBox(dna(i,:));
    gc=gccontent(dna(i,:));
    hair=Hairpin(dna(i,:));
    if con==0
        goal(i)=goal(i)+20;
    elseif con==1
        goal(i)=goal(i)+10;
    elseif con==2
        goal(i)=goal(i)+5;
    else
        goal(i)=goal(i)+0;
    end
    if Tm<=59
        goal(i)=goal(i)+40;
    elseif Tm>59 && Tm<=69
        goal(i)=goal(i)+((69-Tm)/10)*40;
    else
        goal(i)=goal(i)+0;
    end
    
    %If you want to obtain low Tm variance sequences, replace program
    %Tm=zeros(1,length(Virus_list));
    %for ii=Virus_list
    %    [TM,~]=TmBioBox(dna(ii,:));
    %    Tm(ii)=TM;
    %end
    %VARTm=var(Tm);
    %if VARTm<=1
    %    goal(i)=goal(i)+40;
    %elseif VARTm>1 && VARTm<=3
    %    goal(i)=goal(i)+((3-VARTm)/2)*40;
    %else
    %    goal(i)=goal(i)+0;
    %end
    
   
    if sim/(size(Virus_list,2)-1)<=7
        goal(i)=goal(i)+30;
    elseif sim/(size(Virus_list,2)-1)>7 && sim/(size(Virus_list,2)-1)<13
        goal(i)=goal(i)+(13-sim/(size(Virus_list,2)-1))*5;
    elseif sim/(size(Virus_list,2)-1)>=13 
        goal(i)=goal(i)+0;
    end
    if hme/(size(Virus_list,2)-1)<=7
        goal(i)=goal(i)+30;
    elseif hme/(size(Virus_list,2)-1)>7 && hme/(size(Virus_list,2)-1)<13
        goal(i)=goal(i)+(13-hme/size(Virus_list,2))*5;
    elseif hme/(size(Virus_list,2)-1)>=13 || sim/(size(Virus_list,2)-1)>=13 
        goal(i)=goal(i)+0;
    end
    if gc==0.5 && hair==0
        goal(i)=goal(i);
    else
        goal(i)=0;
    end
    DNAInfo(i).Sequence=dna(i,:);
    DNAInfo(i).Similarity=sim;
    DNAInfo(i).H_measure=hme;
    DNAInfo(i).Continuity=con;
    DNAInfo(i).GCcontent=gc;
    DNAInfo(i).Hairpin=hair;
    DNAInfo(i).Tm=Tm;
    DNAInfo(i).Goal=goal(i);
    DNAInfo(i).Transmission=DNAInfo(i).Goal;
end
end

function GTotal=goaltotal(dna,i,Virus_list)
goal=0;
goal=double(goal);
if continuity(dna(i,:))==0
    goal=goal+20;
elseif continuity(dna(i,:))==1
    goal=goal+10;
elseif continuity(dna(i,:))==2
    goal=goal+5;
else
    goal=goal+0;
end
if TmBioBox(dna(i,:))<=59
    goal=goal+40;
elseif  TmBioBox(dna(i,:))>59&& TmBioBox(dna(i,:))<=69
    goal=goal+((69-TmBioBox(dna(i,:)))/10)*40;
else
    goal=goal+0;
end
[sim,~]=similarity(dna(i,:),dna(Virus_list,:),Virus_list,i);
[hme,~]=h_measure(dna(i,:),dna(Virus_list,:),Virus_list,i);

if sim/(size(Virus_list,2)-1)<=7
    goal=goal+30;
elseif sim/(size(Virus_list,2)-1)>7 && sim/(size(Virus_list,2)-1)<12
    goal=goal+(12-sim/(size(Virus_list,2)-1))*5;
elseif sim/(size(Virus_list,2)-1)>=12
    goal=goal+0;
end
if hme/(size(Virus_list,2)-1)<=7
    goal=goal+30;
elseif hme/(size(Virus_list,2)-1)>7 && hme/(size(Virus_list,2)-1)<12
    goal=goal+(12-hme/size(Virus_list,2))*5;
elseif hme/(size(Virus_list,2)-1)>=12 || sim/(size(Virus_list,2)-1)>=12
    goal=goal+0;
end
if gccontent(dna(i,:))~=0.5 || Hairpin(dna(i,:))~=0
    goal=0;
end
GTotal=goal;
end


function EvaluateNewSimHm(dna,Virus_list)

goal=[];
goal=repmat(goal,1,size(Virus_list,2));
global DNAInfo
for i=Virus_list
    goal(i)=0;
end
for i=Virus_list
    [hme,sim]=HmSm(dna(i,:),dna(Virus_list,:),Virus_list,i);
    con=continuity(dna(i,:));
    Tm=TmBioBox(dna(i,:));
    gc=gccontent(dna(i,:));
    hair=Hairpin(dna(i,:));
    if con==0
        goal(i)=goal(i)+20;
    elseif con==1
        goal(i)=goal(i)+10;
    elseif con==2
        goal(i)=goal(i)+5;
    else
        goal(i)=goal(i)+0;
    end
    if Tm<=59
        goal(i)=goal(i)+40;
    elseif Tm>59 && Tm<=69
        goal(i)=goal(i)+((69-Tm)/10)*40;
    else
        goal(i)=goal(i)+0;
    end
    %If you want to obtain low Tm variance sequences, replace program
    %Tm=zeros(1,length(Virus_list));
    %for ii=Virus_list
    %    [TM,~]=TmBioBox(dna(ii,:));
    %    Tm(ii)=TM;
    %end
    %VARTm=var(Tm);
    %if VARTm<=1
    %    goal(i)=goal(i)+40;
    %elseif VARTm>1 && VARTm<=3
    %    goal(i)=goal(i)+((3-VARTm)/2)*40;
    %else
    %    goal(i)=goal(i)+0;
    %end
    
    if sim/(size(Virus_list,2)-1)<=7
        goal(i)=goal(i)+30;
    elseif sim/(size(Virus_list,2)-1)>7 && sim/(size(Virus_list,2)-1)<13
        goal(i)=goal(i)+(13-sim/(size(Virus_list,2)-1))*5;
    elseif sim/(size(Virus_list,2)-1)>=13 
        goal(i)=goal(i)+0;
    end
    if hme/(size(Virus_list,2)-1)<=7
        goal(i)=goal(i)+30;
    elseif hme/(size(Virus_list,2)-1)>7 && hme/(size(Virus_list,2)-1)<13
        goal(i)=goal(i)+(13-hme/size(Virus_list,2))*5;
    elseif hme/(size(Virus_list,2)-1)>=13 || sim/(size(Virus_list,2)-1)>=13 
        goal(i)=goal(i)+0;
    end
    if gc==0.5 && hair==0
        goal(i)=goal(i);
    else
        goal(i)=0;
    end
    DNAInfo(i).Sequence=dna(i,:);
    DNAInfo(i).Similarity=sim;
    DNAInfo(i).H_measure=hme;
    DNAInfo(i).Continuity=con;
    DNAInfo(i).GCcontent=gc;
    DNAInfo(i).Hairpin=hair;
    DNAInfo(i).Tm=Tm;
    DNAInfo(i).Goal=goal(i);
    DNAInfo(i).Transmission=DNAInfo(i).Goal;
end
end


