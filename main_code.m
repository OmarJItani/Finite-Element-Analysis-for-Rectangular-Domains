clear all
clc
NH=-1;
NV=-1;
format long
Tv=[0 -1 0 0 0 0;1 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 -1 0;0 0 0 1 0 0;0 0 0 0 0 1];
Th=[1 0 0 0 0 0;0 1 0 0 0 0 ;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
DOFS_Node=3;

%%%%%%%%%%%%%%  GENERAL INPUTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = 'Specify the vertical dimension of the domain: ';
V= input(prompt);
%V=2;
prompt = 'Specify the horizontal dimension of the domain: ';
H= input(prompt);
%H=2;
prompt = 'To how many elements is the vertical side divided: ';
NV= input(prompt);
prompt = 'To how many elements is the horizontal side divided: ';
NH= input(prompt);
%NH=2;
%NV=2;
prompt = 'Now specify the cross sectional dimensions: b= ';
b= input(prompt);
%b=0.1;
prompt = '                                            h= ';
h= input(prompt);
%h=0.3;
prompt = 'Specify Youngs Modulus of Elasticity E: ';
 E(1:2*NH*NV+NH+NV) = input(prompt);
%E(1:2*NH*NV+NH+NV)=200e9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(NH==-1)
    if(V<=H)
        lh=V/5;
    else lh=H/5;
    end
    NH=H/lh;
else
  lh=H/NH;  
end
if(NV==-1)
    if(V<=H)
        lv=V/5;
    else lv=H/5;
    end
    NV=V/lv;
else
  lv=V/NV; 
end

A=b*h;
Iz=h^3*b/12;

elem=trail_elements(NH,NV) ;
Node=GenerateNodes_omar(H,NH,V,NV);
Node(:,2)=Node(:,2)*100;
Node(:,3)=Node(:,3)*100;
Node=round(Node);
F=zeros(3*size(Node,1),1);

%%%%%%%%%%%  CONCENTRATED LOAD INPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = 'How many nodes is there forces on: ';
numberofcloads= input(prompt);
%numberofcloads=1;
if (numberofcloads>=1)
prompt = 'Specify the nodes coordinates in a matrix form: ';
[lcoordinates]= input(prompt);
lcoordinates=lcoordinates*100;
%lcoordinates=[ 2*100 2*100];
prompt = 'Specify the nodes loads in a matrix form: ';
[clvalue]= input(prompt);
%clvalue=[0 1000 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lnodes=zeros(numberofcloads,2);
for cl=1:numberofcloads
    for i=1:(NH+1)*(NV+1)
        if (Node(i,2)==lcoordinates(cl,1) && Node(i,3)==lcoordinates(cl,2))
           lnodes(cl,1)=cl;
           lnodes(cl,2)=i; 
        end
    end
end
for cl=1:numberofcloads
F(DOFS_Node*lnodes(cl,2)-2,1)=F(DOFS_Node*lnodes(cl,2)-2,1)+clvalue(cl,1);
F(DOFS_Node*lnodes(cl,2)-1,1)=F(DOFS_Node*lnodes(cl,2)-1,1)+clvalue(cl,2);
F(DOFS_Node*lnodes(cl,2),1)=F(DOFS_Node*lnodes(cl,2),1)+clvalue(cl,3);
end
%%%%%%%%%   DISTRIBUTED LOAD INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = 'How many elements is there distributed forces on? ';
numberofdloads= input(prompt);
%numberofdloads=0;
if(numberofdloads>=1)
prompt = 'Specify the node coordinates of the elements that the distributed loads are applied on in a matrix form: ';
[dlcoordinates]= input(prompt);
dlcoordinates=dlcoordinates*100;
%dlcoordinates=[0 0 0 0 ;0 0 0 0 ];
prompt= 'Specify the value of the distributed loads: ';
q= input(prompt);
end
%q=[0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for dl=1:numberofdloads
    for j=1:(NH+1)*(NV+1)
        if (Node(j,2)==dlcoordinates(dl,1) && Node(j,3)==dlcoordinates(dl,2))
           dlnodes(dl,1)=dl;
           dlnodes(dl,2)=j; 
        end
        if (Node(j,2)==dlcoordinates(dl,3) && Node(j,3)==dlcoordinates(dl,4)) 
           dlnodes(dl,3)=j;
        end 
    end
end
for dl=1:numberofdloads
    if( dlnodes(dl,2)+1 == dlnodes(dl,3))
    force=(q(dl)*lh)/2;
    moment=(q(dl)*lh^2)/12;
    F(DOFS_Node*dlnodes(dl,2)-2,1)=F(DOFS_Node*dlnodes(dl,2)-2,1)+0;
    F(DOFS_Node*dlnodes(dl,2)-1,1)=F(DOFS_Node*dlnodes(dl,2)-1,1)+force;
    F(DOFS_Node*dlnodes(dl,2),1)=F(DOFS_Node*dlnodes(dl,2),1)-moment;
    F(DOFS_Node*dlnodes(dl,3)-2,1)=F(DOFS_Node*dlnodes(dl,3)-2,1)+0;
    F(DOFS_Node*dlnodes(dl,3)-1,1)=F(DOFS_Node*dlnodes(dl,3)-1,1)+force;
    F(DOFS_Node*dlnodes(dl,3),1)= F(DOFS_Node*dlnodes(dl,3),1)+moment;
    else
    force=(q(dl)*lv)/2;
    moment=(q(dl)*lv^2)/12;
    F(DOFS_Node*dlnodes(dl,2)-2,1)= F(DOFS_Node*dlnodes(dl,2)-2,1)+force;
    F(DOFS_Node*dlnodes(dl,2)-1,1)= F(DOFS_Node*dlnodes(dl,2)-1,1)+0;
    F(DOFS_Node*dlnodes(dl,2),1)=  F(DOFS_Node*dlnodes(dl,2),1)-moment;
    F(DOFS_Node*dlnodes(dl,3)-2,1)= F(DOFS_Node*dlnodes(dl,3)-2,1)+force;
    F(DOFS_Node*dlnodes(dl,3)-1,1)= F(DOFS_Node*dlnodes(dl,3)-1,1)+0;
    F(DOFS_Node*dlnodes(dl,3),1)= F(DOFS_Node*dlnodes(dl,3),1)+moment;  
    end
end
%%%%%%%%%%%%%%    VOIDS INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 prompt = 'How many voids? ';
 voidnum = input(prompt);

for i=1:voidnum
prompt = 'Specify the coordinates of the bottom-left node of the void: ';
[void] = input(prompt);
 
prompt = 'What is the horizantal length of the void? ';
L = input(prompt);
 
prompt = 'What is the vertical length of the void? ';
B = input(prompt);

void=void*100;

x=ismember( void(1,1),Node(:,2) );
y=ismember( void(1,2),Node(:,3) );
if x==1 && y==1
    [i]=find(Node(:,2)==void(1,1) );
    [j]=find(Node(:,3)==void(1,2) );
    C=intersect(i,j);    
end
%%%%%%%%%%%%forward

L=L/lh;

HH=zeros(round(L),1);
F1=C;
for z=1:L
%to move forward
if z==1
Cc=F1+z; 
end

x=ismember(F1,elem(:,2));% see where starting node is 
y=ismember(Cc,elem(:,3));% see where next node is 
if x==1 && y==1
    [i]=find(elem(:,2)==F1);% which elemnt both node are on
    [j]=find(elem(:,3)==Cc,1); %which element both node are on 
    HH(z,1)=intersect(i,j);     
end 

E(1,HH(z,1))=0;
F1=Cc;
Cc=Cc+1; 
end 
%%%upward Voids

B=B/lv;
VV=zeros(round(B),1);

G=C;
for z=1:B  
   Cc1=[ void(1,1) void(1,2)+((V/NV)*100)];
   
x=ismember( Cc1(1,1) , Node(:,2));
y=ismember( Cc1(1,2),Node(:,3));
if x==1 && y==1
    [i]=find(Node(:,2)==Cc1(1,1));
    [j]=find(Node(:,3)==Cc1(1,2));
    C1=intersect(i,j); 
end 
    
x=ismember(G,elem(:,2));% see where starting node is 
y=ismember(C1,elem(:,3));% see where next node is 
if x==1 && y==1
    [i]=find(elem(:,2)==G);% which elemnt both node are on
    [j]=find(elem(:,3)==C1); %which element both node are on 
    VV(z,1)=intersect(i,j);     
end 

E(1,VV(z,1))=0;
 G=C1;
 void(1,2)=void(1,2)+(V/NV)*100; 
end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=zeros(3*size(Node,1));
Total_Dof=1:1:DOFS_Node*size(Node,1);

 prompt = 'Is there any wall-supported nodes? (If yes enter 1, if no enter 0)? ';
 y1 = input(prompt);
 if(y1==1)
count=1;
prompt = 'Specify the wall-supported nodes: ';
 [restrained]= input(prompt);
%restrained=[ 1  ]; 
res_dof=zeros(3*size(restrained,2),1);

for n=1:3:3*size(restrained,2)
    
res_dof(n)=3*restrained(count)-2;
res_dof(n+1)=3*restrained(count)-1;
res_dof(n+2)=3*restrained(count);
count=count+1;
end
res_dof=res_dof';
 end

%%%%% ROLLER (NO V-DISP) %%%%%%%%%%
prompt = 'Is there any NO V-DISP roller-supported nodes? (If yes enter 1, if no enter 0)? ';
 y2 = input(prompt);
if (y2==1)
 countrv=1;
 prompt = 'Specify the NO V-DISP roller-supported nodes: ';
 [restrainedrv]= input(prompt);
 %restrainedrv=[0];
 res_dofrv=zeros(size(restrainedrv,2),1);
 
 for n=1:size(restrainedrv,2)
     
     res_dofrv(n)=3*restrainedrv(countrv)-1;
     countrv=countrv+1;
 end
 res_dofrv=res_dofrv';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%% ROLLER (NO H-DISP) %%%%%%%%%%
prompt = 'Is there any NO H-DISP roller-supported nodes? (If yes enter 1, if no enter 0)? ';
 y3 = input(prompt);
if (y3==1) 
countrh=1;
prompt = 'Specify the NO H-DISP roller-supported nodes: ';
 [restrainedrh]= input(prompt);
% restrainedrh=[4];
 res_dofrh=zeros(size(restrainedrh,2),1);
 
 for n=1:size(restrainedrh,2)
     
     res_dofrh(n)=3*restrainedrh(countrh)-2;
     countrh=countrh+1;
 end
 res_dofrh=res_dofrh';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%   PIN    %%%%%%%%%%%%%%
prompt = 'Is there any pin-supported nodes? (If yes enter 1, if no enter 0)? ';
 y4 = input(prompt);
 if (y4==1)
countp=1;
prompt = 'Specify the pin-supported nodes: ';
 [restrainedp]= input(prompt);
%restrainedp=[];
res_dofp=zeros(2*size(restrainedp,2),1);

for n=1:2:2*size(restrainedp,2)
    res_dofp(n)=3*restrainedp(countp)-2;
    res_dofp(n+1)=3*restrainedp(countp)-1;
    countp=countp+1;
end
res_dofp=res_dofp';
 end

for m=1:NH*(NV+1)
    
ke(1:6,1:6,m)=E(m).*[A/lh 0 0  -A/lh 0 0 ; ...
       0 12*Iz/lh^3   6*Iz/lh^2 0 -12*Iz/lh^3  6*Iz/lh^2; ...
       0 6*Iz/lh^2  4*Iz/lh  0 -6*Iz/lh^2  2*Iz/lh ; ...
       -A/lh 0 0 A/lh 0 0 ; ... 
       0 -12*Iz/lh^3 -6*Iz/lh^2 0 12*Iz/lh^3  -6*Iz/lh^2 ; ...
       0 6*Iz/lh^2  2*Iz/lh 0 -6*Iz/lh^2  4*Iz/lh ];
   
   if NV>=1
      ke(:,:,m)=Th*ke(1:6,1:6,m)*Th';
   end
end

for m=NH*(NV+1)+1:2*NH*NV+NH+NV
    
ke(1:6,1:6,m)=E(m).*[A/lv 0 0  -A/lv 0 0 ; ...
       0 12*Iz/lv^3   6*Iz/lv^2 0 -12*Iz/lv^3  6*Iz/lv^2; ...
       0 6*Iz/lv^2  4*Iz/lv  0 -6*Iz/lv^2  2*Iz/lv ; ...
       -A/lv 0 0 A/lv 0 0 ; ... 
       0 -12*Iz/lv^3 -6*Iz/lv^2 0 12*Iz/lv^3  -6*Iz/lv^2 ; ...
       0 6*Iz/lv^2  2*Iz/lv 0 -6*Iz/lv^2  4*Iz/lv ]; 
   
         ke(:,:,m)=Tv*ke(1:6,1:6,m)*Tv';  
end

for m=1:2*NH*NV+NH+NV
    
   dof(m,:)=[ DOFS_Node*elem(m,2)-2 DOFS_Node*elem(m,2)-1 DOFS_Node*elem(m,2) ...
             DOFS_Node*elem(m,3)-2 DOFS_Node*elem(m,3)-1 DOFS_Node*elem(m,3)]; %DOF's of each element

for i1 =1:2*DOFS_Node 
        for j1= 1:2*DOFS_Node     
            K(dof(m,i1),dof(m,j1))=K(dof(m,i1),dof(m,j1))+ke(i1,j1,m);      
        end
end
end
if (y1==1)
Free=setxor(Total_Dof,res_dof);
else
    Free=Total_Dof;
end
if (y2==1)
Free=setxor(Free,res_dofrv);
end
if (y3==1)
Free=setxor(Free,res_dofrh);
end
if (y4==1)
Free=setxor(Free,res_dofp);
end

kk=K(Free,Free);
ff=F(Free);
disp(Free)=inv(kk)*ff;

