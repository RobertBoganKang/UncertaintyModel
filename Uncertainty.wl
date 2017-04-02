(* ::Package:: *)

Clear["Global`*"];
(*GA program*)
Module[{a,anydim,calc,CC,Ccoe,child,cnt2orspce,coe,cvg,d,dgtz,dissb,dsetgs,dtheta,err1,err2,fercnr,flnpool10,frml,ii,ini,invmatrixA,iptmp,lamdas,m,mate,matrixA,matrixC,mCC,mid,mincase,Minins,mins,Minxpars,mutation,mxnum,mxpars,nsb,parents,password,pnts,pool,pool10,pos,posvar,prob,resl,run,s,s3d,sb,SCC,schemata,setgs,snd,stg,stl,str,sub20,subs,subs1,subsvar,tfs,theta,theta1,to10,to2,totheta,trial,tsts,utc,utcii,utcII,var,var2o,vars,varsee,wlstrg,xsb,\[Theta]1},
(**************************************************************************************************************************)
(*1. Demo Data*)
(*test Data*)
\[Theta]={{0.03,8.5,17.5},{0.03,5.5,2.5},{0.01,4.5,3.5},{0.01,1.5,0.5},{0.021,5.2,6.4}};
(*n should be an even number; the population number*)n=100;
(*the number of digits in binary list*)num=25;
(*Important Parameters*)
(*the highest order*)rr=2;
(*the number of variables*)mm=2;
(*Untwisted Coefficient*)utc=.3;
(*Untwisted Coefficient II*)utcII=0;
(*Converge Coefficient*)cvg=.8;
(*mutation rate parameter*)a=0.3;
(*schemata for minimum case*)schemata=False;
(*Play sound*)snd=True;
(**************************************************************************************************************************)


(**************************************************************************************************************************)
(*2. The initialization of program*)
ini:=((*n should be even number; others should be integer*)
n=If[OddQ@Round[Abs@n],Round[Abs@n]+1,Round[Abs@n]];{num,rr,mm}=Round@Abs@{num,rr,mm};
(*the digits*)dgtz=Length@IntegerDigits[2^(num-1)]-1;
(*the mid point*)mid=2^(num-1);
(*the maximum number; eliminate machine precision error*)mxnum=(2.^(num-1)-1)/10^dgtz;
(*Normalize C to Ccoe as max*)Ccoe=1.5;
(*convert to nominal form, from 0 to max number*)theta=2(Times[(Max/@Transpose@\[Theta]-Min/@Transpose@\[Theta])^-1,#]&/@((#-Min/@Transpose@\[Theta])&/@\[Theta]))-1;
(*convert back to original space*)totheta[x_]:=(((#+Min/@Transpose@\[Theta])&/@((#*(Max/@Transpose@\[Theta]-Min/@Transpose@\[Theta]))&/@((x+1)/2))));(*how many dimensions in theta*)dtheta=Dimensions[\[Theta]];
(*(*0\[LessEqual]\[Lambda]\[LessEqual]1*)tsts=And@@((0\[LessEqual]var[[#]]\[LessEqual]1)&/@Range[mm]);*)
(*(*\[Lambda]\[GreaterEqual]0*)tsts=And@@((0\[LessEqual]var[[#]])&/@Range[mm]);*)
tsts={};
(*create variables*)var=Table[ToExpression["\[Lambda]"<>ToString[i]],{i,mm}];
varsee=Table[ToExpression["\[Lambda]"<>ToString[i]]->Subscript[\[Lambda],i],{i,mm}];
(*substitute higher dimensions to 0*)sub20=If[Length@var>3,(#->0)&/@var[[4;;]],{}];
(*the coefficient of order*)coe=Select[Tuples[Range[rr+1]-1,mm],Total[#]<rr+1&];
(*Length of coefficient*)m=Length@coe;
(*variable matrix*)vars=Table[Times@@Table[var[[i]]^coe[[j,i]],{i,mm}],{j,m}];
(*the position of single variables*)posvar=Flatten[Position[Coefficient[vars,#],1]&/@var];
(*convert to decimal; on pool;*)
to10[x_]:=Table[((FromDigits[#,2]-mid)/10.^dgtz)&/@Partition[x[[i,j]],num],{i,Dimensions[x][[1]]},{j,Dimensions[x][[2]]}]/(mxnum*Ccoe);
(*convert to binary; on pool*)
to2[x_]:=Table[Flatten@IntegerDigits[Round[x[[i,j]]*mxnum*Ccoe*10.^dgtz]+mid,2,num],{i,Dimensions[x][[1]]},{j,Dimensions[x][[2]]}];
(*Mutation Function*)
mutation[x_]:=If[RandomReal[]<a*0.1*Minins,RandomInteger[],x];
(*the pool of selection*)
(*initialize the pool*)pool=Table[RandomInteger[],{n},{m},{num*dtheta[[2]]}];
ii=0;tfs=0;
mincase={10,0};);
ini;
(**************************************************************************************************************************)


(**************************************************************************************************************************)
(*3. GA Iterations*)run:=While[If[tfs==1&&If[trial<5,True,ii<100],True,False],(*convert pool into decimal*)pool10=to10[pool];
flnpool10=Abs[Flatten/@pool10];
(*F(\[Lambda]):=f(\[Lambda])-\[Theta]*)
err1=(pool10[[#,2;;]])&/@Range[n];err2=Table[(pool10[[i,1]]-#)&/@theta,{i,n}];
(*untwisted coefficient*)utcii=Flatten@Table[vars/.(#->5^utcII&/@var),{dtheta[[2]]}];
(*get series of minimum distance*)
mins=Table[5/dtheta[[1]]*Sqrt[Total[(Quiet@FindMinimum[Flatten@{#,tsts},Transpose@{var,Table[0,{mm}]},AccuracyGoal->dgtz-1]&/@Table[First@Expand@Total@((Transpose[Prepend[err1[[i]],err2[[i,j]]]].Transpose[{vars}])^2),{j,dtheta[[1]]}])[[;;,1]]]+10^-dgtz],{i,n}];
Minins=Min[mins];
(*get the position of min distance*)pos=Position[mins,Minins][[1,1]];
(*The maximum parameters in matrix C; less max parameter would be more convex*)
mxpars=(#/Max@#)&@(Total/@(utcii*#&/@flnpool10));
Minxpars=Min@mxpars;
(*probability parameter to choose parents*)prob=3^-(10*mins*Times@@((#*Max[#]^-1)&/@{mxpars-Minxpars+1000^(.8-utc),mins-10^-(1+10^(1-dgtz)-cvg)*Minins}));
(*select parents*)parents=Table[RandomChoice[prob->pool],{n/2},{2}];mate=(Flatten/@parents[[#]]&/@Range[n/2]);(*the random selection crossover technique and mutation at same time*)
child=Flatten[Table[Table[mutation[mate[[j,RandomInteger[]+1,i]]],{2},{i,num*dtheta[[2]]*m}],{j,n/2}],1];
(*update pool*)pool=Partition[#,num*dtheta[[2]]]&/@child;
(*store the min case ever happened in this program*)mincase=If[Minins<mincase[[1]],s[Minins+1,1];s[1,.05];{Minins,pool10[[pos]]},s[Minins+1,.01];mincase];
ii=ii+1;(*bomb for trial version*)If[trial<5,Null,Pause[2]];(*add schemata for minimum case*)If[schemata,pool10=Append[pool10[[2;;]],mincase[[2]]]];];
(**************************************************************************************************************************)


(**************************************************************************************************************************)
(*4. GA User's Interface*)
(*4.0 Regular Used UI tools*)
stl[x_]:=Style[x,Darker@Blue,Bold,FontFamily->"Times New Roman"];
stg[x_]:=Style[x,Darker[Green,.5],Bold,FontFamily->"Times New Roman"];
str[x_]:=Style[x,Darker@Red,Bold,FontFamily->"Times New Roman"];
(*beep function with 3D sound*)
s[x_,y_]:=If[snd,(s3d=If[x==0,10^-3,1.3 RandomReal[{-1,1}]];EmitSound[Play[{y*2^(-1*^3 t^2) Sin[440 x 2.\[Pi] t],y*2^(-1*^3 t^2) Sin[440 x 2.\[Pi] t+s3d]},{t,-.01,.08},PlayRange->{-1,1}]]),Null];
wlstrg:=("Welcome to Robert Bogan Kang's uncertainty model using Genetic Algorithm!"<>If[trial<5,"","Trial Version!"]);
(*4.1 UI body*)
d:=((*welcome*)If[RandomReal[]<0.001,Speak@wlstrg];(*Welcome Songs*)s[#,.5Sqrt[RandomReal[]]]&/@If[RandomReal[]<0.01,{3/2,0,1,5/3,0,3/2,5/4,0,1,2,0,5/3,3/2,2,3/2,7/5,0,7/5,3/2,0,0,4/3,0,0,6/5,3/2,16/7,2,0,16/9,24/25,1,16/9,8/5,0,3/2,16/15,9/8,6/5,7/5,3/2,8/5,5/4,4/3,3/2,5/3,16/9,8/5,3/2,0,0,6/5,0,3/2,5/4,0,9/8,2,0,8/5,3/2,0,3/2,7/5,8/5,7/5,3/2,0,0,3},{1,0,5/4,9/8,1,0,3/4,0,1,2,9/8,0,1,0,3/4,0,1,0,3/2,5/4,4/3,7/5,3/2,0,3/4,8/10,5/6,0,15/16,0,1}];CreateDialog[Dynamic[Column[{Row@{Framed[Column[{If[ii==0,Button["\:256e(\:256f_\:2570)\:256d",Do[s[.75 n,.5*3.5^-n];,{n,10}];setgs;,Appearance->"None"],Sequence@@{Dynamic[Mean[mins]],Dynamic[Style[Minins,Red]],Dynamic[Mean[prob]],Dynamic[ii]}]}],Background->If[trial<5,LightYellow,LightRed]],Button[If[ii==0,"",Tooltip[Style[".",White],"Reset",TooltipDelay->.5]],s[1,#]&/@{.4,.09,.04};ini;,Appearance->"None"]},Row[{Button["Start",s[1,.8];s[1.5,.1];tfs=1;run;,Method->"Queued"],Button[Tooltip[Style["|",Gray],"Dynamic Settings",TooltipDelay->0.5],Do[s[1.5 1/(n),.3*3^-(n-2)];,{n,2,7}];dsetgs;,Appearance->"None"],If[ii==0,"",Button[Tooltip[Style["|",Darker[Green]],"Results",TooltipDelay->0.5],Do[s[.75 (4-n),3^-(4-n)];,{n,4}];resl;,Appearance->"None"]],Button["Stop",tfs=0;s[1,.8];s[.75,.1];,Method->"Preemptive"]}]}]],Background->LightGray,WindowTitle->"GA"];);
(*settings*)
setgs:=CreateDialog[Column[{Button[Text["Test Data"],iptmp=Input["Please insert here for more data:"];\[Theta]=If[ToString@iptmp=="$Canceled",\[Theta],iptmp];,Method->"Queued",Appearance->None],InputField[Dynamic[\[Theta]],FieldSize->{Automatic,3}],Text["The population number of GA:"],InputField[Dynamic[n],Number],Text["The number of digits in binary list"],InputField[Dynamic[num],Number],Text["The highest order of Polynomial"],InputField[Dynamic[rr]],Text["The number of variables"],InputField[Dynamic[mm]],Button["OK",ini;DialogReturn[];,ImageMargins->5]}],WindowTitle->"Settings"];
(*DynamicSettings*)
dsetgs:=CreateDialog[OpenerView@{Row@{Tooltip[#,"Schemata for Minimum Case",TooltipDelay->.2]&@Checkbox[Dynamic[schemata]],Tooltip[#,"Play Sound?",TooltipDelay->.2]&@Checkbox[Dynamic[snd]]},Row@{Tooltip[#,"Untwisted Coefficient",TooltipDelay->.2]&@VerticalSlider[Dynamic[utcII],ImageSize->{16,70}],Tooltip[#,"Untwisted Coefficient to obvious",TooltipDelay->.2]&@VerticalSlider[Dynamic[utc],ImageSize->{16,70}],Tooltip[#,"Converging Coefficient",TooltipDelay->.2]&@VerticalSlider[Dynamic[cvg],ImageSize->{16,70}],Tooltip[#,"Mutation Coefficient",TooltipDelay->.2]&@VerticalSlider[Dynamic[a],ImageSize->{16,70}]}},WindowTitle->"DS"];
(*4.2 Generate Report Result*)
resl:=(cnt2orspce:=(fercnr=(Union@Flatten[Table[If[If[mm<3,True,(\[Lambda]1==1&&\[Lambda]2==1)||(\[Lambda]2==1&&\[Lambda]3==1)||(\[Lambda]3==1&&\[Lambda]1==1)||(\[Lambda]1==0&&\[Lambda]2==0)||(\[Lambda]2==0&&\[Lambda]3==0)||(\[Lambda]3==0&&\[Lambda]1==0)||(\[Lambda]1==0&&\[Lambda]2==1)||(\[Lambda]2==0&&\[Lambda]3==1)||(\[Lambda]3==0&&\[Lambda]1==1)||(\[Lambda]1==1&&\[Lambda]2==0)||(\[Lambda]2==1&&\[Lambda]3==0)||(\[Lambda]3==1&&\[Lambda]1==0)],(Evaluate@First@totheta@Transpose@SCC)/.sub20,{}],{\[Lambda]1,0,1,1/30},{\[Lambda]2,0,1,1/30},{\[Lambda]3,0,1,1/30}],2])[[2;;]];Which[dtheta[[2]]>=3,Show@{ListPointPlot3D[fercnr[[;;,anydim]],BoxRatios->1,ImageSize->Large,PlotStyle->Directive[PointSize[Small],Black]],ListPointPlot3D[\[Theta][[;;,anydim]],BoxRatios->1,PlotStyle->Directive[PointSize[Large],Red]]},dtheta[[2]]==2,Show[{ListPlot[fercnr,ImageSize->Large,PlotStyle->Directive[PointSize[Small],Black]],ListPlot[\[Theta],PlotStyle->Directive[PointSize[Large],Red]]}],True,{}]);
lamdas:=(Which[mm>2,Graphics3D[{Red,Sphere[pnts[[;;,;;3]],0.05]},ImageSize->Medium,Axes->True],mm==2,Graphics[{PointSize[0.05],Red,Point[pnts]},ImageSize->Medium,Frame->True],mm<2,NumberLinePlot[Flatten[pnts],ImageSize->Medium,PlotStyle->Directive[Red,PointSize->0.05]]]);
(*calculate report parameters*)
calc=((*Converting space from unbounded \[Lambda] to {0,1} domain*)
subs1=((Quiet[FindMinimum[Flatten@{#,tsts},Transpose@{var,Table[0,{mm}]},AccuracyGoal->dgtz-1]]&)/@Table[First[Expand[Total[(Transpose[mincase[[2]]].Transpose[{vars}]-theta[[j]])^2]]],{j,dtheta[[1]]}]);subs=subs1[[;;,2]];
sb=Transpose[subs[[;;,;;,2]]];xsb=Max/@sb;nsb=Min/@sb;
dissb=xsb-nsb;
subsvar=(var[[#]]->Expand[(var-nsb)/dissb][[#]])&/@Range[mm];frml=Expand[vars/.subsvar];
var2o=(#->0&/@var);(*the conversion matrix A*)invmatrixA=Table[If[j==1,frml[[i]],Coefficient[frml[[i]],vars[[j]]]],{i,m},{j,m}]/.var2o;
matrixA=Inverse@invmatrixA;
(*update subs \[Lambda]*)subs[[;;,;;,2]]=Flatten/@(((invmatrixA.Transpose[{vars}])/.#)&/@subs)[[;;,posvar]];
subs1[[;;,2]]=subs;
pnts=subs[[1;;All,1;;All,2]];
(*The matrix C, then convert to [C']=[C].[A]^-1*)
CC=Transpose[mincase[[2]]].matrixA;SCC=Expand[CC.Transpose[{vars}]];
(*Original space [C]*)
mCC=Flatten@Expand@totheta@{Flatten@SCC};
matrixC=Table[If[j==1,mCC[[i]],Coefficient[mCC[[i]],vars[[j]]]],{i,dtheta[[2]]},{j,m}]/.var2o;theta1=(Flatten[SCC/. #]&)/@subs;\[Theta]1=(Flatten[mCC/. #]&)/@subs;
(*choose any 3 dimensions in the surface fitting plot*)anydim=Sort@RandomSample[Range[dtheta[[2]]],3];);
(*Report UI*)
CreateDocument[Column[Panel/@Column/@Partition[#,2]&@Flatten@{str@"Export to D:",Button["Xpt2D",Quiet[CreateDirectory["D:\\#rbk#"]];Export["D:\\#rbk#\\result_"<>IntegerString[Round[10AbsoluteTime[]],36]<>".png",InputNotebook[]];Beep[];,Appearance->"Palette"],str@"{\[Lambda]}",Row[vars,",",Frame->True],stl@"[C']",MatrixForm[CC],stl@"f'(\[Lambda])=[C'].{\[Lambda]}",MatrixForm[SCC],stg@"[C]",MatrixForm[matrixC],stg@"f(\[Lambda])=[C].{\[Lambda]}",MatrixForm@Transpose@{mCC},str@"Error & [\!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\)]->[\[Lambda]]",TableForm[subs1],
stl@"||[\!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\)']-[\[Theta]']|| :>[\!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\)'] vs [\[Theta]']",Row[{MatrixForm[Norm/@(theta1-theta)],stl@" :> ",MatrixForm[theta1],stl@"v.s.",MatrixForm[theta]}],stg@"||[\!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\)]-[\[Theta]]|| :>[\!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\)] vs [\[Theta]]",Row[{MatrixForm[Norm/@(\[Theta]1-\[Theta])],stg@" :> ",MatrixForm[\[Theta]1],stg@"v.s.",MatrixForm[\[Theta]]}],str@"[\[Lambda]]",lamdas,stg@"\[Theta] & f(\[Lambda])",cnt2orspce}]/.varsee,WindowTitle->"Results"];);
(**************************************************************************************************************************)


(**************************************************************************************************************************)(*5. encrypt UI*)CreateDialog[trial=0;Dynamic@Column[{If[trial==5,"Trial Version!",Tooltip[#,"Password!",TooltipDelay->1]&@InputField[Dynamic[password],String,FieldHint->"Enter here!",FieldMasked->True,FieldSize->5.5]],Button[Dynamic@Which[trial<1,"Unlock",trial<5,"Unlock ("<>ToString@(5-trial)<>")",trial==5,"Accept!"],If[password==("rbk"<>StringReverse[StringJoin@@ToString/@DateList[][[2;;5]]])<>"!"||(MousePosition[][[1]]==0&&MousePosition[][[2]]<100&&CurrentValue["ShiftKey"]),DialogReturn[];d;,trial++;If[trial<6,s[.75,.8],d;DialogReturn[]]]]}],WindowTitle->"PW"];];
