*------------------------------------------------------;
*Natalie LePera
*STA-513-01-90
*Final Project
*Due: 09 Dec 2025 (8:45pm)
*------------------------------------------------------;
*------------------------------------------------------;



libname exam '/home/u63760461/STA513/Exams';
proc format;
value tx 1 = "IDC" 2 = "CT" 3 = "SE" 4 = "GDC";
run;

*------------------------------------------------------;
*data prep;
*------------------------------------------------------;

data coke;
set exam.finalcoke;
run;

data coke36;
set coke;
where month in (3,6);
run;

data coke3;
set coke36;
where month = 3;
run;

data coke6;
set coke36;
where month = 6;
run;

*------------------------------------------------------;
*EDA
*------------------------------------------------------;

proc contents data= coke36;
run;

proc sort data = coke36;
by month;
run;

proc freq data = coke36;
by month;
table domi;
run;

proc univariate data = coke36 normal;
var domi;
run;

proc univariate data = coke36 normal;
by month;
var domi;
run;


proc sgplot data = coke36;
title "Distribution of Raw DOMI Scores";
histogram domi / group = month fillattrs= (transparency = 0.5);
run;

proc univariate data = coke36;
by month;
var domi;
histogram / normal;
run;



*------------------------------------------------------;
*model 1
*------------------------------------------------------;

proc mixed data = coke36;
*by month;
class month sites patno;
model domi = tx_cond GSI SUIPR ARS BEH_ARS 
M0ASIDR PSGSI PS_DOMI PS_COLD PS_SOAV PS_HAM27 / solution cl;
random sites / solution cl;
repeated month / sub = patno;
ods output SolutionF = est;
run;


proc mixed data = coke36;
*by month;
class month sites patno;
model domi = tx_cond GSI SUIPR ARS BEH_ARS 
M0ASIDR PSGSI PS_DOMI PS_COLD PS_SOAV PS_HAM27 / solution cl;
repeated month / sub = patno;
ods output SolutionF = est;
run;


proc mixed data = coke36;
*by month;
class month sites ;
model domi = tx_cond GSI SUIPR ARS BEH_ARS 
M0ASIDR PSGSI PS_DOMI PS_COLD PS_SOAV PS_HAM27 / solution cl;
*repeated month / sub = patno;
ods output SolutionF = est;
run;

proc glm data = coke36;
model domi = tx_cond GSI SUIPR ARS BEH_ARS 
M0ASIDR PSGSI PS_DOMI PS_COLD PS_SOAV PS_HAM27 /ss3 solution clm;
run;

proc glm data = coke36 PLOTS=DIAGNOSTICS(UNPACK) ;
class tx_cond;
model domi = tx_cond GSI SUIPR ARS BEH_ARS 
M0ASIDR PSGSI PS_DOMI PS_COLD PS_SOAV PS_HAM27 /ss3 solution clparm;
output out = residm1 p = p r = resid student = studentized rstudent = jknife;
run;

proc univariate data = residm1 normal;
var resid;
run;

*reisudals not normal, time to transreg for box cox;

data coke36d;
set coke36;
tx1 = (tx_cond = 1);
tx2 = (tx_cond = 2);
tx3 = (tx_cond = 3);
run;

/*
proc transreg data = coke36d pbo;
model boxcox(domi) = identity(tx1 tx2 tx3);
run;
*/

proc freq data = coke36d;
tables domi;
run;

data coke36s;
set coke36d;
domi_shift = domi + 0.1;
run;

proc transreg data = coke36s pbo;
model boxcox(domi_shift) = identity(tx1 tx2 tx3);
run;

data coke36s;
set coke36s;
domibc = ((domi_shift**0.25)-1)/0.25;
run;

proc glm data = coke36s PLOTS=DIAGNOSTICS(UNPACK) ;
class tx_cond;
model domibc = tx_cond GSI SUIPR ARS BEH_ARS 
M0ASIDR PSGSI PS_DOMI PS_COLD PS_SOAV PS_HAM27 /ss3 solution clparm;
output out = residm1 p = p r = resid student = studentized rstudent = jknife;
run;

*MODEL 1;
proc glm data = coke36s PLOTS=DIAGNOSTICS(UNPACK) ;
class tx_cond;
model domibc = tx_cond GSI SUIPR PS_DOMI PS_COLD PS_SOAV  /ss3 solution clparm;
output out = residm1 p = p r = resid student = studentized rstudent = jknife;
ods output ParameterEstimates = est1;
run;

proc univariate data = residm1 normal;
var resid;
histogram / normal ;
run;


proc print data = est1;
format estimate 6.3;
run;



proc sort data = coke36s;
by tx_cond;
run;

proc means data = coke36s ;
var domibc;
by tx_cond;
output out = m1means mean = mdomibc;
run;

data coke36s;
merge coke36s m1means;
by tx_cond;
drop _type_ _freq_;
run;

data coke36s;
set coke36s;
where tx_cond ne .;
run;

proc sgplot data = coke36s;
xaxis type = discrete label = "Assigned Treatment";
yaxis label = "Log Transformed & Shifted DOMI Score";
label domibc = "DOMI (transformed)";
label mdomibc = "Mean DOMI (transformed)";
title "Transformed DOMI Score by Treatment Group";
vbox domibc / group = tx_cond;
scatter x = tx_cond y = mdomibc / markerattrs=(symbol = starfilled);
format tx_cond tx.;
run;


proc glm data = coke36s plots=diagnostics(unpack);
class tx_cond crack;
model domi = tx_cond|crack / solution clparm  ;
	lsmeans tx_cond*crack / cl pdiff tdiff;
output out = estm2 p = p r = resid;
run;

proc univariate data = estm2 normal;
var resid;
histogram / normal;
run;

data cokem2;
set coke36d;
crack1 = (crack = 1);
domi_shift = domi + 0.1;
drop crack;
run;

data cokem2;
set cokem2;
where domi ne .;
run;

proc transreg data = cokem2 pbo;
model boxcox(domi_shift) = identity (tx1 tx2 tx3 crack1 tx1*crack1 tx2*crack1 tx3*crack1);
run;

proc glm data = coke36s plots=diagnostics(unpack);
class tx_cond crack;
model domibc = tx_cond|crack / solution clparm  ;
	lsmeans tx_cond*crack / cl pdiff tdiff;
output out = estm3 p = p r = resid;
run;

proc univariate data = estm3 normal;
var resid;
histogram / normal;
run;

*MODEL 2:;
proc glm data = coke36s plots=diagnostics(unpack);
class tx_cond crack;
model domibc = tx_cond|crack / solution clparm  ;
	lsmeans tx_cond*crack / cl pdiff tdiff slice = crack;

estimate 'IDC crack' intercept 1 tx_cond 1 0 0 0 crack 1 0 tx_cond*crack 1 0 0 0 0 0 0 0;
estimate 'CT crack' intercept 1 tx_cond 0 1 0 0 crack 1 0 tx_cond*crack 0 0 1 0 0 0 0 0;
estimate 'SE crack' intercept 1 tx_cond 0 0 1 0 crack 1 0 tx_cond*crack 0 0 0 0 1 0 0 0;
estimate 'GDC crack' intercept 1 tx_cond 0 0 0 1 crack 1 0 tx_cond*crack 0 0 0 0 0 0 1 0;
run;


*----------------------------------;

proc glm data = coke36s plots=diagnostics(unpack);
class mar_stat job gths;
model domi = mar_stat | job | gths / solution clparm ;
output out = estm4 p = p r = resid;
run;

proc univariate data = estm4 normal;
var resid;
histogram / normal;

data cokem3;
set coke36d;
domi_shift = domi + 0.1;
mar1 = (mar_stat = 1);
job1 = (job = 1);
ghse = (gths = 1);
drop mar_stat job gths;
run;

/*
proc transreg data = cokem3 pbo;
model boxcox(domi) = identity(mar1 job1 ghse mar1*job1 mar1*ghse job1*ghse mar1*job1*ghse);
run;
*/

proc transreg data = cokem3 pbo;
model boxcox(domi_shift) = identity(mar1 job1 ghse mar1*job1 mar1*ghse job1*ghse mar1*job1*ghse);
run;

proc glm data = coke36s plots=diagnostics(unpack);
class mar_stat job gths;
model domibc = mar_stat | job | gths / solution clparm ;
output out = estm5 p = p r = resid;
run;
proc univariate data = estm5 normal;
var resid;
histogram / normal;

*model 3;
proc glm data = coke36s plots=diagnostics(unpack);
class mar_stat job gths;
model domibc = mar_stat | job | gths / solution clparm ;
	lsmeans mar_stat*job*gths / pdiff tdiff cl;
	
estimate 'mar emp col' intercept 1 
						mar_stat 1 0 job 1 0 gths 1 0 
						mar_stat*job 1 0 0 0 mar_stat*gths 1 0 0 0 job*gths 1 0 0 0
						mar_stat*job*gths 1 0 0 0 0 0 0 0;
estimate 'mar emp hs'  intercept 1 
						mar_stat 1 0 job 1 0 gths 0 1 
						mar_stat*job 1 0 0 0 mar_stat*gths 0 1 0 0 job*gths 0 1 0 0
						mar_stat*job*gths 0 1 0 0 0 0 0 0;
estimate 'mar une col' intercept 1 
						mar_stat 1 0 job 0 1 gths 1 0 
						mar_stat*job 0 1 0 0 mar_stat*gths 1 0 0 0 job*gths 0 0 1 0
						mar_stat*job*gths 0 0 1 0 0 0 0 0;
estimate 'mar une hs' intercept 1 
						mar_stat 1 0 job 0 1 gths 0 1 
						mar_stat*job 0 1 0 0 mar_stat*gths 0 1 0 0 job*gths 0 0 0 1
						mar_stat*job*gths 0 0 0 1 0 0 0 0;						

							
estimate 'umr emp col' intercept 1 
						mar_stat 0 1 job 1 0 gths 1 0 
						mar_stat*job 0 0 1 0 mar_stat*gths 0 0 1 0 job*gths 1 0 0 0
						mar_stat*job*gths 0 0 0 0 1 0 0 0;
estimate 'umr emp hs'  intercept 1 
						mar_stat 0 1 job 1 0 gths 0 1 
						mar_stat*job 0 0 1 0 mar_stat*gths 0 0 0 1 job*gths 0 1 0 0
						mar_stat*job*gths 0 0 0 0 0 1 0 0;
estimate 'umr une col' intercept 1 
						mar_stat 0 1 job 0 1 gths 1 0 
						mar_stat*job 0 0 0 1 mar_stat*gths 0 0 1 0 job*gths 0 0 1 0
						mar_stat*job*gths 0 0 0 0 0 0 1 0;
estimate 'umr une hs' intercept 1 
						mar_stat 0 1 job 0 1 gths 0 1 
						mar_stat*job 0 0 0 1 mar_stat*gths 0 0 0 1 job*gths 0 0 0 1
						mar_stat*job*gths 0 0 0 0 0 0 0 1;	
						
estimate 'avmarried' intercept 1 
						mar_stat 1 0 job 0.5 0.5 gths 0.5 0.5 
						mar_stat*job 0.5 0.5 0 0 mar_stat*gths 0.5 0.5 0 0 job*gths 0.25 0.25 0.25 0.25
						mar_stat*job*gths 0.25 0.25 0.25 0.25 0 0 0 0;	
estimate 'avunmarried' intercept 1 
						mar_stat 0 1 job 0.5 0.5 gths 0.5 0.5 
						mar_stat*job 0 0 0.5 0.5 mar_stat*gths 0 0 0.5 0.5 job*gths 0.25 0.25 0.25 0.25
						mar_stat*job*gths 0 0 0 0 0.25 0.25 0.25 0.25;	
estimate 'avemp' intercept 1 
						mar_stat 0.5 0.5 job 1 0 gths 0.5 0.5 
						mar_stat*job 0.5 0 0.5 0 mar_stat*gths 0.25 0.25 0.25 0.25 job*gths 0.5 0.5 0 0
						mar_stat*job*gths 0.25 0.25 0 0 0.25 0.25 0 0;	
estimate 'avunemp' intercept 1 
						mar_stat 0.5 0.5 job 0 1 gths 0.5 0.5 
						mar_stat*job 0 0.5 0 0.5 mar_stat*gths 0.25 0.25 0.25 0.25 job*gths 0 0 0.5 0.5
						mar_stat*job*gths 0 0 0.25 0.25 0 0 0.25 0.25;	
estimate 'avcol' intercept 1 
						mar_stat 0.5 0.5 job 0.5 0.5 gths 1 0 
						mar_stat*job 0.25 0.25 0.25 0.25 mar_stat*gths 0.5 0 0.5 0 job*gths 0.5 0 0.5 0
						mar_stat*job*gths 0.25 0 0.25 0 0.25 0 0.25 0;	
estimate 'avhs' intercept 1 
						mar_stat 0.5 0.5 job 0.5 0.5 gths 0 1 
						mar_stat*job 0.25 0.25 0.25 0.25 mar_stat*gths 0 0.5 0 0.5 job*gths 0 0.5 0 0.5
						mar_stat*job*gths 0 0.25 0 0.25 0 0.25 0 0.25;	
												

run;



*model 1;
proc sort data = coke36s;
by month;
run;

data cokem3;
set coke36s;
where month = 3;
run;

data cokem6;
set coke36s;
where month = 6;
run;


proc glm data = coke36s PLOTS=diagnostics(unpack) ;
by month;
class tx_cond;
model domibc = tx_cond GSI SUIPR PS_DOMI PS_COLD PS_SOAV  /ss3 solution clparm;
output out = residm1 p = p r = resid student = studentized rstudent = jknife;
ods output ParameterEstimates = est1;
run;

proc glm data = cokem3 PLOTS=diagnostics(unpack) ;
*by month;
class tx_cond;
model domibc = tx_cond GSI SUIPR PS_DOMI PS_COLD PS_SOAV  /ss3 solution clparm;
output out = residm13 p = p r = resid student = studentized rstudent = jknife;
ods output ParameterEstimates = est13;
run;

proc glm data = cokem6 PLOTS=diagnostics(unpack) ;
*by month;
class tx_cond;
model domibc = tx_cond GSI SUIPR PS_DOMI PS_COLD PS_SOAV  /ss3 solution clparm;
output out = residm16 p = p r = resid student = studentized rstudent = jknife;
ods output ParameterEstimates = est16;
run;


proc sgplot data =residm13;
xaxis type = discrete label = "Treatment";
yaxis label = "Transformed DOMI Score";
title 'Transformed DOMI Score by Treatment Type: Month 3';
scatter x = tx_cond y = domibc / group = tx_cond;
series x = tx_cond y = mdomibc;
scatter x = tx_cond y = mdomibc/ group = tx_cond markerattrs=(symbol = starfilled size = 11);
run;

proc sgplot data =residm16;
xaxis type = discrete label = "Treatment";
yaxis label = "Transformed DOMI Score";
title 'Transformed DOMI Score by Treatment Type: Month 6';
scatter x = tx_cond y = domibc / group = tx_cond;
series x = tx_cond y = mdomibc;
scatter x = tx_cond y = mdomibc/ group = tx_cond markerattrs=(symbol = starfilled size = 11);
run;

proc univariate data = residm13 normal;
var resid;
histogram / normal;
run;

proc univariate data = residm16 normal;
var resid;
histogram / normal;
run;


proc glm data = cokem3 plots=diagnostics(unpack);
class tx_cond crack;
model domibc = tx_cond|crack / solution clparm  ;
	lsmeans tx_cond*crack / cl pdiff tdiff slice = crack adjust = bon;
output out = residm23 p = p r = resid;
run;

proc univariate data = residm23 normal;
var resid;
histogram / normal;
run;

proc glm data = cokem6 plots=diagnostics(unpack);
class tx_cond crack;
model domibc = tx_cond|crack / solution clparm  ;
	lsmeans tx_cond*crack / cl pdiff tdiff slice = crack adjust = bon;
output out = residm26 p = p r = resid;
run;

proc univariate data = residm26 normal;
var resid;
histogram / normal;
run;

proc mixed data = cokem3;
class tx_cond crack;
model domibc = tx_cond|crack / solution cl  ;
	lsmeans tx_cond*crack / cl pdiff  slice = crack adjust = bon;

estimate 'IDC crack' intercept 1 tx_cond 1 0 0 0 crack 1 0 tx_cond*crack 1 0 0 0 0 0 0 0;
estimate 'IDC coke' intercept 1 tx_cond 1 0 0 0 crack 0 1 tx_cond*crack 0 1 0 0 0 0 0 0;
estimate 'CT crack' intercept 1 tx_cond 0 1 0 0 crack 1 0 tx_cond*crack 0 0 1 0 0 0 0 0;
estimate 'CT coke' intercept 1 tx_cond 0 1 0 0 crack 0 1 tx_cond*crack 0 0 0 1 0 0 0 0;
estimate 'SE crack' intercept 1 tx_cond 0 0 1 0 crack 1 0 tx_cond*crack 0 0 0 0 1 0 0 0;
estimate 'SE coke' intercept 1 tx_cond 0 0 1 0 crack 0 1 tx_cond*crack 0 0 0 0 0 1 0 0;
estimate 'GDC crack' intercept 1 tx_cond 0 0 0 1 crack 1 0 tx_cond*crack 0 0 0 0 0 0 1 0;
estimate 'GDC coke' intercept 1 tx_cond 0 0 0 1 crack 0 1 tx_cond*crack 0 0 0 0 0 0 0 1;

estimate 'PS Crack' intercept 1 tx_cond 0 0.5 0.5 0 crack 1 0 tx_cond*crack 0 0 0.5 0 0.5 0 0 0;
estimate 'DC Crack' intercept 1 tx_cond 0.5 0 0 0.5 crack 1 0 tx_cond*crack 0.5 0 0 0 0 0 0.5 0;
estimate 'PS-DC Crack' intercept 0 tx_cond -0.5 0.5 0.5 -0.5 crack 0 0 tx_cond*crack -0.5 0 0.5 0 0.5 0 -0.5 0 / cl;

estimate 'PS Coke' intercept 1 tx_cond 0 0.5 0.5 0 crack 0 1 tx_cond*crack 0 0 0 0.5 0 0.5 0 0;
estimate 'DC Coke' intercept 1 tx_cond 0.5 0 0 0.5 crack 0 1 tx_cond*crack 0 0.5 0 0 0 0 0 0.5;
estimate 'PS-DC Coke' intercept 0 tx_cond -0.5 0.5 0.5 -0.5 crack 0 0 tx_cond*crack 0 -0.5 0 0.5 0 0.5 0 -0.5/ cl;

estimate 'Crack - Coke' intercept 0 tx_cond 0 0 0 0 crack 0 0 tx_cond*crack -0.5 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 / cl;

run;




proc mixed data = cokem6;
class tx_cond crack;
model domibc = tx_cond|crack / solution cl  ;
	lsmeans tx_cond*crack / cl pdiff  slice = crack adjust = bon;

estimate 'IDC crack' intercept 1 tx_cond 1 0 0 0 crack 1 0 tx_cond*crack 1 0 0 0 0 0 0 0;
estimate 'IDC coke' intercept 1 tx_cond 1 0 0 0 crack 0 1 tx_cond*crack 0 1 0 0 0 0 0 0;
estimate 'CT crack' intercept 1 tx_cond 0 1 0 0 crack 1 0 tx_cond*crack 0 0 1 0 0 0 0 0;
estimate 'CT coke' intercept 1 tx_cond 0 1 0 0 crack 0 1 tx_cond*crack 0 0 0 1 0 0 0 0;
estimate 'SE crack' intercept 1 tx_cond 0 0 1 0 crack 1 0 tx_cond*crack 0 0 0 0 1 0 0 0;
estimate 'SE coke' intercept 1 tx_cond 0 0 1 0 crack 0 1 tx_cond*crack 0 0 0 0 0 1 0 0;
estimate 'GDC crack' intercept 1 tx_cond 0 0 0 1 crack 1 0 tx_cond*crack 0 0 0 0 0 0 1 0;
estimate 'GDC coke' intercept 1 tx_cond 0 0 0 1 crack 0 1 tx_cond*crack 0 0 0 0 0 0 0 1;

estimate 'PS Crack' intercept 1 tx_cond 0 0.5 0.5 0 crack 1 0 tx_cond*crack 0 0 0.5 0 0.5 0 0 0;
estimate 'DC Crack' intercept 1 tx_cond 0.5 0 0 0.5 crack 1 0 tx_cond*crack 0.5 0 0 0 0 0 0.5 0;
estimate 'PS-DC Crack' intercept 0 tx_cond -0.5 0.5 0.5 -0.5 crack 0 0 tx_cond*crack -0.5 0 0.5 0 0.5 0 -0.5 0 / cl;

estimate 'PS Coke' intercept 1 tx_cond 0 0.5 0.5 0 crack 0 1 tx_cond*crack 0 0 0 0.5 0 0.5 0 0;
estimate 'DC Coke' intercept 1 tx_cond 0.5 0 0 0.5 crack 0 1 tx_cond*crack 0 0.5 0 0 0 0 0 0.5;
estimate 'PS-DC Coke' intercept 0 tx_cond -0.5 0.5 0.5 -0.5 crack 0 0 tx_cond*crack 0 -0.5 0 0.5 0 0.5 0 -0.5/ cl;

estimate 'Crack - Coke' intercept 0 tx_cond 0 0 0 0 crack 0 0 tx_cond*crack -0.5 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 / cl;

run;




proc sgplot data =residm23;
xaxis type = discrete label = "Treatment";
yaxis label = "Transformed DOMI Score";
title 'Transformed DOMI Score by Treatment Type & Cocaine Type: Month 3';
scatter x = tx_cond y = domibc / group = crack;
series x = tx_cond y = p / group = crack;
scatter x = tx_cond y = p/ group = crack markerattrs=(symbol = starfilled size = 11);
run;

proc sgplot data =residm26;
xaxis type = discrete label = "Treatment";
yaxis label = "Transformed DOMI Score";
title 'Transformed DOMI Score by Treatment Type & Cocaine Type: Month 6';
scatter x = tx_cond y = domibc / group = crack;
series x = tx_cond y = p / group = crack;
scatter x = tx_cond y = p/ group = crack markerattrs=(symbol = starfilled size = 11);
run;





*model 3;
proc glm data = cokem3 plots=diagnostics(unpack);
class mar_stat job gths;
model domibc = mar_stat | job | gths / solution clparm ;
	lsmeans mar_stat*job*gths / pdiff tdiff cl;
output out = residm33 p = p r = resid;
run;

proc univariate data = residm33 normal;
var resid;
run;


proc glm data = cokem6 plots=diagnostics(unpack);
class mar_stat job gths;
model domibc = mar_stat | job | gths / solution clparm ;
	lsmeans mar_stat*job*gths / pdiff tdiff cl;
output out = residm36 p = p r = resid;
run;

proc univariate data = residm36 normal;
var resid;
run;



proc sgpanel data = residm33;
	panelby gths / columns = 2;
	scatter x = mar_stat y = domibc / group = job;
	series x = mar_stat y = p / group = job;
	scatter x = mar_stat y = p / group = job markerattrs=(symbol = starfilled size = 11);
run;

proc sgpanel data = residm36;
	panelby gths / columns = 2;
	scatter x = mar_stat y = domibc / group = job;
	series x = mar_stat y = p / group = job;
	scatter x = mar_stat y = p / group = job markerattrs=(symbol = starfilled size = 11);
run;


proc mixed data = cokem3 ;
class mar_stat job gths;
model domibc = mar_stat | job | gths / solution cl ;
	lsmeans mar_stat*job*gths / pdiff cl adjust = tukey;

estimate 'Secure - Insecure' intercept 0 
						mar_stat 1 -1 job 1 -1 gths 1 -1 
						mar_stat*job 1 0 0 -1 mar_stat*gths 1 0 0 -1 job*gths 1 0 0 -1
						mar_stat*job*gths 1 0 0 0 0 0 0 -1;				
	
estimate 'Mar - Unmar' intercept 0 
						mar_stat 1 -1 job 0 0 gths 0 0
						mar_stat*job 0.5 0.5 -0.5 -0.5 mar_stat*gths 0.5 0.5 -0.5 -0.5 job*gths 0 0 0 0
						mar_stat*job*gths 0.25 0.25 0.25 0.25 -0.25 -0.25 -0.25 -0.25;
			
estimate 'Emp-Unem' intercept 0 
						mar_stat 0 0 job 1 -1 gths 0 0
						mar_stat*job 0.5 -0.5 0.5 -0.5 mar_stat*gths 0 0 0 0 job*gths 0.5 0.5 -0.5 -0.5
						mar_stat*job*gths 0.25 0.25 -0.25 -0.25 0.25 0.25 -0.25 -0.25;	
	
estimate 'Col-HS' intercept 0 
						mar_stat 0 0 job 0 0 gths 1 -1 
						mar_stat*job 0 0 0 0 mar_stat*gths 0.5 -0.5 0.5 -0.5 job*gths 0.5 -0.5 0.5 -0.5
						mar_stat*job*gths 0.25 -0.25 0.25 -0.25 0.25 -0.25 0.25 -0.25 / cl;										
run;



proc mixed data = cokem6 ;
class mar_stat job gths;
model domibc = mar_stat | job | gths / solution cl ;
	lsmeans mar_stat*job*gths / pdiff cl adjust = tukey;

estimate 'Secure - Insecure' intercept 0 
						mar_stat 1 -1 job 1 -1 gths 1 -1 
						mar_stat*job 1 0 0 -1 mar_stat*gths 1 0 0 -1 job*gths 1 0 0 -1
						mar_stat*job*gths 1 0 0 0 0 0 0 -1;				
	
estimate 'Mar - Unmar' intercept 0 
						mar_stat 1 -1 job 0 0 gths 0 0
						mar_stat*job 0.5 0.5 -0.5 -0.5 mar_stat*gths 0.5 0.5 -0.5 -0.5 job*gths 0 0 0 0
						mar_stat*job*gths 0.25 0.25 0.25 0.25 -0.25 -0.25 -0.25 -0.25;
			
estimate 'Emp-Unem' intercept 0 
						mar_stat 0 0 job 1 -1 gths 0 0
						mar_stat*job 0.5 -0.5 0.5 -0.5 mar_stat*gths 0 0 0 0 job*gths 0.5 0.5 -0.5 -0.5
						mar_stat*job*gths 0.25 0.25 -0.25 -0.25 0.25 0.25 -0.25 -0.25;	
	
estimate 'Col-HS' intercept 0 
						mar_stat 0 0 job 0 0 gths 1 -1 
						mar_stat*job 0 0 0 0 mar_stat*gths 0.5 -0.5 0.5 -0.5 job*gths 0.5 -0.5 0.5 -0.5
						mar_stat*job*gths 0.25 -0.25 0.25 -0.25 0.25 -0.25 0.25 -0.25 / cl;										
run;


