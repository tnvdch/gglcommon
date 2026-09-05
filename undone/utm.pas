(****************************************************************************) 
Unit Mapproj;					    (*     21.June 1996     *)
                                                    (************************) 
                                                    (*  ETH ZUERICH GGL-IGP *) 
                                                    (*          (Ge)        *)
(****************************************************************************)
Interface
uses KOORD_TRANS,Ellgeod;
type zone_parameter = record
			 x0	 : extended;	{false origin}			 
			 y0	 : extended;	{false origin}
			 phi0	 : extended;	{center of projection}
			 lambda0 : extended;	{center of projection}
			 dphi	 : extended;	{zone width in phi}
			 dlambda : extended;	{zone width in lambda}
			 k0	 : extended;	{scale}
		      end;
type national_map   = record
			 Name	   : string;
			 proj_name : string;
			 proj_typ  : integer; 
			 ellipsoid : string; 
			 Datum_par : array[1..7] of extended;
			 Num_zones : integer;
			 Zone_Par  : array[1..5] of zone_parameter;
		      end;
const
   Grdeg = 180/200;
   GrRad = pi/200;
   rho	 = pi/180;
	 


Procedure Trans_Mercator_Geo_to_Grid(	 phi0,lambda0,x0,y0,a,e2,phi,lambda : extended;
				     var x,y				    : extended);
Procedure Trans_Mercator_Grid_to_Geo(	 phi0,lambda0,x0,y0,a,e2,x,y : extended;
				     var phi,lambda		     : extended);
Procedure Gauss_Krueger_Geo_to_Grid(	phi,lambda : extended;
				    var east,north : extended;
				    var zone       : integer); 
Procedure Gauss_Krueger_Grid_to_Geo(	east,north : extended; 
				    var phi,lambda : extended; 
				    var zone	   : integer);


Procedure UTM_Geo_to_Grid(    phi,lambda   : extended;
			  var east,north   : extended;
			  var zone	   : integer;
			  var z_designator : char); 
Procedure Geo_to_UTM(	 b,l	    : extended;
		     var east,north : extended;
		     var zone	    : integer);

Implementation

Procedure Trans_Mercator_Geo_to_Grid(	 phi0,lambda0,x0,y0,a,e2,phi,lambda : extended; 
				     var x,y				    : extended);
(****************************************************************************)
(*  ETH Zurich   Geodesy and Geodynamics Lab  (GGL)  , 1996 June 14, (Ge) ***) 
(****************************************************************************)
var  N,M,L,nue2		: extended;
   cp,cp3,cp5,sp,tp,tp2	: extended;
   i			: integer;
   dl			: array[1..6] of extended; c          
			: array[1..6] of extended;
begin			
   dl[1] := lambda-lambda0;
   for i:= 2 to 6 do dl[i]:= dl[i-1]*dl[1];
   M  := Meridian_Radius(a,e2,phi);
   N  := Normal_Radius(a,e2,phi);
   L  := Meridian_length(Phi0,Phi,a,e2);



   nue2 := N/M - 1;
   cp := cos(phi);
   cp3:= cp*cp*cp;
   cp5:= cp3*cp*cp;
   sp := sin(phi);
   tp := sp/cp;
   tp2:= tp*tp;
   
   
   c[1] :=  N    *cp;
   c[3] := -N/6  *cp3*(N/M-tp2);
   c[5] :=  N/120*cp5*(5-18*tp2+tp2*tp2+14*nue2-58*tp2*nue2);
   c[2] := -N/2   *sp*cp;
   c[4] :=  N/24  *sp*cp3*(5    -tp2+9*nue2+4*nue2*nue2);
   c[6] := -N/720 *sp*cp5*(61-58*tp2+tp2*tp2-270*nue2+330*tp2*nue2);
   x  :=     x0 + c[1]*dl[1] - c[3]*dl[3] + c[5]*dl[5] {- c[7]*dl[7]};
   y  := L + y0 - c[2]*dl[2] + c[4]*dl[4] - c[6]*dl[6];
end;

(***************************************************************)

Procedure Trans_Mercator_Grid_to_Geo(	 phi0,lambda0,x0,y0,a,e2,x,y : extended;
				     var phi,lambda		     : extended);
(****************************************************************************) 
(*  ETH Zurich   Geodesy and Geodynamics Lab  (GGL)  , 1996 June 14, (Ge) ***) 
(****************************************************************************)
var  N,M,L,nue2	: extended;
   phii,cp,sp,	
   tp,tp2,tp4,
   cpN,tpN    : extended;
   phi_old    : extended;
   i          : integer;
   dx         :array[1..7] of extended; 
   c          :array[1..7] of extended;
begin
   dx[1] := x-x0;
   

   for i:= 2 to 7 do
      dx[i]:= dx[i-1]*dx[1];
   {Iteration for phii}
   phii := arctan((y-y0)/a)+phi0;
   repeat
      phi_old := phii;
      L    := Meridian_length(Phi0,Phii,a,e2); 
      phii := phi_old + arctan((y-y0-L)/a);
   until abs(phii-phi_old) < 1E-8;
   M  := Meridian_Radius(a,e2,phii);
   N  := Normal_Radius(a,e2,phii);
   nue2 := N/M - 1;
   cp := cos(phii);
   sp := sin(phii);
   tp := sp/cp;
   tp2:= tp*tp;
   tp4:= tp2*tp2;
   cpn:= 1/cp/N;
   tpN:= tp/M/N;
   c[1] := cpN;


   c[3] := cpN/6/N/N            *(N/M +2*tp2);
   c[5] := cpN/120/N/N/N/N      *(5  +28*tp2  +24*tp4);


   c[7] := cpN/5040/N/N/N/N/N/N *(61+662*tp2+1320*tp4+720*tp4*tp2);
   c[2] := tpN/2;
   c[4] := tpN/24/N/N      *(5  +3*tp2+nue2-9*tp2*nue2);
   c[6] := tpN/720/N/N/N/N *(61+90*tp2+45*tp4);
   lambda :=    lambda0 + c[1]*dx[1] - c[3]*dx[3] + c[5]*dx[5] - c[7]*dx[7];
   phi    :=    phii    - c[2]*dx[2] + c[4]*dx[4] - c[6]*dx[6];
end;

(***************************************************************)

Procedure Gauss_Krueger_Geo_to_Grid(	phi,lambda : extended;
				    var east,north : extended;
				    var zone	   : integer);
(****************************************************************************) 
(*  ETH Zurich   Geodesy and Geodynamics Lab  (GGL)  , 1996 June 14, (Ge) ***)
(****************************************************************************)
const
   {BESSEL-ELLIPSOID}		    
   a   = 6377397.155;
   ex2 = 0.006674372231;
var phi_0,North_0,
   dLam,east_0 : extended;
   Lambda_0    : extended;
   i	       : integer;
begin	       
   dLam	:= 3*pi/180;

   
   zone	:= trunc((Lambda+dlam/2)/dlam);
   phi_0	:= 0.0;
   Lambda_0	:= zone*dLam;
   North_0	:= 0.0;


   East_0	:= zone*1e6 + 500000.0;
   Trans_Mercator_Geo_to_Grid(phi_0,Lambda_0,
			      East_0,North_0,
			      Bessel.a,Bessel.ex2,
			      phi,lambda,
			      East,North);
end;

(***************************************************************)

Procedure Gauss_Krueger_Grid_to_Geo(	east,north : extended;
				    var phi,lambda : extended;
				    var zone	   : integer);
(****************************************************************************) 
(*  ETH Zurich   Geodesy and Geodynamics Lab  (GGL)  , 1996 June 14, (Ge) ***)
(****************************************************************************)
var phi_0,North_0,
   dLam,east_0 : extended;
   Lambda_0    : extended;
   i	       : integer;
begin


   dLam	:= 3*pi/180;  {3 Grad-Streifen}
   zone	:= trunc(East/1e6);
   North_0	:= 0.0;
   East_0	:= zone*1e6 + 500000.0;
   phi_0	:= 0.0;
   Lambda_0	:= zone*dLam;



   Trans_Mercator_Grid_to_Geo(phi_0,lambda_0,
			      East_0,North_0, Bessel.a,Bessel.ex2,
			      East,North,
			      phi,lambda);
end; 

(***************************************************************)

Procedure UTM_Geo_to_Grid(    phi,lambda   : extended;
			  var east,north   : extended;
			  var zone	   : integer;
			  var z_designator : char);
const
   {BESSEL-ELLIPSOID}		    
   a   = 6377397.155;
   ex2 = 0.006674372231;
var phi_0,North_0,k0,
   dLam,east_0 : extended;
   Lambda_0    : extended;
   i	       : integer;
begin	       
   
   
   if (phi*180/pi<84) or (phi*180/pi>-80) then
   begin
      {prepare Transverse Mercator Projection} 
      k0	:= 0.9996;
      dLam	:= 6*pi/180;
      zone	:= trunc(Lambda/dlam);
      Lambda_0  := zone*dLam+dlam/2;
      phi_0	:= 0.0;
      zone	:= ((zone+30) mod 60)+1;
      if phi>=	0.0 then  North_0 := 0.0
	 
	 
      else  North_0 := 10000000;
      East_0	:= 500000.0;
      Trans_Mercator_Geo_to_Grid(phi_0,Lambda_0,
				 East_0,North_0, 
				 Bessel.a,Bessel.ex2, phi,lambda,
				 East,North);
   end
   else
   begin
      {prepare Universal Polar Stereographic Projection}
      if phi*180/pi<-80 then
	 if lambda>pi then Z_designator := 'A'
	 else Z_designator := 'B';
      if phi*180/pi> 84 then
	 if lambda>pi then Z_designator := 'Y'
	 else Z_designator := 'Z';
   end;
end;
   (***************************************************************) 
Procedure Geo_to_UTM(	 b,l	    : extended;
		     var east,north : extended;
		     var zone	    : integer);
const
   er	= 6378815.90365;
   pi2	= 1.570796327;
   ald	= 1;
   kf	= 1.0030714396280;
   rho	= 57.2957795;
   a	= 6377397.155{Grosse Halbachse von BESSEL};
   khb	= 6356078.96282{Kleine Halbachse von BESSEL};
   {a	= 6378206.4;
    khb	= 6356583.8;CLARKE}
   eps2	= 0.006674372230614;
   eb	= 2.718281828;
   al	= 1.0007291384304;
   epsk	= 1e-6;
var	
   sgb2,sgb4,sgb6,sgb8,sub1,sub2,bq,pq,cgb2,cgb4,cgb6,l2,l4,sub3, 
   sub5,sub6,sub4,lnull,a0,a1,a2,a3,a4,a5,a6,a8,exzklk,vexz,zexz,sexz,ac, 
   bmin,bmax,lmin,lmax,bmik,bmak,lmik,lmak,balt,lalt, sinphi2,sinphi4,sinphi6, 
   bb,dwdb,x,y,z,k1,r,x1,y1,m,ll,b0,l0,li,la,bi,ba,db,dl,sinphi8,kl1, 
   kl2,kl3,kl4,ck,c,
   f,n0,mv,za,ne,eps,fv,rho1,za1,ne1,cot,tro,sinb0,eps05, 
   ebipm,eb0pm,ebip,ebim,eb0p,eb0m,ebi,eb0,tgi,tg0,n,rho0,tot,toti, 
   ti,t0,m0,teta,tang0,tangi,phi,phi0,phidiff,xs,ys,tgphi2,tgphi02, 
   dll,dlb,xf,yf,zf,uf,cob0,sib0,pi4,ka,kb,kc,kd,ke,kf1,kg,ks, 
   llq,bbq,llq2,llq3,llq4,llq5,llq6,bbq2,
   bbq3,bbq4,bbq5,bbq6,sy11,sy12,sy13, 
   sy14,sy15,sy16,sy31,sy32,sy33,sy34,sy51,sy52,sy1,sy3,sy5,sx01,sx02,sx03, 
   sx04,sx05,sx06,sx0,sx21,sx22,sx23,sx24,sx25,sx41,sx42,sx43,sx2,sx4,sx6 
   :extended; {20}
   i,j,ialt,filter,q,il,ib,sit	:integer;
begin
   {mainprogram utmprojsphaere}
   exzklk:=(sqr(a)-sqr(khb))/sqr(khb);
   ck:=sqr(a)/khb;
   ka:=0.96047;
   kb:=-0.00857;
   kc:=6.41;
   kd:=2.6666;
   ke:=-0.367;
   {kf1:=-0.150;}
   kg:=0.0379;
   ks:=1.0;
   kl1:=1-(0.99*exzklk);
   kl2:=1-(0.984375*exzklk*kl1);
   kl3:=1-(0.9722*exzklk*kl2);
   kl4:=1-(0.9375*exzklk*kl3);
   a0:= 1-(0.75*exzklk*kl4);
   ac:=a0*ck;
   kl1:=1-(1.1413978*exzklk);
   kl2:=1-(1.18892046*exzklk*kl1);
   kl3:=1-(1.28333333*exzklk*kl2);
   kl4:=1-(1.562500*exzklk*kl3);
   a1:= (0.75*exzklk*kl4);
   kl1:=1-(0.984029056*exzklk);
   kl2:=1-(0.977517986*exzklk*kl1);
   kl3:=1-(0.965277778*exzklk*kl2);
   a2:= (0.625*exzklk*kl3);
   
   
   
   kl1:=1-(1.47379333*exzklk);
   kl2:=1-(1.953125*exzklk*kl1);
   a4:= (0.4861111*sqr(exzklk)*kl2);
   kl1:=1-(2.9475*exzklk);
   a6:= (0.41015625*sqr(exzklk)*exzklk*kl1);
   a8:= (0.3609375*sqr(exzklk)*sqr(exzklk));
   vexz:=4*sqr(exzklk);
   zexz:=24-(58*exzklk);
   sexz:=72*exzklk;
   lnull := -75;
   mv:=1/m;
   pi4:=pi/4;
   uf:=1/rho;
   lnull:=lnull/rho;
   lnull := l*rho;
   zone  := trunc((lnull)/6)+31;


   lnull := (6*(zone-31)+3)/rho;
   ll := l-lnull;
   {****}
   sinphi2:=sqr(sin(b));
   sinphi4:=sqr(sinphi2);
   sinphi6:=sinphi2*sinphi4;
   sinphi8:=sqr(sinphi4);
   bq:=ac*b-(a1*ck*sin(b)*cos(b)*(1+(a2*sinphi2)+(a4*sinphi4)+(a6*sinphi6) 
	     +(a8*sinphi8)));
   {bq=Meridianbogenlaenge}
   pq:=ck/(sqrt(sqr(1/cos(b))+exzklk));
   {pq=Radius des Parallelkreises P=n*cos phi}
   cgb2:=cos(b)*cos(b);
   cgb4:=cgb2*cgb2;
   cgb6:=cgb4*cgb2;
   l2:=ll*ll;
   l4:=l2*l2; sub3:=1+(l2/12)*(-
			       1+6*cgb2+(9*exzklk)*cgb4+vexz*cgb6)+(l4/360)*(1-60*cgb2+120*cgb4); 
   {kf1:=1/er;}
   sub5  := exzklk*cgb4;
   sub6  := 2*cgb2;
   sub4  := 1+((l2/6)*(-1+sub6+sub5))+((l4/120)*(1-
				       20*cgb2+zexz*cgb4+sexz*cgb6)); North :=  
   0.9996*(bq+(((pq*l2)/2)*sin(b)*sub3));
   East  :=  0.9996*(pq*ll*sub4)+500000;
   {abb.gl. ,utmprojsphaere}
end;
end; 
   (***************************************************************) 
begin
end.
