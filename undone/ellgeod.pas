(****************************************************************************)
   Unit EllGeod;                                    (*  21.June 1996 ³      *)
                                                    (************************)
                                                    (*  ETH ZUERICH GGL-IGP *)
                                                    (*          (Ge)        *)
(****************************************************************************)
Interface
Function Meridian_Radius(a,e2,phi:extended):extended;
Function Normal_Radius(a,e2,phi:extended):extended;
Function Meridian_length(Phi1,Phi2,a,e2:extended):extended;
Function ARCSIN(X:EXTENDED):EXTENDED;
Function ARCTAN2 ( SI, CO : EXTENDED): EXTENDED;
(****************************************************************************)

Implementation

Function Meridian_Radius(a,e2,phi:extended):extended;
(****************************************************************************)
(*  ETH Zurich   Geodesy and Geodynamics Lab  (GGL)  , 1996 June 14, (Ge) ***)
(****************************************************************************)
var   w,s   : extended;

begin
   s  := sin(phi);
   w  := sqrt(1-e2*s*s);
   Meridian_Radius := a*(1-e2)/w/w/w;
end;
(***************************************************************)
Function Normal_Radius(a,e2,phi:extended):extended;
(****************************************************************************)
(*  ETH Zurich   Geodesy and Geodynamics Lab  (GGL)  , 1996 June 14, (Ge) ***)
(****************************************************************************)
var   w,s   : extended;

begin
   s  := sin(phi);
   w  := sqrt(1-e2*s*s);
   Normal_Radius := a/w;
end;

(***************************************************************)
Function Meridian_length(Phi1,Phi2,a,e2:extended):extended;
(****************************************************************************)
(*  ETH Zurich   Geodesy and Geodynamics Lab  (GGL)  , 1996 June 14, (Ge) ***)
(****************************************************************************)
var dp,
    dphi,phi,sum  : extended;
    i,n_dphi      : longint;

begin
  dphi :=  2/60*pi/180;
  dp := phi2-phi1;
  n_dphi := trunc(dp/dphi);
  dphi   := dp/n_dphi;
  sum := 0.0;
  for i := 1 to n_dphi do
  begin
     phi := phi1+i*dphi-dphi/2;
     sum := sum + dphi*Meridian_Radius(a,e2,phi);
  end;
  Meridian_length := sum;
end;


(****************************************************************************)
    FUNCTION ARCSIN(X:EXTENDED):EXTENDED;
(****************************************************************************)
    BEGIN
       ARCSIN:=ARCTAN(X/(SQRT((1-X*X))+1.0E-15));
    END;

(****************************************************************************)
   FUNCTION ARCTAN2 ( SI, CO : EXTENDED): EXTENDED;
(****************************************************************************)
(* Berechnung des Arcus - tangens QUADRANTENRICHTIG  in RADIAN  [0, 2*pi]   *)
(* si = faktor * sin(Winkel) und co = faktor * cos(Winkel)                  *)
(* ---> (si, co) muss nicht normiert sein                                   *)
(****************************************************************************)
   VAR   X : EXTENDED;
   BEGIN
     IF ABS(CO)>1.0E-12 THEN
     BEGIN
       X:=SI/CO;
       X:=ARCTAN(X);
       if CO < 0.0  THEN X := X + PI;
       if X  < 0.0  THEN X := X + 2*PI;
       ARCTAN2:=X;
     END
     ELSE
     BEGIN
       IF SI<0.0 THEN ARCTAN2:=3/2*PI
                 ELSE ARCTAN2:=PI/2;
     END;
   END;


end.