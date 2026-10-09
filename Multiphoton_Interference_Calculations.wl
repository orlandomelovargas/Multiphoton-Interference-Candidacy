(* ::Package:: *)

(* Companion computational file for

   "Partial Indistinguishability and Sector Structure in
    Multiphoton Interference: From the Effective
    Single-Particle Matrix to the Semiclassical Limit"

   PhD Candidacy Examination, 2026
   Stiven Orlando Melo Vargas

   Mathematical derivations: Mathematical Appendices A--I
   Repository: https://github.com/orlandomelovargas/Multiphoton-Interference-Candidacy.git
   Release: candidacy-2026

   Tested with Wolfram Mathematica 14.5 at machine precision.

   The reference SU(3) numerical configurations use
   overline(Lambda) = lambda = spec(rho).
   This equality is specific to those configurations and is not
   assumed in the general theoretical formulation.
*)


ClearAll["Global`M3*"];
(* sharper multiplicity numbers in the weight diagram. *)
$M3BaseStyle = {FontFamily -> "Times", FontSize -> 12};
$M3WeightMapCache = <||>;

(* Safe power for nonnegative integer exponents *)
M3SafePow[z_, n_Integer?NonNegative] := If[n == 0, 1, z^n];

(* 2x2 determinant. *)
M3Det2[A_?MatrixQ] := Module[{},
  If[Dimensions[A] =!= {2,2},
    Return[Failure["Invalid2x2Matrix", <|"Message" -> "M3Det2 requires a 2x2 matrix."|>]]
  ];
  A[[1,1]] A[[2,2]] - A[[1,2]] A[[2,1]]
];
$M3SectorBlockCache = <||>;
$M3CompositionsCache = <||>;
$M3CompositionIndexCache = <||>;
$M3QBCache = <||>;
(* Keep Dynamic evaluations active in the notebook. *)
If[$FrontEnd =!= Null,
  Quiet @ Check[SetOptions[EvaluationNotebook[], DynamicUpdating -> True], Null]
];

M3PlotLabel[s_] := Style[s, 13, Bold, FontFamily -> "Times"];
M3FrameLabel[s_] := Style[s, 11.5, FontFamily -> "Times"];
M3ScientificMatrixPlot[A_?MatrixQ,title_,labels_:Automatic] := Module[{m,n,cf,vmax,ticksLeft,ticksBottom},
  m = N[Abs[A]];
  n = Length[m];
  vmax = Max[10^-12,Max[m]];
  cf = Function[z,Blend[
    {RGBColor[.99,.99,1.],RGBColor[.82,.88,.94],RGBColor[.42,.60,.77],RGBColor[.12,.29,.47]},z
  ]];
  ticksLeft = If[
    labels === Automatic || Length[labels] =!= n,
    Automatic,
    Table[{i,labels[[n-i+1]]},{i,n}]
  ];
  ticksBottom = If[
    labels === Automatic || Length[labels] =!= n,
    Automatic,
    Table[{i,labels[[i]]},{i,n}]
  ];
  ArrayPlot[
    Reverse[m],
    ColorFunctionScaling -> True,
    ColorFunction -> cf,
    Mesh -> All,
    MeshStyle -> Directive[White,AbsoluteThickness[.55],Opacity[.8]],
    Frame -> True,
    FrameTicks -> {{ticksLeft,None},{ticksBottom,None}},
    PlotLabel -> M3PlotLabel[title],
    PlotLegends -> Placed[
      BarLegend[{cf,{0,vmax}},LegendLabel -> Style["|element|",9,FontFamily -> "Times"]],
      Right
    ],
    Background -> White,
    ImagePadding -> {{48,55},{42,25}},
    BaseStyle -> {FontFamily -> "Times",FontSize -> 11},
    ImageSize -> 430
  ]
];


(* APPENDICES D/E/H: NORMALIZED DATA AND FINITE REPRESENTATION *)
M3NormalizeExact[v_List] := v/Total[v];
(* Functions used directly by the interface. *)
M3LambdaFromIntegerRatios[a_Integer,b_Integer,c_Integer] :=
  Reverse @ Sort @ M3NormalizeExact[{a,b,c}];

M3NbarFromIntegerRatios[a_Integer,b_Integer,c_Integer] :=
  M3NormalizeExact[{a,b,c}];

M3RationalVector[v_List,tol_:10^-10] := Rationalize[N[v],tol];

M3MinimalCompatibleN[lam_List,n_List,tol_:10^-10] := Module[{r},
  r = Join[M3RationalVector[lam,tol],M3RationalVector[n,tol]];
  Apply[LCM,Denominator /@ r]
];

M3MajorizedQ3[lam_List,n_List,tol_:10^-10] := Module[{ls,ns},
  ls = Reverse @ Sort @ N[lam];
  ns = Reverse @ Sort @ N[n];
  Abs[Total[ns]-Total[ls]] <= tol &&
  ns[[1]] <= ls[[1]] + tol &&
  Total[ns[[1;;2]]] <= Total[ls[[1;;2]]] + tol
];

(*Eqs. (D.2)-(D.9), with normalized data as in Eq. (H.12). *)
M3FiniteRepresentationData[lam_List,n_List,Nrep_Integer?Positive,tol_:10^-10] := Module[
  {lr,nr,lamN0,nN0,lamN,nN,nmin,p,q,r,W},
  lr = M3RationalVector[lam,tol];
  nr = M3RationalVector[n,tol];
  nmin = M3MinimalCompatibleN[lr,nr,tol];
  lamN0 = Nrep lr;
  nN0 = Nrep nr;
  If[!VectorQ[lamN0,IntegerQ] || !VectorQ[nN0,IntegerQ],
    Return[Failure["IncompatibleScale",<|"Message"->"N does not make N lambda and N n integers.","SuggestedN"->nmin|>]]
  ];
  lamN = lamN0; nN = nN0;
  If[Total[lamN]=!=Nrep || Total[nN]=!=Nrep || !(lamN[[1]]>=lamN[[2]]>=lamN[[3]]>=0) || Min[nN]<0,
    Return[Failure["InvalidFiniteData",<|"Message"->"The finite data do not form a valid partition/weight."|>]]
  ];
  p=lamN[[1]]-lamN[[2]]; q=lamN[[2]]-lamN[[3]]; r=lamN[[3]]; W=nN-r {1,1,1};
  <|"N"->Nrep,"LambdaN"->lamN,"nN"->nN,"p"->p,"q"->q,"r"->r,"W"->W,
    "Dynkin"->{p,q},"IrrepDimension"->((p+1)(q+1)(p+q+2))/2,"MinimalCompatibleN"->nmin|>
];


(* APPENDIX E \[LongDash] ISOSPIN AND MULTIPLICITY Eqs. (E.25)\[Dash](E.29), (E.38), (E.41). *)
M3IsospinData[data_Association] := Module[{p,q,W,m,deltaJ,jmin,jmax,jvals},
  p = data["p"];
  q = data["q"];
  W = data["W"];
  m = (W[[1]] - W[[2]])/2;
  deltaJ = (p - W[[3]])/2;
  jmin = 1/2 Max[Abs[W[[1]] - W[[2]]],Abs[p - W[[3]]]];
  jmax = 1/2 Min[p + W[[3]],p + 2 q - W[[3]]];
  jvals = If[jmax < jmin,{},Range[jmin,jmax,1]];
  <|"m" -> m,"DeltaJ" -> deltaJ,"JMin" -> jmin,"JMax" -> jmax,"JValues" -> jvals,"Multiplicity" -> Length[jvals]|>
];


(* APPENDIX B: BISTOCHASTIC AND UNISTOCHASTIC GEOMETRY OF SU(3) Eqs. (B.3)-(B.10) and (B.20)-(B.26).*)
M3SectionUV[lam_List,n_List,x_?NumericQ,y_?NumericQ] := Module[{d1,d2},
  d1=lam[[1]]-lam[[3]]; d2=lam[[2]]-lam[[3]];
  If[d2==0,Return[Failure["DegenerateSpectrum",<|"Message"->"The (x,y) chart is singular for lambda2 = lambda3."|>]]];
  {(n[[1]]-lam[[3]]-d1 x)/d2,(n[[2]]-lam[[3]]-d1 y)/d2}
];

M3QSection[lam_List,n_List,x_?NumericQ,y_?NumericQ] := Module[{uv,u,v},
  uv=M3SectionUV[lam,n,x,y]; If[FailureQ[uv],Return[uv]]; {u,v}=uv;
  {{x,u,1-x-u},{y,v,1-y-v},{1-x-y,1-u-v,x+y+u+v-1}}
];

M3BistochasticSectionQ[Q_?MatrixQ,tol_:10^-10] := Min[Flatten[N[Q]]]>=-tol;

M3TriangleLengths[Q_?MatrixQ] := Module[{xi},
  xi=N[Q[[All,1]]Q[[All,2]]]; Sqrt[Map[Max[0.,#]&,xi]]
];

M3TriangleInequalitiesQ[L_List,tol_:10^-10] := Max[L]<=Total[L]-Max[L]+tol;

M3HeronPolynomialQ[Q_?MatrixQ] := Module[{xi},
  xi=Q[[All,1]]Q[[All,2]];
  2(xi[[1]]xi[[2]]+xi[[1]]xi[[3]]+xi[[2]]xi[[3]])-Total[xi^2]
];

M3FrameAdmissiblePointQ[lam_List,n_List,x_?NumericQ,y_?NumericQ,tol_:10^-10] := Module[{Q},
  Q=M3QSection[lam,n,x,y]; If[FailureQ[Q],Return[False]];
  M3BistochasticSectionQ[Q,tol]&&M3HeronPolynomialQ[Q]>=-tol
];

M3PhaseSolution3[L_List,branch_:1,tol_:10^-12] := Module[{l1,l2,l3,c,phi2,phi3,sgn},
  {l1,l2,l3}=N[L]; sgn=If[branch>=0,1,-1];
  If[!M3TriangleInequalitiesQ[{l1,l2,l3},tol],
    Return[Failure["NoTriangle",<|"Message"->"The sides do not close a triangle."|>]]];
  Which[
    l1<tol&&l2<tol&&l3<tol,{0.,0.,0.},
    l1<tol,{0.,0.,Pi},
    l2<tol,{0.,0.,Pi},
    l3<tol,{0.,sgn Pi,0.},
    True,
      c=Clip[(l3^2-l1^2-l2^2)/(2 l1 l2),{-1,1}];
      phi2=sgn ArcCos[c]; phi3=Arg[-(l1+l2 Exp[I phi2])];
      {0.,phi2,phi3}
  ]
];

M3ReconstructUnitary3[Q_?MatrixQ,branch_:1,tol_:10^-10] := Module[{L,ph,c1,c2,c3,cross,V,detPhase},
  If[Min[Flatten[N[Q]]]<-tol,Return[Failure["InvalidQ",<|"Message"->"Q contains negative entries."|>]]];
  L=M3TriangleLengths[Q]; ph=M3PhaseSolution3[L,branch,tol]; If[FailureQ[ph],Return[ph]];
  c1=Sqrt[Map[Max[0.,#]&,N[Q[[All,1]]]]];
  c2=Sqrt[Map[Max[0.,#]&,N[Q[[All,2]]]]]Exp[I ph];
  cross=Conjugate[Cross[c1,c2]];
  If[Norm[cross]<tol,Return[Failure["DegenerateColumns",<|"Message"->"The third column could not be constructed."|>]]];
  c3=cross/Norm[cross]; V=Transpose[{c1,c2,c3}];
  detPhase=Det[V]/Abs[Det[V]]; V[[All,3]]*=Conjugate[detPhase]; Chop[V,tol]
];

(* reconstruction of rho and Gamma. *)
M3RhoGammaFromFrame[lam_List,n_List,V_?MatrixQ,tol_:10^-12] := Module[{rho,Gamma,Dinv},
  If[Min[n]<=0,Return[Failure["ZeroOccupation",<|"Message"->"n must be strictly positive to construct Gamma."|>]]];
  rho=Chop[V . DiagonalMatrix[lam] . ConjugateTranspose[V],tol];
  Dinv=DiagonalMatrix[1/Sqrt[n]]; Gamma=Chop[Dinv . rho . Dinv,tol];
  <|"Rho"->rho,"Gamma"->Gamma|>
];


(* NUMERICAL SAMPLING OF PHYSICALLY ADMISSIBLE FRAMES *)
M3FarthestPointSample[pts_List,k_Integer?Positive] := Module[{center,remaining,selected,first,pick},
  If[pts==={},Return[{}]];
  center=Mean[N[pts]]; first=First@MinimalBy[pts,Norm[N[#]-center]&];
  selected={first}; remaining=DeleteCases[pts,first];
  While[Length[selected]<Min[k,Length[pts]]&&remaining=!={},
    pick=First@MaximalBy[remaining,Function[p,Min[Norm[N[p-#]]&/@selected]]];
    AppendTo[selected,pick]; remaining=DeleteCases[remaining,pick]];
  selected
];

M3AdmissibleGridPoints[lam_List,n_List,gridN_Integer?Positive] := Module[{grid},
  grid=N@Subdivide[0,1,gridN]; Select[Tuples[grid,2],M3FrameAdmissiblePointQ[lam,n,#[[1]],#[[2]]]&]
];


(* APPENDIX E: ISOSPIN AND WIGNER-D COORDINATES Combinatorial weight Q(a3|W,j): Eq. (E.46).*)
M3QWeight[p_Integer?NonNegative,q_Integer?NonNegative,W_List,j_,a3_Integer?NonNegative] := Module[{b3,s,args},
  b3=a3+q-W[[3]]; s=(p+W[[3]])/2-a3; args={a3,b3,s-j,s+j+1};
  If[!And@@(IntegerQ[#]&&#>=0&/@args),0,Factorial[p+q]/Times@@(Factorial/@args)]
];

(* Normalization N(W,j), direct sum: Eq. (E.47). *)
M3NWeightSum[p_Integer?NonNegative,q_Integer?NonNegative,W_List,j_] := Total[M3QWeight[p,q,W,j,#]&/@Range[0,p]];

(* Closed normalization N(W,j): Eq. (E.49). *)
M3NWeight[p_Integer?NonNegative,q_Integer?NonNegative,W_List,j_] := Module[{u,v,args},
  u=(p+W[[3]])/2; v=(p+2 q-W[[3]])/2; args={u-j,u+j+1,v-j,v+j+1};
  If[!And@@(IntegerQ[#]&&#>=0&/@args),0,Factorial[p+q]Factorial[p+q+1]/Times@@(Factorial/@args)]
];


(* Polynomial continuation of Wigner-D for a complex 2x2 matrix: Eqs. (E.57)-(E.60). *)
(*Appendix E.3, polynomial representation of D^j_{m,m'}(g). *)
M3WignerDPolynomial[j_,m_,mp_,g_?MatrixQ] := Module[{kmin,kmax,norm,validLabels},
  If[Dimensions[g]=!={2,2},Return[Failure["InvalidWignerMatrix",<|"Message"->"The matrix g used in Wigner-D must be 2x2."|>]]];
  validLabels=And[IntegerQ[2 j],IntegerQ[j+m],IntegerQ[j-m],IntegerQ[j+mp],IntegerQ[j-mp],j+m>=0,j-m>=0,j+mp>=0,j-mp>=0];
  If[!validLabels,Return[0]];
  kmin=Ceiling[Max[0,mp-m]]; kmax=Floor[Min[j+mp,j-m]];
  If[kmin>kmax,Return[0]];
  norm=Sqrt[Binomial[2 j,j+m]Binomial[2 j,j+mp]];
  If[PossibleZeroQ[norm],Return[0]];
  Total[Table[
    With[{a11=j+mp-k,a12=m-mp+k,a21=k,a22=j-m-k},
      Factorial[2 j]/(Factorial[a11]Factorial[a12]Factorial[a21]Factorial[a22])*
      M3SafePow[g[[1,1]],a11]*M3SafePow[g[[1,2]],a12]*
      M3SafePow[g[[2,1]],a21]*M3SafePow[g[[2,2]],a22]],
    {k,kmin,kmax}]]/norm
];

(* Local matrix g(V): Eqs. (E.54)-(E.56). *)
M3GMatrix[V_?MatrixQ] := Module[{zeta,eta},
  If[Dimensions[V]=!={3,3},Return[Failure["InvalidFrame",<|"Message"->"The frame V must be a 3x3 matrix."|>]]];
  zeta=V[[All,1]]; eta=Conjugate[V[[All,3]]];
  {{zeta[[1]],eta[[2]]},{zeta[[2]],-eta[[1]]}}
];

(* Wigner-D coordinates of the projected state: Eqs. (E.65)-(E.68) and (E.73)-(E.81). *)
M3WignerCoordinates[finite_Association,iso_Association,V_?MatrixQ,tol_:10^-12] := Module[
  {p,q,W,m,deltaJ,jvals,g,detg,gTilde,raw,norm,coords,probs},
  p=finite["p"]; q=finite["q"]; W=finite["W"]; m=iso["m"]; deltaJ=iso["DeltaJ"]; jvals=iso["JValues"];
  If[jvals==={},Return[Failure["EmptyIsospinSector",<|"Message"->"The weight contains no allowed values of j."|>]]];
  g=M3GMatrix[V]; If[FailureQ[g],Return[g]];
  detg=M3Det2[g]; If[FailureQ[detg],Return[detg]];
  detg=Chop[N[detg],tol];

  (* det(g)=0 does NOT mean that the physical state is singular. It only means that this local Wigner chart is no longer valid. *)
  If[Abs[N[detg]]<=tol,Return[Failure["WignerChartSingular",<|"Message"->"det(g) is close to zero. Another chart or the exact ambient construction must be used."|>]]];

  (* Eq. (E.56). *)
  gTilde=g/Sqrt[detg];

  (* Unnormalized coordinates A_j(V). *)
  raw=Table[With[{Nj=M3NWeight[p,q,W,j]},Sqrt[N[(2 j+1)Nj]]*M3WignerDPolynomial[j,m,deltaJ,gTilde]],{j,jvals}];
  norm=Norm[N[raw]];
  If[!NumericQ[norm]||norm<=tol,Return[Failure["ZeroWignerNorm",<|"Message"->"The Wigner coordinate vector has numerically zero norm."|>]]];
  coords=Chop[N[raw/norm],tol]; probs=Chop[Abs[coords]^2,tol];

  <|"g"->g,"DetG"->detg,"gTilde"->gTilde,"JValues"->jvals,"RawCoordinates"->raw,
    (* Z(W,gTilde) = Sum_j |A_j|^2 *)
    "Z"->norm^2,"CoordinateNorm"->norm,"Coordinates"->coords,
    (* Eq. (E.79). *)
    "OrbitalProbabilities"->probs,
    (* Eq. (E.81), normalization diagnostic. *)
    "ProbabilitySum"->Total[probs]|>
];
(* Orbital probability P_orb(j|W,V): Eqs. (E.79)-(E.81). *)


(* APPENDIX F: COORDINATES, KERNEL, GRAM MATRIX, AND RANK; Eqs. (F.20)-(F.48). *)
(* Kernel between two frames in isospin coordinates: Eqs. (F.20)-(F.25). *)
(* Construction of matrix C: Eq. (F.35). *)
M3CoordinateMatrix[coordinates_List] := Module[{dims},
  If[coordinates==={},Return[{}]];
  dims=Length/@coordinates;
  If[Length[DeleteDuplicates[dims]]=!=1,Return[Failure["IncompatibleCoordinates",<|"Message"->"All coordinate vectors must have the same dimension."|>]]];
  Transpose[coordinates]
];

(* Gram matrix: Eqs. (F.33)-(F.40), in particular G = C^dagger C in Eq. (F.36). *)
M3GramMatrix[C_?MatrixQ,tol_:10^-12] := Chop[ConjugateTranspose[N[C]] . N[C],tol];

(* Eigenvalues of the Gram matrix *)
M3GramEigenvalues[G_?MatrixQ] := Sort[Re[Eigenvalues[N[(G+ConjugateTranspose[G])/2]]],Greater];

(* Singular values of C *)
M3CoordinateSingularValues[C_?MatrixQ] := SingularValueList[N[C]];

(* Numerical rank of C: Eqs. (F.41)-(F.48). *)
M3CoordinateRank[C_?MatrixQ,tol_:10^-8] := Module[{sv,scale},
  sv=M3CoordinateSingularValues[C]; If[sv==={},Return[0]];
  scale=Max[sv]; If[scale==0,Return[0]];
  Count[sv,x_/;x>tol scale]
];

(* Rank of G as an independent diagnostic *)
M3GramRank[G_?MatrixQ,tol_:10^-8] := Module[{ev,scale},
  ev=M3GramEigenvalues[G]; If[ev==={},Return[0]];
  scale=Max[Abs[ev]]; If[scale==0,Return[0]];
  Count[ev,x_/;x>tol scale]
];

(* Incremental selection of linearly independent frames *)
M3IndependentFrameIndices[C_?MatrixQ,tol_:10^-8] := Module[{basis={},indices={},v,r,scale},
  Do[
    v=N[C[[All,i]]]; r=v;
    (* First orthogonalization pass. *)
    Do[r=r-(Conjugate[b] . r)b,{b,basis}];
    (* Second pass to improve numerical stability. *)
    Do[r=r-(Conjugate[b] . r)b,{b,basis}];
    scale=Max[1.,Norm[v]];
    If[Norm[r]>tol scale,AppendTo[basis,r/Norm[r]];AppendTo[indices,i]],
    {i,Dimensions[C][[2]]}
  ];
  indices
];

(* Compact Gram / rank diagnostic *)
M3GramDiagnostics[C_?MatrixQ,tol_:10^-8] := Module[{G,sv,ev,rankC,rankG,inds,scaleC,scaleG},
  G=M3GramMatrix[C]; sv=M3CoordinateSingularValues[C]; ev=M3GramEigenvalues[G];
  scaleC=If[sv==={},0,Max[sv]];
  rankC=If[scaleC==0,0,Count[sv,x_/;x>tol scaleC]];
  scaleG=If[ev==={},0,Max[Abs[ev]]];
  rankG=If[scaleG==0,0,Count[ev,x_/;x>tol scaleG]];
  inds=M3IndependentFrameIndices[C,tol];
  <|"Gram"->G,"GramEigenvalues"->ev,"SingularValues"->sv,"RankC"->rankC,"RankG"->rankG,
    "RankAgreement"->(rankC==rankG),"IndependentFrameIndices"->inds|>
];


(* APPENDIX H: STATIONARY POINT, HESSIAN, AND TOROIDAL INTEGRAL; Eqs. (H.9)-(H.69). *)
(* Residual of the stationary condition: Eqs. (H.25)-(H.26). *)
M3StationaryResidual[lam_List,n_List,Q_?MatrixQ] := Module[{res},res=Q . lam-n;Chop[N[res],10^-13]];

(* Complete intensive Hessian data: Eqs. (H.43)-(H.59). *)
M3HessianData[lam_List,Q_?MatrixQ] := Module[{alpha1,alpha2,a,b,w12,w13,w23,H,Hred,detH},
  If[Dimensions[Q]=!={3,3},Return[Failure["InvalidQ",<|"Message"->"Q must be a 3x3 matrix."|>]]];
  alpha1=lam[[1]]-lam[[2]]; alpha2=lam[[2]]-lam[[3]];

  (* Eq. (H.48), with the first and third columns entering Eq. (H.54). *)
  a=Q[[All,1]]; b=Q[[All,3]];

  (* Eq. (H.54), in intensive scaling. *)
  w12=alpha1 a[[1]]a[[2]]+alpha2 b[[1]]b[[2]];
  w13=alpha1 a[[1]]a[[3]]+alpha2 b[[1]]b[[3]];
  w23=alpha1 a[[2]]a[[3]]+alpha2 b[[2]]b[[3]];

  (* Laplacian form of the Hessian. *)
  H={{w12+w13,-w12,-w13},{-w12,w12+w23,-w23},{-w13,-w23,w13+w23}};

  (* Set phi1 = 0. *)
  Hred={{w12+w23,-w23},{-w23,w13+w23}};

  (* Eq. (H.59). *)
  detH=w12 w13+w12 w23+w13 w23;

  <|"Alpha1"->alpha1,"Alpha2"->alpha2,"a"->a,"b"->b,
    "Weights"-><|"w12"->w12,"w13"->w13,"w23"->w23|>,
    "Hessian"->Chop[N[H],10^-13],"ReducedHessian"->Chop[N[Hred],10^-13],
    "ReducedDeterminant"->Chop[N[detH],10^-13]|>
];

(* Intensive toroidal action: Eq. (H.14), specialized to the SU(3) form of Eqs. (D.48)-(D.51). *)
M3TorusActionIntensive[finite_Association,Q_?MatrixQ,phi2_?NumericQ,phi3_?NumericQ] := Module[
  {Nrep,p,q,W,a,b,phi,delta,z1,z2},
  Nrep=finite["N"]; p=finite["p"]; q=finite["q"]; W=finite["W"];
  a=N[Q[[All,1]]]; b=N[Q[[All,3]]]; phi={0.,phi2,phi3};

  (* (W-q 1)/N for the intensive action. *)
  delta=N[(W-q{1,1,1})/Nrep];
  z1=Total[a Exp[I phi]]; z2=Total[b Exp[-I phi]];
  -I delta . phi+(p/Nrep)Log[z1]+(q/Nrep)Log[z2]
];


(* Numerically evaluated toroidal integral: Eqs. (D.49)-(D.52) and (H.9)-(H.14). *)
(* Gaussian stationary approximation: Eqs. (H.64)-(H.69), with M = 3 in Eq. (H.69). *)
(* APPENDICES E/G: EXACT BASIS |W,j> AND PHYSICAL MATRIX sigma_W; Eqs. (E.45)-(E.53), (G.23), and (G.67)-(G.82). *)
(* Compositions of N into three occupations *)
M3Compositions3[N_Integer?NonNegative] := Module[{cached,comps},
  cached=Lookup[$M3CompositionsCache,N,Missing["NotCached"]];
  If[ListQ[cached],Return[cached]];
  comps=Flatten[Table[{a,b,N-a-b},{a,0,N},{b,0,N-a}],1];
  AssociateTo[$M3CompositionsCache,N->comps];
  comps
];

(* Compact key for occupation associations. *)
M3WtKey[v_List] := M3Wt@@v;

M3CompositionIndex3[N_Integer?NonNegative] := Module[{cached,comps,index},
  cached=Lookup[$M3CompositionIndexCache,N,Missing["NotCached"]];
  If[AssociationQ[cached],Return[cached]];
  comps=M3Compositions3[N];
  index=AssociationThread[M3WtKey/@comps,Range[Length[comps]]];
  AssociateTo[$M3CompositionIndexCache,N->index];
  index
];

(* Cofactor matrix: Eqs. (G.67)-(G.69). *)
M3CofactorMatrix3[A_?MatrixQ] := Module[{a},
  If[Dimensions[A]=!={3,3},Return[Failure["InvalidMatrix",<|"Message"->"A must be a 3x3 matrix."|>]]]; a=A;
  {{a[[2,2]]a[[3,3]]-a[[2,3]]a[[3,2]],-(a[[2,1]]a[[3,3]]-a[[2,3]]a[[3,1]]),a[[2,1]]a[[3,2]]-a[[2,2]]a[[3,1]]},
   {-(a[[1,2]]a[[3,3]]-a[[1,3]]a[[3,2]]),a[[1,1]]a[[3,3]]-a[[1,3]]a[[3,1]],-(a[[1,1]]a[[3,2]]-a[[1,2]]a[[3,1]])},
   {a[[1,2]]a[[2,3]]-a[[1,3]]a[[2,2]],-(a[[1,1]]a[[2,3]]-a[[1,3]]a[[2,1]]),a[[1,1]]a[[2,2]]-a[[1,2]]a[[2,1]]}}
];

(* Full column of Sym^N(A): Eqs. (G.73)-(G.74). *)
M3SymmetricPowerColumn[A_?MatrixQ,in_List] := Module[
  {Nrep,An,degree,comps,vec,nextComps,nextIndex,next,j,t,idx,amp,occ,occ2,i,idx2},
  If[Dimensions[A]=!={3,3}||Length[in]=!=3||Min[in]<0||!VectorQ[in,IntegerQ],
    Return[Failure["InvalidSymmetricPowerInput",<|"Message"->"A must be 3x3 and in must be a nonnegative integer occupation."|>]]];
  Nrep=Total[in]; An=N[A]; degree=0; comps={{0,0,0}}; vec={1.+0.I};
  Do[Do[
    nextComps=M3Compositions3[degree+1]; nextIndex=M3CompositionIndex3[degree+1];
    next=ConstantArray[0.+0.I,Length[nextComps]];
    Do[
      amp=vec[[idx]];
      If[amp=!=0,
        occ=comps[[idx]];
        Do[occ2=ReplacePart[occ,i->occ[[i]]+1]; idx2=Lookup[nextIndex,M3WtKey[occ2],0];
          If[idx2>0,next[[idx2]]+=amp*An[[i,j]]*Sqrt[occ[[i]]+1]/Sqrt[t]],{i,3}]
      ],{idx,Length[vec]}];
    vec=next; comps=nextComps; degree++,
    {t,1,in[[j]]}],{j,3}];
  If[degree=!=Nrep,Return[Failure["SymmetricPowerDegreeError",<|"Message"->"The recurrence did not reach the expected degree."|>]]];
  Chop[vec,10^-14]
];


(* Individual element of Sym^N(A): Eq. (G.74). *)
(* Ambient pairs compatible with weight W: Eqs. (D.34)-(D.35) and (G.75)-(G.76). *)
M3WeightOccupationPairs[p_Integer?NonNegative,q_Integer?NonNegative,W_List] := Module[{delta},
  delta=W-q{1,1,1};
  Cases[M3Compositions3[p],k_:>With[{l=k-delta},If[VectorQ[l,IntegerQ]&&Min[l]>=0&&Total[l]==q,<|"K"->k,"L"->l|>,Nothing]]]
];

(* Safe Clebsch-Gordan evaluation. Avoids evaluating 3j symbols outside the triangle rules. Mathematically, those coefficients are exactly zero. *)
M3ClebschGordanSafe[j1_,m1_,j2_,m2_,j_,m_,tol_:10^-12] := Module[{valid},
  valid = And[
    And@@(IntegerQ[2 #]& /@ {j1,j2,j,m1,m2,m}),
    j1>=0, j2>=0, j>=0,
    Abs[m1]<=j1, Abs[m2]<=j2, Abs[m]<=j,
    TrueQ[m1+m2==m],
    Abs[j1-j2]<=j<=j1+j2,
    IntegerQ[j1+j2+j],
    IntegerQ[j1+m1], IntegerQ[j1-m1],
    IntegerQ[j2+m2], IntegerQ[j2-m2],
    IntegerQ[j+m], IntegerQ[j-m]
  ];
  If[!TrueQ[valid], 0., Quiet@Check[Chop[N[ClebschGordan[{j1,m1},{j2,m2},{j,m}]],tol],0.]]
];

(* Physical isospin basis B^(W): Eqs. (E.45)-(E.53) and (G.77). *)
M3IsospinBasisData[finite_Association,iso_Association,tol_:10^-10] := Module[
  {p,q,r,W,m,jvals,pairs,nWeights,basis,k,l,j1,m1,j2,m2,cg,qw,coeff,gram,err},
  p=finite["p"]; q=finite["q"]; r=finite["r"]; W=finite["W"]; m=iso["m"]; jvals=iso["JValues"];
  If[jvals==={},Return[Failure["ZeroMultiplicity",<|"Message"->"The weight does not appear in the irrep."|>]]];
  pairs=M3WeightOccupationPairs[p,q,W];
  If[pairs==={},Return[Failure["NoWeightPairs",<|"Message"->"There are no ambient pairs compatible with W."|>]]];
  nWeights=AssociationThread[jvals,M3NWeight[p,q,W,#]&/@jvals];
  basis=Transpose[Table[Table[
    k=pair["K"]; l=pair["L"]; j1=(k[[1]]+k[[2]])/2; m1=(k[[1]]-k[[2]])/2;
    j2=(l[[1]]+l[[2]])/2; m2=(l[[2]]-l[[1]])/2;
    cg=M3ClebschGordanSafe[j1,m1,j2,m2,j,m,tol]; qw=M3QWeight[p,q,W,j,k[[3]]];
    coeff=If[qw>0&&nWeights[j]>0,(-1)^l[[1]]*Sqrt[N[qw/nWeights[j]]]*cg,0.]; coeff,
    {pair,pairs}],{j,jvals}]];
  gram=Chop[ConjugateTranspose[basis] . basis,tol]; err=Max[Abs[Flatten[N[gram-IdentityMatrix[Length[jvals]]]]]];
  <|"p"->p,"q"->q,"r"->r,"W"->W,"JValues"->jvals,"OccupationPairs"->pairs,"BasisMatrix"->basis,
    "BasisGram"->gram,"BasisOrthogonalityError"->err,"NWeights"->nWeights|>
];

(* Ambient operator block: Eqs. (G.78)-(G.79); transition analogue in Eqs. (G.83)-(G.84). *)
M3AmbientOperatorBlock[A_?MatrixQ,outPairs_List,inPairs_List] := Module[
  {An,cofA,p,q,kIndex,lIndex,kInputs,lInputs,kCols,lCols,ko,lo,ki,li,iko,ilo,kcol,lcol},
  If[outPairs=={}||inPairs=={},Return[ConstantArray[0.+0.I,{Length[outPairs],Length[inPairs]}]]];
  An=N[A]; cofA=M3CofactorMatrix3[An]; If[FailureQ[cofA],Return[cofA]];
  p=Total[inPairs[[1]]["K"]]; q=Total[inPairs[[1]]["L"]];
  kIndex=M3CompositionIndex3[p]; lIndex=M3CompositionIndex3[q];
  kInputs=DeleteDuplicates[Lookup[inPairs,"K"]]; lInputs=DeleteDuplicates[Lookup[inPairs,"L"]];

  (* LOCAL cache: it exists only during this evaluation. This is useful because it avoids recomputing the same column. *)
  kCols=AssociationThread[M3WtKey/@kInputs,M3SymmetricPowerColumn[An,#]&/@kInputs];
  lCols=AssociationThread[M3WtKey/@lInputs,M3SymmetricPowerColumn[cofA,#]&/@lInputs];

  Table[
    ko=outPairs[[a]]["K"]; lo=outPairs[[a]]["L"]; ki=inPairs[[b]]["K"]; li=inPairs[[b]]["L"];
    iko=Lookup[kIndex,M3WtKey[ko],0]; ilo=Lookup[lIndex,M3WtKey[lo],0];
    kcol=Lookup[kCols,M3WtKey[ki],Missing["NotFound"]]; lcol=Lookup[lCols,M3WtKey[li],Missing["NotFound"]];
    If[iko==0||ilo==0||MissingQ[kcol]||MissingQ[lcol],0.+0.I,kcol[[iko]]*lcol[[ilo]]],
    {a,Length[outPairs]},{b,Length[inPairs]}]
];

(* Physical matrix sigma_W(A): Eqs. (G.79)-(G.82). *)
M3SigmaFromMatrix[finite_Association,iso_Association,A_?MatrixQ,tol_:10^-10] := Module[
  {bd,basis,pairs,Aphys,ambient,raw,rawScale,hermError,hermErrorRel,herm,tr,sigma},

  bd=M3IsospinBasisData[finite,iso,tol];
  If[FailureQ[bd],Return[bd]];

  basis=bd["BasisMatrix"];
  pairs=bd["OccupationPairs"];

  (* rho and Gamma are Hermitian operators. Before inducing the representation, we remove only anti-Hermitian numerical noise from the 3x3 matrix. *)
  Aphys=Chop[(N[A]+ConjugateTranspose[N[A]])/2,tol];

  ambient=M3AmbientOperatorBlock[Aphys,pairs,pairs];
  If[FailureQ[ambient],Return[ambient]];

  raw=ConjugateTranspose[basis] . ambient . basis;
  rawScale=Max[10^-300,Max[Abs[Flatten[N[raw]]]]];
  hermError=Max[Abs[Flatten[N[raw-ConjugateTranspose[raw]]]]];
  hermErrorRel=N[hermError/rawScale];

  (* The test must be RELATIVE. In D^(Lambda)(rho), entries can be physically very small (for example 10^-15) even when the normalized density is perfectly regular. *)
  If[hermErrorRel>10^-6,
    Return[Failure["NonHermitianPhysicalBlock",<|
      "Message"->"The physical block has an excessively large relative Hermiticity error.",
      "HermiticityError"->hermError,
      "RelativeHermiticityError"->hermErrorRel,
      "RawScale"->rawScale|>]]
  ];

  herm=Chop[(raw+ConjugateTranspose[raw])/2,Max[10^-300,tol rawScale]];
  tr=Re[Tr[N[herm]]];

  (* The trace is not compared with tol in absolute terms: for rho it can be extremely small for physical reasons. Only a nonnumeric, nonpositive trace or one below machine range is rejected. *)
  If[!NumericQ[tr]||tr<=0||Abs[tr]<10^-300,
    Return[Failure["ZeroTrace",<|
      "Message"->"The projected block has zero/nonpositive trace or lies outside the numerical range.",
      "Trace"->tr,"RawScale"->rawScale|>]]
  ];

  sigma=Chop[N[herm/tr],10^-12];

  <|"Sigma"->sigma,
    "UnnormalizedBlock"->herm,
    "TraceBeforeNormalization"->tr,
    "HermiticityError"->hermError,
    "RelativeHermiticityError"->hermErrorRel,
    "RawScale"->rawScale,
    "BasisData"->bd|>
]

(* sigma diagnostics *)
M3SigmaDiagnostics[sigma_?MatrixQ,tol_:10^-10] := Module[{herm,ev,probs,purity},
  herm=(sigma+ConjugateTranspose[sigma])/2; ev=Sort[Re[Eigenvalues[N[herm]]],Greater];
  probs=Chop[Re[Diagonal[N[herm]]],tol]; purity=Re[Tr[N[herm . herm]]];
  <|"Eigenvalues"->ev,"PhysicalJProbabilities"->probs,"ProbabilitySum"->Total[probs],"Purity"->purity,
    "TraceError"->Abs[Tr[N[herm]]-1],"MinimumEigenvalue"->Min[ev]|>
];

(* Exact orbital coordinates in the same basis |W,j>: Eqs. (E.73)-(E.81), using the ambient realization (D.31)-(D.36). *)
M3ExactOrbitCoordinates[finite_Association,iso_Association,V_?MatrixQ,tol_:10^-10] := Module[
  {bd,p,q,pairs,basis,highestWeightPair,block,vec,rawCoords,norm,coords},
  If[Dimensions[V]=!={3,3},Return[Failure["InvalidFrame",<|"Message"->"V must be a 3x3 matrix."|>]]];
  bd=M3IsospinBasisData[finite,iso,tol]; If[FailureQ[bd],Return[bd]];
  p=bd["p"]; q=bd["q"]; pairs=bd["OccupationPairs"]; basis=bd["BasisMatrix"];

  (* Eq. (D.31): highest-weight monomial |p,0,0> tensor |0,0,q>. *)
  highestWeightPair={<|"K"->{p,0,0},"L"->{0,0,q}|>};
  block=M3AmbientOperatorBlock[V,pairs,highestWeightPair]; If[FailureQ[block],Return[block]];
  vec=Flatten[block]; rawCoords=ConjugateTranspose[basis] . vec; norm=Norm[rawCoords];
  If[!NumericQ[norm]||norm<=tol,Return[Failure["ZeroProjection",<|"Message"->"The projection of the orbital state has zero norm."|>]]];
  coords=Chop[N[rawCoords/norm],tol];
  <|"Coordinates"->coords,"ProjectionNormSquared"->N[norm^2],"OrbitalProbabilities"->Chop[Abs[coords]^2,tol],
    "JValues"->bd["JValues"],"BasisData"->bd|>
];


(* Independent Wigner-D vs. ambient-construction comparison: Eqs. (E.73)-(E.80). *)
(* APPENDIX G: TRANSITIONS AND DETECTION PROBABILITIES; Eqs. (G.39)-(G.50) and (G.83)-(G.87). *)
(* Finite data for another weight of the same irrep *)
M3FiniteDataAtWeight[finite_Association,nOut_List] := Module[{out,Nrep,r,W},
  Nrep=finite["N"]; r=finite["r"];
  If[Length[nOut]=!=3||!VectorQ[nOut,IntegerQ]||Min[nOut]<0||Total[nOut]=!=Nrep,
    Return[Failure["InvalidOutputWeight",<|"Message"->"nOut must be a nonnegative integer occupation summing to N."|>]]];
  W=nOut-r{1,1,1}; out=Join[finite,<|"nN"->nOut,"W"->W|>]; out
];

(* All output weights present in the irrep *)
M3AllIrrepWeights[finite_Association] := Module[{p,q,r,ks,ls,candidates,valid},
  p=finite["p"]; q=finite["q"]; r=finite["r"]; ks=M3Compositions3[p]; ls=M3Compositions3[q];
  candidates=DeleteDuplicates[Flatten[Table[k+q{1,1,1}-l+r{1,1,1},{k,ks},{l,ls}],1]];
  valid=Select[candidates,Function[nOut,Module[{fd,iso},
    fd=M3FiniteDataAtWeight[finite,nOut]; If[FailureQ[fd],Return[False]];
    iso=M3IsospinData[fd]; AssociationQ[iso]&&iso["Multiplicity"]>0]]];
  SortBy[valid,{#[[3]],#[[2]],#[[1]]}&]
];

(* Interferometer unitarity test *)
M3UnitaryQ3[U_?MatrixQ,tol_:10^-10] := Module[{err},
  If[Dimensions[U]=!={3,3},Return[False]];
  err=Max[Abs[Flatten[N[ConjugateTranspose[U] . U-IdentityMatrix[3]]]]]; err<=tol
];

(* Ambient context for an input weight *)
M3AmbientOperatorContext[A_?MatrixQ,inPairs_List] := Module[{An,cofA,p,q,kIndex,lIndex,kInputs,lInputs,kCols,lCols},
  If[inPairs==={},Return[Failure["EmptyInputWeight",<|"Message"->"The input weight contains no ambient pairs."|>]]];
  An=N[A]; cofA=M3CofactorMatrix3[An]; If[FailureQ[cofA],Return[cofA]];
  p=Total[inPairs[[1]]["K"]]; q=Total[inPairs[[1]]["L"]];
  kIndex=M3CompositionIndex3[p]; lIndex=M3CompositionIndex3[q];
  kInputs=DeleteDuplicates[Lookup[inPairs,"K"]]; lInputs=DeleteDuplicates[Lookup[inPairs,"L"]];
  kCols=AssociationThread[M3WtKey/@kInputs,M3SymmetricPowerColumn[An,#]&/@kInputs];
  lCols=AssociationThread[M3WtKey/@lInputs,M3SymmetricPowerColumn[cofA,#]&/@lInputs];
  If[AnyTrue[Values[kCols],FailureQ]||AnyTrue[Values[lCols],FailureQ],
    Return[Failure["SymmetricPowerFailure",<|"Message"->"Failed to construct the columns of the ambient operator."|>]]];
  <|"A"->An,"p"->p,"q"->q,"InputPairs"->inPairs,"KIndex"->kIndex,"LIndex"->lIndex,"KColumns"->kCols,"LColumns"->lCols|>
];

(* Rectangular block from an already constructed context *)
M3AmbientOperatorBlockFromContext[context_Association,outPairs_List] := Module[
  {inPairs,kIndex,lIndex,kCols,lCols,ko,lo,ki,li,iko,ilo,kcol,lcol},
  inPairs=context["InputPairs"]; kIndex=context["KIndex"]; lIndex=context["LIndex"];
  kCols=context["KColumns"]; lCols=context["LColumns"];
  If[outPairs==={},Return[ConstantArray[0.+0.I,{0,Length[inPairs]}]]];
  Table[
    ko=outPairs[[a]]["K"]; lo=outPairs[[a]]["L"]; ki=inPairs[[b]]["K"]; li=inPairs[[b]]["L"];
    iko=Lookup[kIndex,M3WtKey[ko],0]; ilo=Lookup[lIndex,M3WtKey[lo],0];
    kcol=Lookup[kCols,M3WtKey[ki],Missing["NotFound"]]; lcol=Lookup[lCols,M3WtKey[li],Missing["NotFound"]];
    If[iko==0||ilo==0||MissingQ[kcol]||MissingQ[lcol],0.+0.I,kcol[[iko]]*lcol[[ilo]]],
    {a,Length[outPairs]},{b,Length[inPairs]}]
];

(* TRANSITION MATRICES *)
M3TransitionContext[finite_Association,isoIn_Association,U_?MatrixQ,tol_:10^-10] := Module[{bIn,ambientContext},
  If[!M3UnitaryQ3[U,100 tol],Return[Failure["NonUnitaryInterferometer",<|"Message"->"The interferometer U is not unitary within tolerance."|>]]];
  bIn=M3IsospinBasisData[finite,isoIn,tol]; If[FailureQ[bIn],Return[bIn]];
  ambientContext=M3AmbientOperatorContext[U,bIn["OccupationPairs"]]; If[FailureQ[ambientContext],Return[ambientContext]];
  <|"FiniteData"->finite,"InputIsospinData"->isoIn,"InputBasis"->bIn,"AmbientContext"->ambientContext,"U"->U|>
];

(* Transition to an output weight: Eqs. (G.40)-(G.44) and (G.83)-(G.84). *)
M3TransitionToWeight[context_Association,nOut_List,tol_:10^-10] := Module[{finite,bIn,outFinite,isoOut,bOut,block,T},
  finite=context["FiniteData"]; bIn=context["InputBasis"];
  outFinite=M3FiniteDataAtWeight[finite,nOut]; If[FailureQ[outFinite],Return[outFinite]];
  isoOut=M3IsospinData[outFinite];
  If[FailureQ[isoOut]||isoOut["Multiplicity"]==0,
    Return[Failure["ZeroOutputMultiplicity",<|"Message"->"The output weight does not appear in the irrep."|>]]];
  bOut=M3IsospinBasisData[outFinite,isoOut,tol]; If[FailureQ[bOut],Return[bOut]];
  block=M3AmbientOperatorBlockFromContext[context["AmbientContext"],bOut["OccupationPairs"]]; If[FailureQ[block],Return[block]];
  T=ConjugateTranspose[bOut["BasisMatrix"]] . block . bIn["BasisMatrix"];
  <|"nOut"->nOut,"WOut"->outFinite["W"],"TransitionMatrix"->Chop[N[T],tol],"InputBasis"->bIn,
    "OutputBasis"->bOut,"OutputFiniteData"->outFinite,"OutputIsospinData"->isoOut|>
];


(* Individual transition: simple interface *)
(* Construction of ALL transition matrices *)
M3AllTransitions[finite_Association,isoIn_Association,U_?MatrixQ,tol_:10^-10] := Module[{context,outWeights,transitions},
  context=M3TransitionContext[finite,isoIn,U,tol]; If[FailureQ[context],Return[context]];
  outWeights=M3AllIrrepWeights[finite];
  transitions=Table[M3TransitionToWeight[context,nOut,tol],{nOut,outWeights}];
  If[AnyTrue[transitions,FailureQ],
    Return[Failure["TransitionConstructionFailure",<|"Message"->"At least one output weight could not construct its transition matrix.",
      "Failures"->Select[transitions,FailureQ]|>]]];
  transitions
];

(* Physical probabilities: Eqs. (G.43)-(G.45) and (G.85). *)
M3PhysicalDetectionProbability[T_?MatrixQ,sigma_?MatrixQ,tol_:10^-10] := Module[{z,p},
  If[Dimensions[sigma][[1]]=!=Dimensions[T][[2]]||Dimensions[sigma][[2]]=!=Dimensions[T][[2]],
    Return[Failure["DetectionDimensionMismatch",<|"Message"->"The dimensions of sigma and T are not compatible."|>]]];
  z=Tr[N[T . sigma . ConjugateTranspose[T]]];
  If[Abs[Im[z]]>100 tol Max[1.,Abs[Re[z]]],
    Return[Failure["ComplexProbability",<|"Message"->"The probability has a significant imaginary part.","Value"->z|>]]];
  p=Re[z];
  If[p< -100 tol,Return[Failure["NegativeProbability",<|"Message"->"A significantly negative physical probability was obtained.","Value"->p|>]]];
  If[Abs[p]<tol,0.,N[p]]
];

(* Orbital probability for a transition matrix T: Eqs. (G.46)-(G.48) and (G.86). *)
M3OrbitalDetectionProbability[T_?MatrixQ,c_List] := Module[{v},
  If[Length[c]=!=Dimensions[T][[2]],
    Return[Failure["DetectionDimensionMismatch",<|"Message"->"The dimensions of c and T are not compatible."|>]]];
  v=T . c; N[Re[Conjugate[v] . v]]
];

(* Physical distribution over all output weights: Eqs. (G.45), (G.85), and normalization (G.87). *)
M3PhysicalOutputDistribution[transitions_List,sigma_?MatrixQ,tol_:10^-9] := Module[{results,total},
  results=Table[With[{prob=M3PhysicalDetectionProbability[tr["TransitionMatrix"],sigma,tol]},
    If[FailureQ[prob],prob,<|"nOut"->tr["nOut"],"WOut"->tr["WOut"],
      "Multiplicity"->tr["OutputIsospinData"]["Multiplicity"],"Probability"->prob|>]],{tr,transitions}];
  If[AnyTrue[results,FailureQ],
    Return[Failure["PhysicalProbabilityFailure",<|"Message"->"Failed to evaluate at least one physical probability.","Failures"->Select[results,FailureQ]|>]]];
  total=Total[Lookup[results,"Probability"]];
  <|"Results"->results,"TotalProbability"->N[total],"NormalizationError"->N[Abs[total-1]]|>
];

(* Orbital distribution over all output weights: Eqs. (G.48), (G.86), and normalization (G.87). *)
M3OrbitalOutputDistribution[transitions_List,c_List] := Module[{results,total},
  results=Table[With[{prob=M3OrbitalDetectionProbability[tr["TransitionMatrix"],c]},
    If[FailureQ[prob],prob,<|"nOut"->tr["nOut"],"WOut"->tr["WOut"],
      "Multiplicity"->tr["OutputIsospinData"]["Multiplicity"],"Probability"->prob|>]],{tr,transitions}];
  If[AnyTrue[results,FailureQ],
    Return[Failure["OrbitalProbabilityFailure",<|"Message"->"Failed to evaluate at least one orbital probability.","Failures"->Select[results,FailureQ]|>]]];
  total=Total[Lookup[results,"Probability"]];
  <|"Results"->results,"TotalProbability"->N[total],"NormalizationError"->N[Abs[total-1]]|>
];

(* Physical vs. orbital comparison *)


(* REFERENCE INTERFEROMETERS *)
M3Fourier3[] := Module[{w},
  w=Exp[2 Pi I/3]; 1/Sqrt[3]*{{1,1,1},{1,w,w^2},{1,w^2,w}}
];

M3Interferometer[name_,theta_:Pi/4,phase_:0] := Switch[
  name,
  "Identity",IdentityMatrix[3],
  "Fourier F3",M3Fourier3[],
  "Beam splitter 1-2",{{Cos[theta],-Exp[I phase] Sin[theta],0},{Exp[-I phase] Sin[theta],Cos[theta],0},{0,0,1}},
  _,IdentityMatrix[3]
];


(* APPENDIX F.5: ADJOINT (1,1) EXAMPLE AND FOURIER FRAMES; Eqs. (F.61)-(F.74). *)
(* GT, ISOSPIN, AND GLOBAL SCHUR-WEYL EXTENSIONS *)
(* APPENDIX E.1: INDEPENDENT GELFAND-TSETLIN CONSTRUCTION; Eqs. (E.2)-(E.29). *)
M3GTPatternsU3[lamN_List,nN_List] := Module[{l1,l2,l3,Nrep,nu1,muSum,mu1,mu2,weight,patterns},
  If[Length[lamN]=!=3||Length[nN]=!=3||Total[lamN]=!=Total[nN],Return[Failure["InvalidGTData",<|"Message"->"Lambda and n must have three components and the same sum."|>]]];
  If[!VectorQ[lamN,IntegerQ]||!VectorQ[nN,IntegerQ]||Min[nN]<0,Return[Failure["InvalidGTData",<|"Message"->"Lambda and n must be valid integer vectors."|>]]];
  {l1,l2,l3}=lamN;
  If[!(l1>=l2>=l3>=0),Return[Failure["InvalidPartition",<|"Message"->"Lambda must satisfy Lambda1 >= Lambda2 >= Lambda3 >= 0."|>]]];
  Nrep=Total[lamN];

  (* For a fixed weight: nu1 = n1 mu1+mu2 = n1+n2. *)
  nu1=nN[[1]]; muSum=nN[[1]]+nN[[2]];
  patterns=DeleteCases[Table[
    mu2=muSum-mu1; weight={nu1,mu1+mu2-nu1,Nrep-mu1-mu2};
    If[l1>=mu1>=l2&&l2>=mu2>=l3&&mu1>=nu1>=mu2&&weight===nN,
      <|"Mu1"->mu1,"Mu2"->mu2,"Nu1"->nu1,"Weight"->weight,"j"->(mu1-mu2)/2|>,Nothing],
    {mu1,l2,l1}],Nothing];
  SortBy[patterns,#["j"]&]
];

M3GTMultiplicity[lamN_List,nN_List] := Module[{patterns},
  patterns=M3GTPatternsU3[lamN,nN]; If[FailureQ[patterns],patterns,Length[patterns]]
];


M3GTData[finite_Association] := Module[{patterns,jvals},
  patterns=M3GTPatternsU3[finite["LambdaN"],finite["nN"]]; If[FailureQ[patterns],Return[patterns]];
  jvals=Lookup[patterns,"j",{}]; <|"Patterns"->patterns,"JValues"->jvals,"Multiplicity"->Length[patterns]|>
];

(* GT / ISOSPIN / N(W,j) DIAGNOSTICS *)
M3MultiplicityDiagnostics[finite_Association,iso_Association] := Module[{gt,p,q,W,jIso,jGT,nErrors},
  gt=M3GTData[finite]; If[FailureQ[gt],Return[gt]];
  p=finite["p"]; q=finite["q"]; W=finite["W"]; jIso=iso["JValues"]; jGT=gt["JValues"];
  nErrors=If[jIso==={},{0},Abs[M3NWeightSum[p,q,W,#]-M3NWeight[p,q,W,#]]&/@jIso];
  <|"MultiplicityGT"->gt["Multiplicity"],"MultiplicityIsospin"->iso["Multiplicity"],"GTJValues"->jGT,
    "IsospinJValues"->jIso,"JValuesAgreement"->(Sort[jGT]===Sort[jIso]),
    "MultiplicityAgreement"->(gt["Multiplicity"]==iso["Multiplicity"]),"NClosedSumMaxError"->N[Max[nErrors]]|>
];
(* APPENDICES D/H: PROJECTED NORM IN SU(3); Eqs. (D.37)-(D.52) and (H.3)-(H.14). *)


(* APPENDICES C/D: GENERALIZATION BY EXTERIOR POWERS TO U(M) AND SPECIALIZATION TO SU(3). *)
(* ADDITIONAL OBSERVABLES OF THE PHYSICAL DENSITY *)
M3VonNeumannEntropy[sigma_?MatrixQ,tol_:10^-12] := Module[{herm,ev},
  herm=(sigma+ConjugateTranspose[sigma])/2; ev=Select[Re[Eigenvalues[N[herm]]],#>tol&]; N[-Total[ev Log[ev]]]
];


(* AUXILIARY SCHUR-WEYL DATA *)
M3FiniteDataFromIntegerData[lamN_List,nN_List] := Module[{Nrep},
  If[Length[lamN]=!=3||Length[nN]=!=3||!VectorQ[lamN,IntegerQ]||!VectorQ[nN,IntegerQ]||
    Total[lamN]=!=Total[nN]||Total[lamN]<=0,
    Return[Failure["InvalidIntegerFiniteData",<|"Message"->"Lambda and n must be compatible integer vectors."|>]]];
  Nrep=Total[lamN]; M3FiniteRepresentationData[lamN/Nrep,nN/Nrep,Nrep]
];

(* Dimension d_[Lambda] of the S_N module via hook length. *)
M3SymmetricGroupDimension[lamN_List] := Module[{Nrep,hooks,leg},
  If[!VectorQ[lamN,IntegerQ]||Min[lamN]<0||!OrderedQ[Reverse[lamN]],
    Return[Failure["InvalidPartition",<|"Message"->"Lambda must be a nonincreasing integer partition."|>]]];
  Nrep=Total[lamN];
  hooks=Flatten[Table[leg=Count[Drop[lamN,i],x_/;x>=j];lamN[[i]]-j+leg+1,{i,Length[lamN]},{j,1,lamN[[i]]}]];
  If[hooks==={},1,Factorial[Nrep]/Times@@hooks]
];

M3Partitions3[N_Integer?NonNegative] := PadRight[#,3]&/@IntegerPartitions[N,{1,3}];

(* EXACT PROBABILITY OF A SCHUR-WEYL SECTOR: Eq. (A.33), equivalently Eq. (G.49). *)
M3SectorBlockData[lamN_List,nN_List,Gamma_?MatrixQ,tol_:10^-10] := Module[
  {multiplicity,finite,iso,r,dS,Cn,gammaHerm,detGamma,detFactor,sigmaData,traceReduced,prob},
  multiplicity=M3GTMultiplicity[lamN,nN]; If[FailureQ[multiplicity],Return[multiplicity]];
  If[multiplicity==0,Return[<|"Lambda"->lamN,"Multiplicity"->0,"Probability"->0.|>]];
  finite=M3FiniteDataFromIntegerData[lamN,nN]; If[FailureQ[finite],Return[finite]];
  iso=M3IsospinData[finite]; If[FailureQ[iso],Return[iso]];
  r=finite["r"]; dS=M3SymmetricGroupDimension[lamN]; If[FailureQ[dS],Return[dS]];
  Cn=Multinomial@@nN; gammaHerm=(Gamma+ConjugateTranspose[Gamma])/2; detGamma=Re[Det[N[gammaHerm]]];
  If[detGamma< -100 tol,
    Return[Failure["InvalidGammaDeterminant",<|"Message"->"Gamma has a significantly negative determinant.","Determinant"->detGamma|>]]];
  If[Abs[detGamma]<tol,detGamma=0.]; detFactor=If[r==0,1.,detGamma^r];

  (* If the determinant factor makes the sector exactly zero, sigma does not need to be constructed. *)
  If[detFactor==0.,
    Return[<|"Lambda"->lamN,"Multiplicity"->multiplicity,"FiniteData"->finite,"IsospinData"->iso,
      "DeterminantFactor"->detFactor,"Probability"->0.|>]];

  sigmaData=M3SigmaFromMatrix[finite,iso,Gamma,tol]; If[FailureQ[sigmaData],Return[sigmaData]];
  traceReduced=sigmaData["TraceBeforeNormalization"]; prob=N[(dS/Cn)*detFactor*traceReduced];
  If[prob< -100 tol,
    Return[Failure["NegativeSectorProbability",<|"Message"->"A significantly negative sector probability was obtained.","Probability"->prob|>]]];
  If[Abs[prob]<tol,prob=0.];
  <|"Lambda"->lamN,"Multiplicity"->multiplicity,"FiniteData"->finite,"IsospinData"->iso,
    "SymmetricGroupDimension"->dS,"MultinomialCn"->Cn,"DeterminantGamma"->detGamma,
    "DeterminantFactor"->detFactor,"TraceReduced"->traceReduced,"SigmaData"->sigmaData,"Probability"->prob|>
];


(* DISTRIBUTION OVER ALL IRREPS *)
(* GLOBAL SCHUR-WEYL DETECTION *)
M3FullSchurWeylDetection[nN_List,Gamma_?MatrixQ,U_?MatrixQ,tol_:10^-9] := Module[
  {parts,sectorRows,aggregate,sector,weight,finite,iso,sigma,transitions,detection,key,old,nOut,pOut,results,total},
  parts=Select[M3Partitions3[Total[nN]],TrueQ[M3GTMultiplicity[#,nN]>0]&];
  sectorRows={}; aggregate=<||>;

  Do[
    sector=M3SectorBlockDataCached[lamN,nN,Gamma,tol]; If[FailureQ[sector],Return[sector]];
    weight=sector["Probability"];
    AppendTo[sectorRows,<|"Lambda"->lamN,"Multiplicity"->sector["Multiplicity"],"Probability"->weight|>];

    If[NumericQ[weight]&&weight>0,
      finite=sector["FiniteData"]; iso=sector["IsospinData"]; sigma=sector["SigmaData"]["Sigma"];
      transitions=M3AllTransitions[finite,iso,U,tol]; If[FailureQ[transitions],Return[transitions]];
      detection=M3PhysicalOutputDistribution[transitions,sigma,tol]; If[FailureQ[detection],Return[detection]];

      Do[
        nOut=row["nOut"]; pOut=row["Probability"]; key=M3WtKey[nOut];
        old=Lookup[aggregate,key,<|"nOut"->nOut,"Probability"->0.|>];
        aggregate=Join[aggregate,<|key-><|"nOut"->nOut,"Probability"->N[old["Probability"]+weight pOut]|>|>],
        {row,detection["Results"]}]
    ],
    {lamN,parts}];

  results=SortBy[Values[aggregate],-#["Probability"]&]; total=Total[Lookup[results,"Probability",{}]];
  <|"SectorData"->sectorRows,"Results"->results,"TotalProbability"->N[total],"NormalizationError"->N[Abs[total-1]]|>
];

(* SCHUR-WEYL BLOCK CACHE: avoids reconstructing sigma for the same sector multiple times. *)
M3SectorCacheKey[lamN_List,nN_List,Gamma_?MatrixQ,tol_] :=
  ToString[InputForm[{lamN,nN,Chop[N[Gamma],10^-11],tol}]];

M3SectorBlockDataCached[lamN_List,nN_List,Gamma_?MatrixQ,tol_:10^-9] := Module[{key,val},
  key=M3SectorCacheKey[lamN,nN,Gamma,tol];
  val=Lookup[$M3SectorBlockCache,key,Missing["NotCached"]];
  If[AssociationQ[val]||FailureQ[val],Return[val]];
  val=M3SectorBlockData[lamN,nN,Gamma,tol];
  $M3SectorBlockCache=Join[$M3SectorBlockCache,<|key->val|>];
  val
];

(* BOUNDARIES AND DEGENERATE SPECTRA IN THE FIBER *)
(* Fiber representative for lambda2 = lambda3 *)
M3DegenerateRepresentativeFrame[lam_List,n_List,tol_:10^-10] := Module[
  {d1,d2,q1,c1,complement,c2,c3,V,detV,phase,Q,residual},
  If[Length[lam]=!=3||Length[n]=!=3,Return[Failure["InvalidDegenerateData",<|"Message"->"lambda and n must have three components."|>]]];
  d1=lam[[1]]-lam[[3]]; d2=lam[[2]]-lam[[3]];

  (* This routine only handles the case lambda2 = lambda3. *)
  If[Abs[N[d2]]>tol,Return[Failure["NonDegenerateSpectrum",<|"Message"->"This routine is used only when lambda2 = lambda3."|>]]];

  (* Case 1: fully degenerate spectrum *)
  If[Abs[N[d1]]<=tol,
    If[Max[Abs[N[n-lam]]]>100 tol,Return[Failure["IncompatibleFullyDegenerateSpectrum",<|"Message"->"For a fully degenerate spectrum, n=lambda must hold."|>]]];
    V=IdentityMatrix[3]; Q=Abs[V]^2;
    Return[<|"Point"->Missing["DegenerateChart"],"Branch"->"degenerate","Q"->Q,"V"->V,"DegenerateSpectrum"->True|>]
  ];

  (* Case 2: lambda1 > lambda2 = lambda3 *)
  q1=N[(n-lam[[3]]{1,1,1})/d1];

  (* First-column compatibility. *)
  If[Min[q1]<-tol||Max[q1]>1+tol||Abs[Total[q1]-1]>100 tol,
    Return[Failure["IncompatibleDegenerateFiber",<|"Message"->"The data do not produce a compatible unit-norm first column.","FirstColumnProbabilities"->q1|>]]];

  (* We only remove small numerical boundary errors. *)
  q1=Clip[q1,{0,1}];

  (* Fix the gauge of the first column to be real and nonnegative. *)
  c1=Sqrt[q1];

  (* Complete c1 to an orthonormal basis. This NullSpace is evaluated only in the degenerate case, so it is not a bottleneck. *)
  complement=Orthogonalize[NullSpace[{Conjugate[c1]}]];
  If[Length[complement]<2,Return[Failure["ComplementFailure",<|"Message"->"The first column could not be completed to an orthonormal basis."|>]]];

  {c2,c3}=complement[[1;;2]]; V=Transpose[{c1,c2,c3}];

  (* Fix only one column phase to obtain det(V)=1. The modulus of no entry is changed. *)
  detV=Det[N[V]];
  If[Abs[detV]<=tol,Return[Failure["SingularCompletion",<|"Message"->"The orthonormal completion produced a singular matrix."|>]]];

  phase=detV/Abs[detV]; V[[All,3]]=Conjugate[phase]*V[[All,3]];
  V=Chop[N[V],tol]; Q=Abs[V]^2;

  (* Final test of the physical condition. *)
  residual=Q . N[lam]-N[n];
  If[Norm[residual]>100 tol,Return[Failure["DegenerateFrameResidual",<|"Message"->"The constructed frame does not satisfy Q lambda = n to the required precision.","Residual"->residual|>]]];

  <|"Point"->Missing["DegenerateChart"],"Branch"->"degenerate","Q"->Q,"V"->V,"DegenerateSpectrum"->True|>
];

(* SAMPLING OF CANDIDATE FRAMES *)
M3CandidateFrames[lam_List,n_List,maxFrames_Integer?Positive,gridN_Integer?Positive,includeConjugate_:True,tol_:10^-10] := Module[
  {d2,degenerate,allPts,baseCount,pts,records,Q,heron,branches,V},
  d2=lam[[2]]-lam[[3]];

  (* DEGENERATE SPECTRUM *)
  If[Abs[N[d2]]<=tol,
    degenerate=M3DegenerateRepresentativeFrame[lam,n,tol];
    If[FailureQ[degenerate],Return[degenerate]];
    Return[{degenerate}]
  ];

  (* NONDEGENERATE SPECTRUM *)
  allPts=M3AdmissibleGridPoints[lam,n,gridN]; If[allPts==={},Return[{}]];

  (* If both conjugate branches are included, approximately half of the requested number corresponds to distinct points in the (x,y) plane. *)
  baseCount=Ceiling[maxFrames/If[TrueQ[includeConjugate],2,1]];
  pts=M3FarthestPointSample[allPts,baseCount];

  records=Flatten[Table[
    Q=M3QSection[lam,n,p[[1]],p[[2]]];
    If[FailureQ[Q],Nothing,
      heron=M3HeronPolynomialQ[Q];

      (* H_Heron > 0: unistochastic interior. *)
      (* there are two distinct conjugate branches. *)
      (* *)
      (* H_Heron = 0: orthostochastic boundary. *)
      (* the two branches collapse to a single direction. *)

      branches=If[TrueQ[includeConjugate]&&N[heron]>tol,{1,-1},{1}];
      Table[
        V=M3ReconstructUnitary3[Q,branch,tol];
        If[MatrixQ[V],
          <|"Point"->p,"Branch"->If[branch==1,"+","-"],"Q"->Q,"V"->V,"Heron"->heron,"DegenerateSpectrum"->False|>,
          Nothing],
        {branch,branches}]
    ],
    {p,pts}],1];

  Take[records,UpTo[maxFrames]]
];


(* APPENDIX B / FIGS. 2--4: GEOMETRY FIGURES *)
$M3GeometryBlue = RGBColor[.12,.42,.70];
$M3GeometryFill = RGBColor[.78,.87,.94];
$M3GeometryGreen = RGBColor[.05,.50,.30];
$M3GeometryOrange = RGBColor[.90,.45,.08];
$M3GeometryGray = GrayLevel[.55];
$M3GeometryImageSize = {500,380};

M3GeometryFrameLabel[s_] := Style[s,12,FontFamily->"Times",Black];

(* bistochastic region -- unchanged criterion *)
M3FiberBistochasticPointQ[lam_List,n_List,x_?NumericQ,y_?NumericQ,tol_:10^-9] := Module[{Q},
  Q=M3QSection[lam,n,x,y]; If[FailureQ[Q],Return[False]]; Min[Flatten[N[Q]]]>=-tol
];

(* The three expensive background objects are built only once in the Geometry tab
   and reused by the two fiber plots. *)
M3GeometryBasePlots[lam_List,n_List] := Module[{pB,pU,pH,heronAt},
  heronAt[x_?NumericQ,y_?NumericQ] := Module[{Q=M3QSection[lam,n,x,y]},
    If[FailureQ[Q],Indeterminate,N[M3HeronPolynomialQ[Q]]]
  ];

  pB=RegionPlot[
    M3FiberBistochasticPointQ[lam,n,x,y],{x,0,1},{y,0,1},
    PlotStyle->Directive[White,Opacity[0]],
    BoundaryStyle->Directive[$M3GeometryGray,Dashed,AbsoluteThickness[1.35]],
    PlotPoints->42,MaxRecursion->1,PerformanceGoal->"Quality"
  ];

  pU=RegionPlot[
    M3FrameAdmissiblePointQ[lam,n,x,y],{x,0,1},{y,0,1},
    PlotStyle->Directive[$M3GeometryFill,Opacity[.62]],
    BoundaryStyle->None,
    PlotPoints->44,MaxRecursion->1,PerformanceGoal->"Quality"
  ];

  pH=ContourPlot[
    heronAt[x,y],{x,0,1},{y,0,1},
    Contours->{0},ContourShading->False,
    ContourStyle->Directive[$M3GeometryBlue,AbsoluteThickness[2.0]],
    RegionFunction->Function[{xx,yy,z},TrueQ[M3FiberBistochasticPointQ[lam,n,xx,yy]]],
    PlotPoints->56,MaxRecursion->1,PerformanceGoal->"Quality"
  ];

  <|"Bistochastic"->pB,"Unistochastic"->pU,"HeronBoundary"->pH|>
];

M3GeometryRangeFromRecords[records_List] := Module[{pts,xRange,yRange,dx,dy},
  pts=Select[Lookup[records,"Point",{}],ListQ[#]&&Length[#]==2&];
  If[pts==={},Return[{{0,1},{0,1}}]];
  xRange=MinMax[N[pts[[All,1]]]]; yRange=MinMax[N[pts[[All,2]]]];
  dx=Max[.075,.23 Max[10^-6,xRange[[2]]-xRange[[1]]]];
  dy=Max[.075,.23 Max[10^-6,yRange[[2]]-yRange[[1]]]];
  {
    {Max[0.,xRange[[1]]-dx],Min[1.,xRange[[2]]+dx]},
    {Max[0.,yRange[[1]]-dy],Min[1.,yRange[[2]]+dy]}
  }
];

M3BistochasticLegend[] := Framed[
  Grid[{
    {
      Graphics[{Directive[$M3GeometryGray,Dashed,AbsoluteThickness[1.35]],Line[{{0,.5},{1,.5}}]},ImageSize->{30,10},PlotRange->{{0,1},{0,1}}],
      Style["bistochastic boundary",9.4,FontFamily->"Times"]
    },
    {
      Graphics[{Directive[$M3GeometryFill,Opacity[.62]],Rectangle[{0,0},{1,1}]},ImageSize->{30,10},PlotRange->{{0,1},{0,1}}],
      Style["unistochastic region",9.4,FontFamily->"Times"]
    },
    {
      Graphics[{Directive[$M3GeometryBlue,AbsoluteThickness[2.0]],Line[{{0,.5},{1,.5}}]},ImageSize->{30,10},PlotRange->{{0,1},{0,1}}],
      Style[TraditionalForm[Subscript["H","Heron"]==0],9.4,FontFamily->"Times"]
    }
  },Alignment->Left,Spacings->{.55,.32}],
  Background->White,
  FrameStyle->Directive[GrayLevel[.82],AbsoluteThickness[.6]],
  FrameMargins->{{7,7},{5,5}},RoundingRadius->4
];

M3FrameLegend[] := Framed[
  Grid[{
    {
      Graphics[{Directive[$M3GeometryGray,Dashed,AbsoluteThickness[1.35]],Line[{{0,.5},{1,.5}}]},ImageSize->{27,10},PlotRange->{{0,1},{0,1}}],
      Style["bistochastic boundary",9.0,FontFamily->"Times"],
      Graphics[{Directive[$M3GeometryFill,Opacity[.62]],Rectangle[{0,0},{1,1}]},ImageSize->{18,10},PlotRange->{{0,1},{0,1}}],
      Style["unistochastic region",9.0,FontFamily->"Times"]
    },
    {
      Graphics[{$M3GeometryGreen,Disk[{.5,.5},.25]},ImageSize->{18,12},PlotRange->{{0,1},{0,1}}],
      Style["independent frame",9.0,FontFamily->"Times"],
      Graphics[{$M3GeometryOrange,Disk[{.5,.5},.25]},ImageSize->{18,12},PlotRange->{{0,1},{0,1}}],
      Style["redundant frame",9.0,FontFamily->"Times"]
    },
    {
      Graphics[{GrayLevel[.28],Disk[{.5,.5},.25]},ImageSize->{18,12},PlotRange->{{0,1},{0,1}}],
      Style["branch +",9.0,FontFamily->"Times"],
      Graphics[{GrayLevel[.28],RegularPolygon[{.5,.5},{.28,Pi/4},4]},ImageSize->{18,12},PlotRange->{{0,1},{0,1}}],
      Style["branch -",9.0,FontFamily->"Times"]
    }
  },Alignment->Left,Spacings->{.42,.30}],
  Background->White,
  FrameStyle->Directive[GrayLevel[.82],AbsoluteThickness[.6]],
  FrameMargins->{{7,7},{5,5}},RoundingRadius->4
];

(* Clean bistochastic/unistochastic overview added to Geometry. *)
M3BistochasticUnistochasticFigure[lam_List,n_List,records_List,base_:Automatic] := Module[
  {plots,range,legend},
  plots=If[AssociationQ[base],base,M3GeometryBasePlots[lam,n]];
  range=M3GeometryRangeFromRecords[records];
  legend=M3BistochasticLegend[];

  Show[
    plots["Unistochastic"],plots["Bistochastic"],plots["HeronBoundary"],
    Frame->True,Axes->False,AspectRatio->1,
    PlotRange->range,PlotRangePadding->Scaled[.018],
    FrameStyle->Directive[Black,AbsoluteThickness[.9]],
    FrameTicksStyle->Directive[Black,10],
    FrameLabel->{
      M3GeometryFrameLabel[TraditionalForm[Subscript["q",11]]],
      M3GeometryFrameLabel[TraditionalForm[Subscript["q",21]]]
    },
    Background->White,ImageSize->$M3GeometryImageSize,
    ImagePadding->{{52,18},{42,15}},BaseStyle->{FontFamily->"Times",FontSize->10.2},
    Epilog->{Inset[legend,Scaled[{.035,.035}],{Left,Bottom}]}
  ]
];

(* Candidate frames on the fiber: same data and branch information as the original plot. *)
M3FiberFramesFigure[lam_List,n_List,records_List,independentIndices_List:{},base_:Automatic] := Module[
  {validIndices,localRecords,pts,branches,range,xRange,yRange,xSpan,ySpan,branchShift,branchDir,
   dispPts,groups,markerR,centerR,connectors,markers,independentQ,plots,legend},

  validIndices=Select[Range[Length[records]],Function[i,ListQ[records[[i]]["Point"]]&&Length[records[[i]]["Point"]]==2]];
  If[validIndices==={},Return[Panel["The (x,y) chart is not available for these frames."]]];

  localRecords=records[[validIndices]];
  pts=Lookup[localRecords,"Point"];
  branches=Lookup[localRecords,"Branch"];
  independentQ[i_]:=MemberQ[independentIndices,validIndices[[i]]];

  range=M3GeometryRangeFromRecords[localRecords];
  xRange=range[[1]]; yRange=range[[2]];
  xSpan=Max[10^-6,xRange[[2]]-xRange[[1]]];
  ySpan=Max[10^-6,yRange[[2]]-yRange[[1]]];
  branchShift=.018 xSpan;
  branchDir[b_]:=Which[
    StringContainsQ[ToString[b],"+"],-1,
    StringContainsQ[ToString[b],"-"],1,
    True,0
  ];
  dispPts=Table[pts[[i]]+{branchShift branchDir[branches[[i]]],0},{i,Length[pts]}];
  groups=GatherBy[Range[Length[pts]],Round[pts[[#]],10^-9]&];
  markerR=.0105 Min[xSpan,ySpan];
  centerR=.0030 Min[xSpan,ySpan];

  connectors=Table[
    With[{sg=SortBy[g,dispPts[[#,1]]&],qq=pts[[First[g]]]},
      {
        Directive[GrayLevel[.56],AbsoluteThickness[.8]],
        If[Length[sg]>=2,Line[{dispPts[[First[sg]]],dispPts[[Last[sg]]]}],{}],
        GrayLevel[.32],Disk[qq,centerR]
      }
    ],
    {g,groups}
  ];

  markers=Table[
    With[{
      c=If[independentQ[i],$M3GeometryGreen,$M3GeometryOrange],
      pp=dispPts[[i]],
      plusQ=StringContainsQ[ToString[branches[[i]]],"+"]
    },
      {
        EdgeForm[Directive[White,AbsoluteThickness[1.0]]],c,
        If[plusQ,Disk[pp,1.10 markerR],RegularPolygon[pp,{1.15 markerR,Pi/4},4]],
        Text[
          Style[Row[{validIndices[[i]],branches[[i]]}],8.7,Bold,FontFamily->"Times",FontColor->c],
          pp+{0,.016 ySpan}
        ]
      }
    ],
    {i,Length[pts]}
  ];

  plots=If[AssociationQ[base],base,M3GeometryBasePlots[lam,n]];
  legend=M3FrameLegend[];

  Show[
    plots["Unistochastic"],plots["Bistochastic"],plots["HeronBoundary"],
    Graphics[{connectors,markers}],
    Frame->True,Axes->False,PlotRange->range,AspectRatio->1,
    PlotRangePadding->Scaled[.018],
    FrameStyle->Directive[Black,AbsoluteThickness[.9]],
    FrameTicksStyle->Directive[Black,10],
    FrameLabel->{
      M3GeometryFrameLabel[TraditionalForm[Subscript["q",11]]],
      M3GeometryFrameLabel[TraditionalForm[Subscript["q",21]]]
    },
    Background->White,ImageSize->$M3GeometryImageSize,
    ImagePadding->{{52,18},{42,15}},BaseStyle->{FontFamily->"Times",FontSize->10.2},
    Epilog->{Inset[legend,Scaled[{.035,.035}],{Left,Bottom}]}
  ]
];

(* Improved unitarity triangle.  The side lengths and Heron quantity are exactly
   those computed by M3TriangleLengths and M3HeronPolynomialQ. *)
M3UnitarityTriangleFigure[Q_?MatrixQ] := Module[
  {L,l1,l2,l3,xC,yC,pts,H,tri,summary,blue,fill},

  L=M3TriangleLengths[Q];
  {l1,l2,l3}=N[L];
  H=N[M3HeronPolynomialQ[Q]];
  blue=$M3GeometryBlue;
  fill=$M3GeometryFill;

  If[Max[L]<=10^-12,Return[Panel["Completely degenerate triangle."]]];
  If[l1<=10^-12,Return[Panel[Row[{"Degenerate triangle;  H_Heron = ",NumberForm[H,{6,4}]}]]]];

  xC=(l1^2+l3^2-l2^2)/(2 l1);
  yC=Sqrt[Max[0,l3^2-xC^2]];
  pts={{0,0},{l1,0},{xC,yC}};

  tri=Graphics[
    {
      Directive[fill,Opacity[.62]],Polygon[pts],
      Directive[blue,AbsoluteThickness[2.25],JoinForm["Round"]],Line[Append[pts,First[pts]]],
      Directive[GrayLevel[.28],PointSize[.016]],Point[pts],

      Text[
        Style[Row[{TraditionalForm[Subscript["L",1]]," = ",NumberForm[l1,{5,3}]}],10.4,FontFamily->"Times",Black],
        Mean[pts[[{1,2}]]]+{0,-.050}
      ],
      Text[
        Style[Row[{TraditionalForm[Subscript["L",2]]," = ",NumberForm[l2,{5,3}]}],10.4,FontFamily->"Times",Black],
        Mean[pts[[{2,3}]]]+{.040,.018}
      ],
      Text[
        Style[Row[{TraditionalForm[Subscript["L",3]]," = ",NumberForm[l3,{5,3}]}],10.4,FontFamily->"Times",Black],
        Mean[pts[[{3,1}]]]+{-.040,.018}
      ]
    },
    PlotRange->All,PlotRangePadding->Scaled[.16],
    ImagePadding->10,ImageSize->355,Background->White
  ];

  summary=Column[
    {
      Style[TraditionalForm[Subscript["H","Heron"]==16 \[ScriptCapitalA]^2],11.4,FontFamily->"Times"],
      Style[Row[{TraditionalForm[Subscript["H","Heron"]]," = ",NumberForm[H,{8,6}]}],10.8,FontFamily->"Times"]
    },
    Alignment->Center,Spacings->.28
  ];

  Pane[
  Column[
    {tri,summary},
    Alignment->Center,
    Spacings->.20
  ],
  ImageSize->$M3GeometryImageSize,
  Alignment->Center
]
];


(* Compact selected-frame summary for the Geometry tab.  Presentation only:
   all numerical quantities reuse the existing definitions without altering the physics. *)
M3SelectedFrameSummaryPanel[lam_List,n_List,Q_?MatrixQ,V_?MatrixQ] := Module[
  {H,j2,residual,reconstructionError,matrixCard,diagnosticsCard,titleColor,
   matrixDisplay,residualLabel,reconstructionLabel},

  H=N[M3HeronPolynomialQ[Q]];
  j2=H/4;
  residual=Norm[Q . N[lam]-N[n]];
  reconstructionError=Max[Abs[Flatten[Abs[V]^2-Q]]];
  titleColor=RGBColor[.10,.28,.50];

  matrixDisplay=MatrixForm[
    Map[NumberForm[#,{6,4}]&,N[Q],{2}]
  ];

  residualLabel=Style[
    Row[{
      "|| ",Style["Q",Italic]," ",TraditionalForm[Overscript[\[Lambda],"_"]],
      " - ",TraditionalForm[Overscript["n","_"]]," ||"
    }],
    10.2,Bold,FontFamily->"Times"
  ];

  reconstructionLabel=Style[
    Row[{
      "max  || ",Superscript[Row[{"|",Style["V",Italic],"|"}],2],
      " - ",Style["Q",Italic]," ||"
    }],
    10.2,Bold,FontFamily->"Times"
  ];

  matrixCard=Framed[
    Column[{
      Style["Bistochastic matrix  Q(V)",11.2,Bold,FontFamily->"Times",FontColor->titleColor],
      Style[matrixDisplay,11.0,FontFamily->"Times"]
    },Alignment->Center,Spacings->.55],
    Background->White,
    FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.7]],
    FrameMargins->{{18,18},{11,11}},
    RoundingRadius->5
  ];

  diagnosticsCard=Framed[
    Grid[
      {
        {
          Style[TraditionalForm[Subscript["H","Heron"]],10.5,Bold,FontFamily->"Times"],
          Style[NumberForm[H,{9,7}],10.2,FontFamily->"Times"]
        },
        {
          Style[TraditionalForm[Superscript["J",2]],10.5,Bold,FontFamily->"Times"],
          Style[NumberForm[j2,{9,7}],10.2,FontFamily->"Times"]
        },
        {
          residualLabel,
          Style[ScientificForm[residual,3],10.2,FontFamily->"Times"]
        },
        {
          reconstructionLabel,
          Style[ScientificForm[reconstructionError,3],10.2,FontFamily->"Times"]
        }
      },
      Alignment->{{Left,Right}},
      Spacings->{1.3,.72}
    ],
    Background->White,
    FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.7]],
    FrameMargins->{{15,15},{10,10}},
    RoundingRadius->5
  ];

  Framed[
    Column[{
      Style["Selected frame",12.5,Bold,FontFamily->"Times",FontColor->titleColor],
      Grid[
        {{matrixCard,diagnosticsCard}},
        Alignment->{Center,Center},
        Spacings->{1.5,0}
      ]
    },Alignment->Center,Spacings->.70],
    Background->White,
    FrameStyle->Directive[GrayLevel[.88],AbsoluteThickness[.75]],
    FrameMargins->{{16,16},{12,12}},
    RoundingRadius->6
  ]
];

(* Conjugate-branch fiber figure.  This uses only the branch structure V_+ and V_-. *)
M3ConjugateBranchSampleGroups[lam_List,n_List,records_List,maxPoints_Integer:6,tol_:10^-10] := Module[
  {recordGroups,usableRecordGroups,needed,auxGridN,auxPts,auxInteriorPts,existingPts,auxGroups},

  recordGroups=GatherBy[
    Select[records,KeyExistsQ[#,"Point"]&&ListQ[#Point]&&Length[#Point]==2&],
    Round[#Point,10^-10]&
  ];

  usableRecordGroups=Select[
    recordGroups,
    Sort[DeleteDuplicates[Lookup[#,"Branch",""]]]=={"+","-"}&
  ];
  usableRecordGroups=Take[usableRecordGroups,UpTo[maxPoints]];
  needed=maxPoints-Length[usableRecordGroups];
  If[needed<=0,Return[usableRecordGroups]];

  (* If the current frame bank does not contain enough explicit conjugate pairs,
     complete the display sample from the same admissible fiber without changing
     the fiber definition or the unitary reconstruction. *)
  auxGridN=Max[24,4 maxPoints+8];
  auxPts=M3AdmissibleGridPoints[lam,n,auxGridN];
  existingPts=Round[Lookup[First /@ usableRecordGroups,"Point",{}],10^-10];

  auxInteriorPts=Select[
    M3FarthestPointSample[auxPts,4 maxPoints],
    Module[{Q=M3QSection[lam,n,#[[1]],#[[2]]]},
      MatrixQ[Q]&&N[M3HeronPolynomialQ[Q]]>tol&&
      !MemberQ[existingPts,Round[#,10^-10]]
    ]&
  ];
  auxInteriorPts=Take[DeleteDuplicatesBy[auxInteriorPts,Round[#,10^-10]&],UpTo[needed]];

  auxGroups=Map[
    Function[p,
      Module[{Q,heron,Vp,Vm},
        Q=M3QSection[lam,n,p[[1]],p[[2]]];
        heron=M3HeronPolynomialQ[Q];
        Vp=M3ReconstructUnitary3[Q,1,tol];
        Vm=M3ReconstructUnitary3[Q,-1,tol];
        If[MatrixQ[Vp]&&MatrixQ[Vm],
          {
            <|"Point"->p,"Branch"->"+","Q"->Q,"V"->Vp,"Heron"->heron,"DegenerateSpectrum"->False|>,
            <|"Point"->p,"Branch"->"-","Q"->Q,"V"->Vm,"Heron"->heron,"DegenerateSpectrum"->False|>
          },
          Nothing
        ]
      ]
    ],
    auxInteriorPts
  ];

  Take[Join[usableRecordGroups,DeleteCases[auxGroups,Nothing]],UpTo[maxPoints]]
];

M3ConjugateBranchFiberFigure[lam_List,n_List,records_List,base_:Automatic,maxPoints_Integer:6] := Module[
  {plots,range,xRange,yRange,xSpan,ySpan,shift,sampleGroups,samplePts,legend,markers},

  plots=If[AssociationQ[base],base,M3GeometryBasePlots[lam,n]];
  range=M3GeometryRangeFromRecords[records];
  xRange=range[[1]]; yRange=range[[2]];
  xSpan=Max[10^-6,xRange[[2]]-xRange[[1]]];
  ySpan=Max[10^-6,yRange[[2]]-yRange[[1]]];
  shift=.016 xSpan;

  sampleGroups=M3ConjugateBranchSampleGroups[lam,n,records,maxPoints];
  If[sampleGroups==={},Return[Panel["No interior conjugate branch pairs are available for the current fiber."]]];
  samplePts=Lookup[#,"Point"]& /@ (First /@ sampleGroups);

  legend=Framed[
    Grid[{{
      Graphics[{Directive[$M3GeometryBlue],Disk[{.5,.5},.22]},ImageSize->{16,12},PlotRange->{{0,1},{0,1}}],
      Style[TraditionalForm[Subscript["V","+"]],9.2,FontFamily->"Times"],
      Graphics[{Directive[$M3GeometryOrange],RegularPolygon[{.5,.5},{.27,Pi/4},4]},ImageSize->{16,12},PlotRange->{{0,1},{0,1}}],
      Style[TraditionalForm[Subscript["V","-"]],9.2,FontFamily->"Times"]
    }},Alignment->Left,Spacings->{.45,.25}],
    Background->White,
    FrameStyle->Directive[GrayLevel[.82],AbsoluteThickness[.6]],
    FrameMargins->{{6,6},{4,4}},RoundingRadius->4
  ];

  markers=Flatten@Table[
    Module[{pt=samplePts[[i]],pp,pm},
      pp=pt+{-shift,0};
      pm=pt+{shift,0};
      {
        {GrayLevel[.32],Disk[pt,.0028 Min[xSpan,ySpan]]},
        {EdgeForm[Directive[White,AbsoluteThickness[.95]]],$M3GeometryBlue,Disk[pp,.0105 Min[xSpan,ySpan]]},
        Text[Style[Row[{i,"+"}],8.8,Bold,FontFamily->"Times",FontColor->$M3GeometryBlue],pp+{0,.018 ySpan}],
        {EdgeForm[Directive[White,AbsoluteThickness[.95]]],$M3GeometryOrange,RegularPolygon[pm,{.0112 Min[xSpan,ySpan],Pi/4},4]},
        Text[Style[Row[{i,"-"}],8.8,Bold,FontFamily->"Times",FontColor->$M3GeometryOrange],pm+{0,.018 ySpan}]
      }
    ],
    {i,Length[samplePts]}
  ,1];

  Show[
    plots["Unistochastic"],plots["Bistochastic"],plots["HeronBoundary"],
    Graphics[markers],
    Frame->True,Axes->False,AspectRatio->1,
    PlotRange->range,PlotRangePadding->Scaled[.018],
    FrameStyle->Directive[Black,AbsoluteThickness[.9]],
    FrameTicksStyle->Directive[Black,10],
    FrameLabel->{
      M3GeometryFrameLabel[TraditionalForm[Subscript["q",11]]],
      M3GeometryFrameLabel[TraditionalForm[Subscript["q",21]]]
    },
    Background->White,ImageSize->$M3GeometryImageSize,
    ImagePadding->{{52,18},{42,15}},BaseStyle->{FontFamily->"Times",FontSize->10.2},
    Epilog->{Inset[legend,Scaled[{.035,.035}],{Left,Bottom}]}
  ]
];

M3GramHeatmap[G_?MatrixQ] := Module[
  {A,K,cf,cells,barX0,barX1,barY0,barY1,nBar,barRects,tickVals,tickPrimitives},

  A=Reverse[Chop[Abs[N[G]]^2,10^-14]];
  K=Length[A];
  If[K==0,Return[Graphics[{},Background->White,ImageSize->470]]];

  cf=ColorData["LakeColors"];

  cells=Flatten@Table[
    {
      EdgeForm[Directive[White,AbsoluteThickness[.55]]],
      FaceForm[cf[Clip[A[[i,j]],{0,1}]]],
      Rectangle[{j-1,i-1},{j,i}]
    },
    {i,1,K},{j,1,K}
  ,1];

  barX0=K+0.30;
  barX1=K+0.70;
  barY0=0.18 K;
  barY1=0.82 K;
 nBar = 512;

barRects = Raster[
  Table[
    {
      List @@ ColorConvert[cf[v], "RGB"]
    },
    {v, N@Subdivide[0, 1, nBar - 1]}
  ],
  {{barX0, barY0}, {barX1, barY1}}
];

  tickVals=N@Range[0,1,.2];
  tickPrimitives=Flatten@Table[
    {
      Directive[GrayLevel[.25],AbsoluteThickness[.65]],
      Line[{{barX1,barY0+v(barY1-barY0)},{barX1+.10,barY0+v(barY1-barY0)}}],
      Text[
        Style[NumberForm[v,{2,1}],8.5,FontFamily->"Times"],
        {barX1+.15,barY0+v(barY1-barY0)},
        {-1,0}
      ]
    },
    {v,tickVals}
  ,1];

  Graphics[
    {
      cells,
      Directive[Black,AbsoluteThickness[.85]],
      FaceForm[None],EdgeForm[Directive[Black,AbsoluteThickness[.85]]],
      Rectangle[{0,0},{K,K}],

      barRects,
      EdgeForm[Directive[GrayLevel[.35],AbsoluteThickness[.7]]],
      FaceForm[None],
      Rectangle[{barX0,barY0},{barX1,barY1}],
      tickPrimitives,

      Text[M3FrameLabel["frame a"],{K/2,-.20}],
      Text[Rotate[M3FrameLabel["frame b"],Pi/2],{-.20,K/2}],
      Text[
        Style[
 Row[{
   "|",
   Subscript[
    Style["G", Italic],
    Row[{Style["a", Italic], Style["b", Italic]}]
   ],
   "|",
   Superscript["", 2]
 }],
 9.5,
 Bold,
 FontFamily -> "Times"
],
        {(barX0+barX1)/2,barY1+.28}
      ]
    },
    PlotRange->{{-.42,K+1.15},{-.42,K+.18}},
    PlotRangePadding->Scaled[.01],
    AspectRatio->Automatic,
    Background->White,
    ImageSize->{500,380},
    BaseStyle->{FontFamily->"Times",FontSize->10.2}
  ]
];


(* APPENDIX H.6 / FIG. 14 \[LongDash] HESSIAN FIGURES: parameters, input data, and visual style; reuses the stationary point, 
reduced Hessian, its diagonalization, and the exact toroidal action defined elsewhere. *)

$M3HessianExactColor = RGBColor[.12,.34,.62];
$M3HessianQuadraticColor = RGBColor[.86,.43,.10];
$M3HessianDirection1Color = RGBColor[.12,.34,.62];
$M3HessianDirection2Color = RGBColor[.16,.56,.42];
$M3HessianStationaryColor = RGBColor[.72,.15,.10];

(* APPENDIX H.2 \[LongDash] STATIONARY POINT, Eqs. (H.25)\[Dash](H.26): the toroidal action uses phi_1 = 0, with stationary point (phi_2, phi_3) = (0, 0). *)
M3HessianStationaryPoint[] := {0.,0.};
(* APPENDIX H.2: verification of the stationary condition*)
M3HessianStationaryCheck[lam_List,n_List,Q_?MatrixQ,tol_:10^-9] := Module[
  {res},
  res = M3StationaryResidual[lam,n,Q];
  <|
    "Residual" -> res,
    "ResidualNorm" -> Norm[N[res]],
    "SatisfiedQ" -> TrueQ[Norm[N[res]] <= 100 tol]
  |>
];
(* APPENDIX H.4: CONSTRUCTION OF H_red, Eqs. (H.54)-(H.59)                  *)
M3HessianReducedData[lam_List,Q_?MatrixQ,tol_:10^-10] := Module[
  {hData,Hred,evals,evecs,ord},

  hData = M3HessianData[lam,Q];
  If[FailureQ[hData], Return[hData]];

  Hred = N[(hData["ReducedHessian"] + ConjugateTranspose[hData["ReducedHessian"]])/2];
  {evals,evecs} = Eigensystem[Hred];
  ord = Ordering[Re[evals]];
  evals = Chop[Re[evals[[ord]]],10^-13];
  evecs = Normalize /@ Chop[Re[evecs[[ord]]],10^-13];

  If[!VectorQ[evals,NumericQ] || Length[evals] =!= 2,
    Return[Failure["InvalidReducedHessian", <|
      "Message" -> "The reduced Hessian could not be diagonalized numerically."
    |>]]
  ];

If[Min[evals]<=tol,
  Return[Failure["NonPositiveReducedHessian",<|
    "Message"->"The reduced Hessian is not positive and nondegenerate at the selected stationary frame.",
    "Eigenvalues"->evals
  |>]]
];

<|"FullHessianData"->hData,"ReducedHessian"->Hred,"Eigenvalues"->evals,
  "Eigenvectors"->evecs,"ReducedDeterminant"->N[hData["ReducedDeterminant"]]|>
];
(* APPENDIX H.6 / FIG. 14(a) \[LongDash] DIAGONALIZATION OF H_red AND Delta s_exact: diagonalization is 
included in M3HessianReducedData; Delta s_exact follows Eq. (H.75). *)

M3HessianExactDeltaFunction[finite_Association,Q_?MatrixQ,phiStar_List] := Module[
  {sStar},
  sStar = Re[M3TorusActionIntensive[finite,Q,phiStar[[1]],phiStar[[2]]]];
  Function[{x,y},
    Quiet @ Check[
      Re[M3TorusActionIntensive[finite,Q,x,y]] - sStar,
      Indeterminate
    ]
  ]
];

(* APPENDIX H.6 \[LongDash] Delta s_quad, Eq. (H.75): Delta s_quad = -1/2 delta_phi^T H_red delta_phi, using the existing sign convention. *)

M3HessianQuadraticDeltaFunction[Hred_?MatrixQ,phiStar_List] := Function[{x,y},
  -0.5 ({x,y} - phiStar) . Hred . ({x,y} - phiStar)
];

(* APPENDIX H.6 / FIG. 14 \[LongDash] LOCAL WINDOW: centered at the stationary point and scaled from a quadratic reference ellipse so the 
outer exact contour fills about 70\[Dash]80% of the useful plotting area. *)

M3HessianArticleData[lam_List,n_List,finite_Association,Q_?MatrixQ,tol_:10^-9] := Module[
  {phiStar,stationary,hessian,Hred,evals,evecs,exactDelta,quadDelta,probeHalfWidth,
   probeGrid,probeVals,finiteProbeVals,outerDrop,semiaxes,localHalfWidth,contourLevels},

  phiStar=M3HessianStationaryPoint[];
  stationary=M3HessianStationaryCheck[lam,n,Q,tol];
  hessian=M3HessianReducedData[lam,Q,Max[10^-12,tol/100]];
  If[FailureQ[hessian],Return[hessian]];

  Hred=hessian["ReducedHessian"]; evals=hessian["Eigenvalues"]; evecs=hessian["Eigenvectors"];
  exactDelta=M3HessianExactDeltaFunction[finite,Q,phiStar];
  quadDelta=M3HessianQuadraticDeltaFunction[Hred,phiStar];

  probeHalfWidth=Pi/4;
  probeGrid=Subdivide[-probeHalfWidth,probeHalfWidth,20];
  probeVals=Flatten@Table[exactDelta[x,y],{x,probeGrid},{y,probeGrid}];
  finiteProbeVals=Select[probeVals,NumericQ];
  outerDrop=If[finiteProbeVals==={},0.02,
    Clip[0.22 Abs[Min[finiteProbeVals]],{5*10^-4,0.12}]
  ];

  semiaxes=Sqrt[2 outerDrop/evals];
  localHalfWidth=Clip[(4/3) Max[semiaxes],{0.24,Pi/2}];
  contourLevels=Sort[-Rest@Subdivide[0,outerDrop,6]];

  <|"StationaryPoint"->phiStar,"StationaryCheck"->stationary,"ReducedHessian"->Hred,
    "Eigenvalues"->evals,"Eigenvectors"->evecs,"ReducedDeterminant"->hessian["ReducedDeterminant"],
    "ExactDeltaFunction"->exactDelta,"QuadraticDeltaFunction"->quadDelta,"OuterDrop"->outerDrop,
    "EllipseSemiaxes"->semiaxes,"LocalHalfWidth"->localHalfWidth,"ContourLevels"->contourLevels|>
];


M3HessianFrameOptions[imageSize_:430,padding_:{{54,18},{42,16}}] := Sequence[
  Frame->True,Axes->False,AspectRatio->1,
  FrameStyle->Directive[Black,AbsoluteThickness[.9]],
  FrameTicksStyle->Directive[Black,9.5],
  FrameLabel->{M3FrameLabel[TraditionalForm[Subscript[\[Phi],2]]],
    M3FrameLabel[TraditionalForm[Subscript[\[Phi],3]]]},
  Background->White,ImageSize->imageSize,ImagePadding->padding,BaseStyle->$M3BaseStyle
];

M3HessianStationaryMarker[phiStar_List] := Graphics[{
  Directive[White,PointSize[.024]],Point[phiStar],
  Directive[$M3HessianStationaryColor,PointSize[.013]],Point[phiStar]
}];

M3HessianPrincipalDirectionOverlay[phiStar_List,evecs_List,scales_List,evals_List] := Module[
  {infoBox,vminText,vmaxText,lamMinText,lamMaxText},
  vminText=Row[{"(",NumberForm[N[evecs[[1,1]]],{5,4},NumberPadding->{"","0"}],", ",
    NumberForm[N[evecs[[1,2]]],{5,4},NumberPadding->{"","0"}],")"}];
  vmaxText=Row[{"(",NumberForm[N[evecs[[2,1]]],{5,4},NumberPadding->{"","0"}],", ",
    NumberForm[N[evecs[[2,2]]],{5,4},NumberPadding->{"","0"}],")"}];
  lamMinText=NumberForm[N[evals[[1]]],{5,4},NumberPadding->{"","0"}];
  lamMaxText=NumberForm[N[evals[[2]]],{5,4},NumberPadding->{"","0"}];

  infoBox=Framed[
    Grid[{
      {Style[TraditionalForm[Subscript["v","min"]],9.0,Bold,FontFamily->"Times",
        FontColor->$M3HessianDirection1Color],
       Style[Row[{" = ",vminText}],9.0,FontFamily->"Times"],
       Style[Row[{TraditionalForm[Subscript["\[Lambda]","min"]]," = ",lamMinText}],9.0,FontFamily->"Times"]},
      {Style[TraditionalForm[Subscript["v","max"]],9.0,Bold,FontFamily->"Times",
        FontColor->$M3HessianDirection2Color],
       Style[Row[{" = ",vmaxText}],9.0,FontFamily->"Times"],
       Style[Row[{TraditionalForm[Subscript["\[Lambda]","max"]]," = ",lamMaxText}],9.0,FontFamily->"Times"]}
    },Alignment->{{Left,Left,Left}},Spacings->{.45,.30}],
    Background->White,FrameStyle->Directive[GrayLevel[.82],AbsoluteThickness[.65]],
    FrameMargins->{{7,7},{5,5}},RoundingRadius->4
  ];

  Graphics[{
    Directive[$M3HessianDirection1Color,AbsoluteThickness[2.0]],Arrowheads[.026],
    Arrow[{phiStar-scales[[1]] evecs[[1]],phiStar+scales[[1]] evecs[[1]]}],
    Directive[$M3HessianDirection2Color,AbsoluteThickness[2.0]],Arrowheads[.026],
    Arrow[{phiStar-scales[[2]] evecs[[2]],phiStar+scales[[2]] evecs[[2]]}],
    Inset[infoBox,Scaled[{.965,.045}],{Right,Bottom}]
  }]
];

(* Internal renderer reused by the public contour figure. *)
M3HessianExactContoursPanelFromData[data_Association,imageSize_:430,padding_:{{54,18},{42,16}}] := Module[
  {exactDelta,phiStar,evecs,evals,halfWidth,levels,dirScales,contours,overlay},
  exactDelta=data["ExactDeltaFunction"]; phiStar=data["StationaryPoint"];
  evecs=data["Eigenvectors"]; evals=data["Eigenvalues"];
  halfWidth=data["LocalHalfWidth"]; levels=data["ContourLevels"];
  dirScales=0.98 data["EllipseSemiaxes"];

  contours=ContourPlot[
    exactDelta[phi2,phi3],
    {phi2,phiStar[[1]]-halfWidth,phiStar[[1]]+halfWidth},
    {phi3,phiStar[[2]]-halfWidth,phiStar[[2]]+halfWidth},
    Contours->levels,ContourShading->None,
    ContourStyle->Directive[GrayLevel[.22],AbsoluteThickness[1.08]],
    PlotPoints->56,MaxRecursion->2,Exclusions->None,
    Evaluate@M3HessianFrameOptions[imageSize,padding]
  ];

  overlay=Show[
    M3HessianPrincipalDirectionOverlay[phiStar,evecs,dirScales,evals],
    M3HessianStationaryMarker[phiStar]
  ];

  Show[contours,overlay]
];

(* Fig. 14(a) \[LongDash] exact contours and principal directions *)

M3HessianExactContoursFigure[
  lam_List,n_List,finite_Association,Q_?MatrixQ,tol_:10^-9,imageSize_:430
] := Module[{data},
  data = M3HessianArticleData[lam,n,finite,Q,tol];
  If[FailureQ[data], Return[data]];
  M3HessianExactContoursPanelFromData[data,imageSize]
];
(* Fig. 14(b) \[LongDash] exact action versus quadratic approximation; contour-comparison format. *)

M3HessianExactVsQuadraticFigure[lam_List,n_List,finite_Association,Q_?MatrixQ,tol_:10^-9,imageSize_:430] := Module[
  {data,exactDelta,quadDelta,phiStar,halfWidth,levels,exactContours,quadContours,overlay,legendBox},

  data=M3HessianArticleData[lam,n,finite,Q,tol];
  If[FailureQ[data],Return[data]];

  exactDelta=data["ExactDeltaFunction"]; quadDelta=data["QuadraticDeltaFunction"];
  phiStar=data["StationaryPoint"]; halfWidth=data["LocalHalfWidth"]; levels=data["ContourLevels"];

  exactContours=ContourPlot[
    exactDelta[phi2,phi3],
    {phi2,phiStar[[1]]-halfWidth,phiStar[[1]]+halfWidth},
    {phi3,phiStar[[2]]-halfWidth,phiStar[[2]]+halfWidth},
    Contours->levels,ContourShading->None,
    ContourStyle->Directive[$M3HessianExactColor,AbsoluteThickness[1.85]],
    PlotPoints->58,MaxRecursion->2,Exclusions->None,
    Evaluate@M3HessianFrameOptions[imageSize,{{54,18},{42,16}}]
  ];

  quadContours=ContourPlot[
    quadDelta[phi2,phi3],
    {phi2,phiStar[[1]]-halfWidth,phiStar[[1]]+halfWidth},
    {phi3,phiStar[[2]]-halfWidth,phiStar[[2]]+halfWidth},
    Contours->levels,ContourShading->None,
    ContourStyle->Directive[$M3HessianQuadraticColor,Dashed,AbsoluteThickness[1.65]],
    PlotPoints->42,MaxRecursion->1,Exclusions->None
  ];

  overlay=M3HessianStationaryMarker[phiStar];

  legendBox=Framed[
    Grid[{
      {Graphics[{Directive[$M3HessianExactColor,AbsoluteThickness[1.85]],Line[{{0,.5},{1,.5}}]},
         ImageSize->{28,10},PlotRange->{{0,1},{0,1}},Background->White],
       Style["exact action",9,FontFamily->"Times"]},
      {Graphics[{Directive[$M3HessianQuadraticColor,Dashed,AbsoluteThickness[1.65]],Line[{{0,.5},{1,.5}}]},
         ImageSize->{28,10},PlotRange->{{0,1},{0,1}},Background->White],
       Style["quadratic approximation",9,FontFamily->"Times"]}
    },Alignment->Left,Spacings->{.45,.25}],
    Background->White,FrameStyle->Directive[GrayLevel[.82],AbsoluteThickness[.65]],
    FrameMargins->{{7,7},{5,5}},RoundingRadius->4
  ];

  Show[exactContours,quadContours,overlay,
    Graphics[{Inset[legendBox,Scaled[{.965,.045}],{Right,Bottom}]}]]
];


(* real and imaginary parts of sigma with a common diverging scale. *)
M3SigmaComponentHeatmap[sigma_?MatrixQ,jvals_List,part_String,vmax_:Automatic]:=Module[
{A,K,scale,cf,leftTicks,bottomTicks,epilog,title,subtitle,legendLabel},
A=Switch[part,"Re",Re[N[sigma]],"Im",Im[N[sigma]],_,Re[N[sigma]]]; K=Length[A];

If[K==0,Return[Panel["The sigma matrix is empty."]]]; scale=If[vmax===Automatic,Max[10^-12,Max[Abs[Flatten[A]]]],Max[10^-12,N[vmax]]];
cf=Function[z,Blend[ {RGBColor[.14,.34,.62],RGBColor[.96,.97,.985],RGBColor[.76,.18,.14]}, Clip[(z+scale)/(2 scale),{0,1}] ]];
leftTicks=If[Length[jvals]===K, Table[{i,Style[TraditionalForm[jvals[[K-i+1]]],9.2,FontFamily->"Times"]},{i,K}],Automatic];
bottomTicks=If[Length[jvals]===K, Table[{i,Style[TraditionalForm[jvals[[i]]],9.2,FontFamily->"Times"]},{i,K}],Automatic];
epilog=If[K<=8, Flatten@Table[ Text[ Style[ NumberForm[A[[i,j]],If[Abs[A[[i,j]]]<10^-3,{3,2},{3,2}]], 7.8,FontFamily->"Times",
FontColor->If[Abs[A[[i,j]]]>0.55 scale,White,GrayLevel[.12]] ], {j,K-i+1} ],{i,K},{j,K}],{}];
title=If[part==="Re","Real part of the density matrix","Imaginary part of the density matrix"];
subtitle=If[part==="Re","Color = Re sigma_jj'","Color = Im sigma_jj'"]; legendLabel=If[part==="Re","Re sigma_jj'","Im sigma_jj'"];

M3ScientificResultPanel[ title, subtitle<>"; common symmetric scale around zero.", ArrayPlot[ Reverse[A], PlotRange->{-scale,scale},
ColorFunctionScaling->False, ColorFunction->cf, Mesh->All, MeshStyle->Directive[White,Opacity[.72],AbsoluteThickness[.55]],
Frame->True,Axes->False, FrameTicks->{{leftTicks,None},{bottomTicks,None}}, FrameLabel->{M3FrameLabel["j'"],M3FrameLabel["j"]},
FrameStyle->Directive[GrayLevel[.18],AbsoluteThickness[1.0]], Epilog->epilog, PlotLegends->Placed[
BarLegend[{cf,{-scale,scale}},LegendLabel->Style[legendLabel,9.2,FontFamily->"Times"],
LabelStyle->Directive[8.5,FontFamily->"Times"]],Right], Background->White, ImagePadding->{{48,64},{42,22}},
ImageSize->$M3PaperCanvasImageSize, ImageMargins->0, BaseStyle->{FontFamily->"Times",FontSize->10.2} ] ] ];

(* Orbital vs. physical isospin probabilities for the selected frame. *)
M3IsospinProbabilityFigure[jvals_List,pOrb_List,pPhys_List]:=Module[ { blue,orange,green,red,nj,xTicks,yMax,barHalf,barGap,
pOrbN,pPhysN,deltaN,deltaRawMax,deltaTickStep,deltaPlotMax, deltaYTicks,probBars,deltaBars,legendBox,topGraphic,bottomGraphic,
combinedGraphic }, If[Length[jvals]=!=Length[pOrb] || Length[jvals]=!=Length[pPhys], Return[Failure["ProbabilityDimensionMismatch",<|
"Message"->"j, P_orb, and P_phys must have the same length." |>]] ]; nj=Length[jvals];

If[nj==0,Return[Panel["There are no allowed isospin channels."]]]; pOrbN=N[pOrb]; pPhysN=N[pPhys]; deltaN=pPhysN - pOrbN;
blue=RGBColor[.16,.38,.64]; orange=RGBColor[.84,.40,.10]; green=RGBColor[.24,.58,.38]; red=RGBColor[.78,.27,.23];

(* Leave only the small amount of vertical room needed for the legend. *)
yMax=1.18 Max[Append[Join[pOrbN,pPhysN],10^-6]];

(* Symmetric and readable vertical scale for Delta P_j. *)
deltaRawMax=Max[Append[Abs[deltaN],10^-6]]; deltaTickStep=Which[ deltaRawMax<=.025,.01, deltaRawMax<=.06,.02, deltaRawMax<=.12,.025,
deltaRawMax<=.25,.05, deltaRawMax<=.50,.10, True,.20 ]; deltaPlotMax=deltaTickStep Ceiling[1.04 deltaRawMax/deltaTickStep];
deltaYTicks=Table[ {v,Style[NumberForm[Chop[v],{Infinity,2}],9.2,FontFamily->"Times"]}, {v,-deltaPlotMax,deltaPlotMax,deltaTickStep} ];
xTicks=Table[ {i,Style[TraditionalForm[jvals[[i]]],10,FontFamily->"Times"]}, {i,nj} ];

(* Slightly narrower bars with a visible separation inside each pair. *)
barHalf=.125; barGap=.052; probBars=Flatten @ Table[ { { EdgeForm[Directive[Darker[blue,.18],AbsoluteThickness[.65]]],
FaceForm[Directive[blue,Opacity[.92]]], Rectangle[{i - 2 barHalf - barGap,0},{i - barGap,pOrbN[[i]]}] }, {
EdgeForm[Directive[Darker[orange,.18],AbsoluteThickness[.65]]], FaceForm[Directive[orange,Opacity[.90]]],
Rectangle[{i + barGap,0},{i + 2 barHalf + barGap,pPhysN[[i]]}] } }, {i,nj} ]; deltaBars=Table[ {
EdgeForm[Directive[GrayLevel[.25],AbsoluteThickness[.5]]], FaceForm[Directive[If[deltaN[[i]]>=0,green,red],Opacity[.95]]],
Rectangle[{i - .22,0},{i + .22,deltaN[[i]]}] }, {i,nj} ]; legendBox=Framed[ SwatchLegend[ {blue,orange}, {
TraditionalForm[Subscript[P,"orb"][j]], TraditionalForm[Subscript[P,"phys"][j]] }, LegendLayout->"Row", LegendMarkerSize->11,
LabelStyle->Directive[9.5,FontFamily->"Times"], Spacings->.45 ], Background->White, FrameStyle->None, FrameMargins->{{5,5},{2,2}} ];
topGraphic=Graphics[ { probBars, Inset[legendBox,Scaled[{.50,.985}],{Center,Top}] }, Frame->True, Axes->False,
FrameStyle->Directive[GrayLevel[.18],AbsoluteThickness[.85]], FrameTicks->{{Automatic,None},{None,None}},
FrameTicksStyle->Directive[Black,9.4], FrameLabel->{None,M3FrameLabel[TraditionalForm[P[j]]]}, GridLines->None,
PlotRange->{{.5,nj + .5},{0,yMax}}, PlotRangePadding->{{0,0},{0,0}}, Background->White, AspectRatio->.46,
ImageSize->$M3PaperCanvasImageSize[[1]], ImagePadding->{{60,18},{8,8}}, BaseStyle->$M3BaseStyle ]; bottomGraphic=Graphics[ {
Directive[GrayLevel[.45],AbsoluteThickness[.8],Dashing[{.012,.012}]], Line[{{.5,0},{nj + .5,0}}], deltaBars }, Frame->True,
Axes->False, FrameStyle->Directive[GrayLevel[.18],AbsoluteThickness[.85]], FrameTicks->{{deltaYTicks,None},{xTicks,None}},
FrameTicksStyle->Directive[Black,9.4], FrameLabel->{ M3FrameLabel[TraditionalForm[j]],
M3FrameLabel[TraditionalForm[Subscript[\[CapitalDelta]P,j]]] }, PlotRange->{{.5,nj + .5},{-deltaPlotMax,deltaPlotMax}},
PlotRangePadding->{{0,0},{0,0}}, Background->White, AspectRatio->.20, ImageSize->$M3PaperCanvasImageSize[[1]],
ImagePadding->{{60,18},{40,6}}, BaseStyle->$M3BaseStyle ];

(* Build the two-panel figure first. *)
combinedGraphic=GraphicsColumn[ {topGraphic,bottomGraphic}, Spacings->-.02, Alignment->Center, Background->White,
ImageSize->$M3PaperCanvasImageSize, ImageMargins->0 ];

Graphics[ { White, Rectangle[{0,0},$M3PaperCanvasImageSize], Inset[ combinedGraphic, $M3PaperCanvasImageSize/2, Center,
$M3PaperCanvasImageSize ] }, PlotRange->{ {0,$M3PaperCanvasImageSize[[1]]}, {0,$M3PaperCanvasImageSize[[2]]} }, PlotRangePadding->0,
PlotRangeClipping->False, ImagePadding->0, ImageMargins->0, ImageSize->$M3PaperCanvasImageSize, AspectRatio->
N[$M3PaperCanvasImageSize[[2]]/$M3PaperCanvasImageSize[[1]]], Background->White, BaseStyle->$M3BaseStyle ] ];

M3SigmaHeatmapFigure[sigma_?MatrixQ,jvals_List]:=Module[ {m,n,vmax,colorStops,cf,ticksLeft,ticksBottom,legend,heat}, m=N[Abs[sigma]];
n=Length[m]; vmax=Max[10^-12,Max[m]];
colorStops={ RGBColor[.99,.995,1.00], RGBColor[.86,.91,.96], RGBColor[.62,.74,.85], RGBColor[.36,.56,.74], RGBColor[.12,.30,.50] };
cf=Function[z, Blend[colorStops,Clip[N[z]/vmax,{0,1}]] ]; ticksLeft=If[ Length[jvals]=!=n, Automatic, Table[{i,jvals[[n-i+1]]},{i,n}]
]; ticksBottom=If[ Length[jvals]=!=n, Automatic, Table[{i,jvals[[i]]},{i,n}] ]; legend=BarLegend[ {cf,{0,vmax}}, LegendLabel->Style[
Row[{"|",Subscript["\[Sigma]",Row[{"j",Superscript["j","\[Prime]"]}]],"|"}], 10, FontFamily->"Times" ],
LabelStyle->Directive[9,FontFamily->"Times"], LegendMarkerSize->{16,255} ];
heat=ArrayPlot[ Reverse[Transpose[m]], ColorFunctionScaling->False, ColorFunction->cf, Mesh->All,
MeshStyle->Directive[White,AbsoluteThickness[.55],Opacity[.8]], Frame->True, FrameTicks->{{ticksLeft,None},{ticksBottom,None}},
FrameLabel->{ {M3FrameLabel[Superscript["j","\[Prime]"]],None}, {M3FrameLabel["j"],None} }, Background->White,
PlotRangeClipping->False, ImagePadding->{{48,112},{42,18}}, Epilog->{ Inset[legend,Scaled[{1.105,.50}],Center] },
BaseStyle->{FontFamily->"Times",FontSize->11}, ImageSize->$M3PaperCanvasImageSize, ImageMargins->0 ]; heat ];

(* sigma spectrum with compact spectral observables. *)
M3SigmaSpectrumFigure[sigma_?MatrixQ,diag_Association]:=Module[ {ev,pts,tr,purity,entropy,minEig,maxEig,p,specSummary,blue},
ev=Sort[Re@N[Lookup[diag,"Eigenvalues",Eigenvalues[(sigma+ConjugateTranspose[sigma])/2]]],Greater]; tr=Re@Tr[N[sigma]];
purity=N[Lookup[diag,"Purity",Re@Tr[N[sigma . sigma]]]]; entropy=M3VonNeumannEntropy[sigma];
minEig=If[ev==={},0.,Min[ev]]; maxEig=If[ev==={},1.,Max[ev]]; pts=Transpose[{Range[Length[ev]],ev}]; blue=RGBColor[.15,.38,.65];
p=ListLinePlot[ pts,Joined->True, PlotStyle->Directive[blue,AbsoluteThickness[2.1]],
PlotMarkers->{Graphics[{EdgeForm[White],FaceForm[blue],Disk[]}],.026},
Filling->Axis,FillingStyle->Directive[RGBColor[.73,.83,.92],Opacity[.35]], Frame->True,Axes->False,
FrameLabel->{M3FrameLabel["index  k"],M3FrameLabel["lambda_k(sigma)"]},
GridLines->{Automatic,{0}},GridLinesStyle->Directive[GrayLevel[.92],AbsoluteThickness[.45]],
FrameStyle->Directive[GrayLevel[.18],AbsoluteThickness[1.0]],
PlotRange->{All,{Min[-.02 Max[1,maxEig],1.2 minEig],1.10 Max[10^-8,maxEig]}},
Background->White,ImageSize->$M3PaperCanvasImageSize,ImagePadding->{{54,16},{44,22}},ImageMargins->0,BaseStyle->{FontFamily->"Times",FontSize->10.2}
]; specSummary=Framed[ Grid[{ {Style["Tr sigma",Bold],NumberForm[tr,{8,6}],Style["Tr sigma^2",Bold],NumberForm[purity,{7,5}]},
{Style["S(sigma)",Bold],NumberForm[entropy,{7,5}],Style["lambda_min(sigma)",Bold],ScientificForm[minEig,3]}
},Alignment->{{Left,Right,Left,Right}},Spacings->{1.0,.55},BaseStyle->{FontFamily->"Times",FontSize->9.3}],
Background->RGBColor[.985,.991,.998],FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.65]],
FrameMargins->{{8,8},{5,5}},RoundingRadius->5 ];

M3ScientificResultPanel[ "Density-matrix spectrum", "Eigenvalues of sigma ordered from largest to smallest.",
Column[{p,specSummary},Alignment->Center,Spacings->.32] ] ];

M3DetectionComparisonFigure[comparison_Association,top_:18]:=Module[ {
phys,orb,physRows,orbRows,weights,pPhys,orbAssoc,pOrb,delta,score,order, selWeights,selPhys,selOrb,selDelta,nSel,labelAngle,xTicks,
blue,orange,red,negBlue,gray,barHalf,barGap,yMax,dMaxSel, barsMain,barsDelta,legendBoxMain,panelA,panelB,
coords,xr,yr,xSpan,ySpan,maxP,cfPhys,rPhys,cartanPhysPrims,legendPhys,panelC, dMax,cfDelta,rDelta,cartanDeltaPrims,legendDelta,panelD,
dTV,fcl,sumPhys,sumOrb,errPhys,errOrb,summary,labelWeight }, phys=comparison["Physical"]; orb=comparison["Orbital"];
physRows=Lookup[phys,"Results",{}]; orbRows=Lookup[orb,"Results",{}];

If[physRows==={},Return[Panel["There are no output probabilities."]]]; weights=Lookup[physRows,"nOut",{}];
pPhys=N@Lookup[physRows,"Probability",{}];
orbAssoc=AssociationThread[M3WtKey/@Lookup[orbRows,"nOut",{}],N@Lookup[orbRows,"Probability",{}]];
pOrb=N@Lookup[orbAssoc,M3WtKey/@weights,0.]; delta=pPhys-pOrb; score=MapThread[Max,{pPhys,pOrb}];
order=Take[Reverse@Ordering[score],Min[top,Length[score]]]; selWeights=weights[[order]]; selPhys=pPhys[[order]]; selOrb=pOrb[[order]];
selDelta=delta[[order]]; nSel=Length[order]; labelAngle=If[nSel>10,Pi/4,0];

labelWeight[w_List]:=Row[{"(",w[[1]],",",w[[2]],",",w[[3]],")"}];
xTicks=Table[{i,Rotate[Style[labelWeight[selWeights[[i]]],7.6,FontFamily->"Times"],labelAngle]},{i,nSel}]; blue=RGBColor[.19,.42,.68];
orange=RGBColor[.85,.39,.10]; red=RGBColor[.78,.18,.12]; negBlue=RGBColor[.16,.38,.67]; gray=GrayLevel[.78];

(* bars grouped by event, selected with max(P_orb,P_phys). *)
barHalf=.13; barGap=.05; yMax=1.10 Max[Append[Join[selPhys,selOrb],10^-8]]; barsMain=Flatten@Table[ { Tooltip[
{EdgeForm[Directive[Darker[blue,.22],AbsoluteThickness[.65]]],FaceForm[Directive[blue,Opacity[.95]]],
Rectangle[{i-2 barHalf-barGap/2,0},{i-barGap/2,selOrb[[i]]}]},
Column[{Style["Orbital",Bold],Row[{"n' = ",selWeights[[i]]}],Row[{"P = ",NumberForm[selOrb[[i]],{7,5}]}]}] ], Tooltip[
{EdgeForm[Directive[Darker[orange,.22],AbsoluteThickness[.65]]],FaceForm[Directive[orange,Opacity[.95]]],
Rectangle[{i+barGap/2,0},{i+2 barHalf+barGap/2,selPhys[[i]]}]},
Column[{Style["Physical",Bold],Row[{"n' = ",selWeights[[i]]}],Row[{"P = ",NumberForm[selPhys[[i]],{7,5}]}]}] ] },{i,nSel}];
legendBoxMain=Framed[ Grid[{ {
Graphics[{EdgeForm[Directive[Darker[blue,.20],AbsoluteThickness[.55]]],FaceForm[blue],Rectangle[{0,0},{1,1}]},ImageSize->{12,12},PlotRange->{{0,1},{0,1}},Background->None],
Style[TraditionalForm[Subscript["P","orb"]["n'"]],9.8,FontFamily->"Times"] }, {
Graphics[{EdgeForm[Directive[Darker[orange,.20],AbsoluteThickness[.55]]],FaceForm[orange],Rectangle[{0,0},{1,1}]},ImageSize->{12,12},PlotRange->{{0,1},{0,1}},Background->None],
Style[TraditionalForm[Subscript["P","phys"]["n'"]],9.8,FontFamily->"Times"] } },Alignment->{Left,Center},Spacings->{.45,.22}],
Background->White, FrameStyle->Directive[GrayLevel[.80],AbsoluteThickness[.55]], FrameMargins->{{6,6},{4,4}},RoundingRadius->3 ];
panelA=M3ScientificResultPanel[ "Orbital and physical detection distributions",
"The outputs with the largest max(P_orb,P_phys) are shown, without favoring either distribution.", Graphics[ { barsMain,
Inset[legendBoxMain,Scaled[{.86,.83}],Center] }, Frame->True,Axes->False,PlotRange->{{.45,nSel+.55},{0,yMax}},
FrameTicks->{{Automatic,None},{xTicks,None}}, FrameLabel->{M3FrameLabel["output occupation  n'"],M3FrameLabel["P(n')"]},
GridLines->None, FrameStyle->Directive[GrayLevel[.18],AbsoluteThickness[1.0]],
FrameTicksStyle->Directive[GrayLevel[.18],8.8,FontFamily->"Times"], Background->White, ImageSize->$M3PaperCanvasImageSize,
AspectRatio->.62, ImagePadding->{{58,16},{If[nSel>10,72,52],16}}, ImageMargins->0, BaseStyle->{FontFamily->"Times",FontSize->10.2} ] ];

(* event-by-event difference with a diverging scale. *)
dMaxSel=Max[10^-12,Max[Abs[selDelta]]]; barsDelta=Table[
With[{d=selDelta[[i]],c=Which[selDelta[[i]]>10^-13,red,selDelta[[i]]<-10^-13,negBlue,True,gray]}, Tooltip[
{EdgeForm[Directive[Darker[c,.25],AbsoluteThickness[.72]]],FaceForm[Directive[c,Opacity[.95]]],
Rectangle[{i-.26,Min[0,d]},{i+.26,Max[0,d]}]}, Column[{Row[{"n' = ",selWeights[[i]]}],Row[{"Delta P = ",ScientificForm[d,3]}]}] ]
],{i,nSel}]; panelB=M3ScientificResultPanel[ "Probability difference by event",
"Red: P_phys>P_orb. Blue: P_phys<P_orb. The horizontal line corresponds to Delta P=0.", Graphics[
{barsDelta,Directive[GrayLevel[.35],Dashed,AbsoluteThickness[1.0]],Line[{{.45,0},{nSel+.55,0}}]},
Frame->True,Axes->False,PlotRange->{{.45,nSel+.55},1.16 {-dMaxSel,dMaxSel}}, FrameTicks->{{Automatic,None},{xTicks,None}},
FrameLabel->{M3FrameLabel["output occupation  n'"],M3FrameLabel["Delta P(n')"]}, GridLines->None,
FrameStyle->Directive[GrayLevel[.18],AbsoluteThickness[1.0]], FrameTicksStyle->Directive[GrayLevel[.18],8.8,FontFamily->"Times"],
Background->White, ImageSize->$M3PaperCanvasImageSize, AspectRatio->.62, ImagePadding->{{58,16},{If[nSel>10,72,52],16}},
ImageMargins->0, BaseStyle->{FontFamily->"Times",FontSize->10.2} ] ];

(* Cartan coordinates of all outputs. *)
coords=N[{(#[[1]]-#[[2]])/Sqrt[2],(#[[1]]+#[[2]]-2 #[[3]])/Sqrt[6]}&/@weights]; xr=MinMax[coords[[All,1]]]; yr=MinMax[coords[[All,2]]];
xSpan=Max[10^-8,xr[[2]]-xr[[1]]]; ySpan=Max[10^-8,yr[[2]]-yr[[1]]]; xr=xr+{-0.07 xSpan,0.07 xSpan}; yr=yr+{-0.07 ySpan,0.07 ySpan};
(* Cartan map of P_phys. *)
maxP=Max[10^-14,Max[pPhys]]; cfPhys=ColorData["LakeColors"];

rPhys[p_?NumericQ]:=.050+.105 Sqrt[Max[0,p]/maxP]; cartanPhysPrims=MapThread[ Function[{pt,w,p}, Tooltip[
{EdgeForm[Directive[GrayLevel[.30],AbsoluteThickness[.65]]], FaceForm[cfPhys[.10+.84 Clip[p/maxP,{0,1}]]],Disk[pt,rPhys[p]]},
Column[{Row[{"n' = ",w}],Row[{"P_phys = ",NumberForm[p,{7,5}]}]}] ] ],{coords,weights,pPhys}];
legendPhys=BarLegend[{"LakeColors",{0,maxP}},LegendLabel->Style["P_phys(n')",8.8,FontFamily->"Times"],LabelStyle->Directive[8.3,FontFamily->"Times"]];
panelC=M3ScientificResultPanel[ "Physical distribution on the Cartan diagram", "The color and size of each node represent P_phys(n').",
Legended[ Graphics[ cartanPhysPrims, Frame->True,Axes->False,AspectRatio->1,PlotRange->{xr,yr},
FrameLabel->{M3FrameLabel["X_C"],M3FrameLabel["Y_C"]},
GridLines->Automatic,GridLinesStyle->Directive[GrayLevel[.93],AbsoluteThickness[.42]],
FrameStyle->Directive[GrayLevel[.18],AbsoluteThickness[1.0]], FrameTicksStyle->Directive[GrayLevel[.18],8.7,FontFamily->"Times"],
Background->White,ImageSize->430,ImagePadding->{{50,18},{42,18}},BaseStyle->{FontFamily->"Times",FontSize->10.2} ],
Placed[legendPhys,Right] ] ];

(* Cartan map of Delta P.  One physical color map is used everywhere:
negative -> red, zero -> white, positive -> blue. *)
dMax=Max[10^-14,Max[Abs[delta]]]; cfDelta=Function[d, Which[ d<=0, Blend[ {RGBColor[.78,.16,.12],GrayLevel[.985]},
Clip[(d+dMax)/dMax,{0,1}] ], True, Blend[ {GrayLevel[.985],RGBColor[.13,.36,.68]}, Clip[d/dMax,{0,1}] ] ] ];

(* Slightly larger Cartan markers for article readability.
Their relative size still scales with Sqrt[|Delta P|]. *)
rDelta[d_?NumericQ]:=.075+.155 Sqrt[Abs[d]/dMax]; cartanDeltaPrims=MapThread[ Function[{pt,w,d}, Tooltip[ {
FaceForm[White],Disk[pt,rDelta[d]+.012], EdgeForm[Directive[GrayLevel[.25],AbsoluteThickness[.55]]],
FaceForm[cfDelta[d]],Disk[pt,rDelta[d]] }, Column[{Row[{"n' = ",w}],Row[{"Delta P = ",ScientificForm[d,3]}]}] ]
],{coords,weights,delta}];

(* Use exactly the same physical color function for both the Cartan points
and the legend, with a symmetric readable scale around zero. *)
legendDelta=BarLegend[ { (cfDelta[#]&), {-dMax,dMax} }, Ticks->Table[ { val, If[ val==0, Style["0",8.9,FontFamily->"Times"],
Style[ScientificForm[val,2],8.9,FontFamily->"Times"] ] }, {val,{-dMax,-dMax/2,0,dMax/2,dMax}} ], LegendLabel->Style[
TraditionalForm@HoldForm[\[CapitalDelta] P[Superscript[n,"\[Prime]"]]], 9.8,FontFamily->"Times" ],
LabelStyle->Directive[8.9,FontFamily->"Times"], LegendMarkerSize->{17,245} ]; panelD=M3ScientificResultPanel[
"Physical redistribution on the Cartan diagram", "Red: physical decrease; blue: physical increase; white: Delta P = 0.", Graphics[ {
cartanDeltaPrims, Inset[legendDelta,Scaled[{1.15,.50}],Center] }, Frame->True,Axes->False,AspectRatio->1,PlotRange->{xr,yr},
PlotRangePadding->Scaled[.025], FrameLabel->{ Style[TraditionalForm@HoldForm[Subscript[X,C]],11.2,FontFamily->"Times"],
Style[TraditionalForm@HoldForm[Subscript[Y,C]],11.2,FontFamily->"Times"] }, GridLines->None,
FrameStyle->Directive[GrayLevel[.18],AbsoluteThickness[.95]], FrameTicksStyle->Directive[GrayLevel[.18],8.9,FontFamily->"Times"],
Background->White, ImageSize->455, ImagePadding->{{54,102},{40,16}}, ImageMargins->0, BaseStyle->{FontFamily->"Times",FontSize->10.2} ]
]; dTV=.5 Total[Abs[delta]]; fcl=Total[Sqrt[MapThread[Max[0,#1] Max[0,#2]&,{pPhys,pOrb}]]];
sumPhys=N@Lookup[phys,"TotalProbability",Total[pPhys]]; sumOrb=N@Lookup[orb,"TotalProbability",Total[pOrb]];
errPhys=N@Lookup[phys,"NormalizationError",Abs[sumPhys-1]]; errOrb=N@Lookup[orb,"NormalizationError",Abs[sumOrb-1]];
summary=M3ScientificResultPanel[ "Global comparison summary",
"D_TV=0 corresponds to identical distributions; F_cl=1 corresponds to perfect classical overlap.", Grid[{
{Style[Subscript["D","TV"],Bold],NumberForm[dTV,{8,6}],Style[Subscript["F","cl"],Bold],NumberForm[fcl,{8,6}]},
{Style["Sum P_phys",Bold],NumberForm[sumPhys,{10,8}],Style["physical error",Bold],ScientificForm[errPhys,3]},
{Style["Sum P_orb",Bold],NumberForm[sumOrb,{10,8}],Style["orbital error",Bold],ScientificForm[errOrb,3]}
},Alignment->{{Left,Right,Left,Right}},Spacings->{1.2,.60}, Background->{1->{RGBColor[.975,.985,.996]}},Frame->All,
BaseStyle->{FontFamily->"Times",FontSize->9.6}] ];

Column[{ Grid[{{panelA,panelB},{panelC,panelD}},Alignment->Top,Spacings->{.65,.55}], summary },Alignment->Center,Spacings->.75] ];

M3RunningRanks[C_?MatrixQ,tol_:10^-8]:=Module[{K},K=Dimensions[C][[2]]; Table[M3CoordinateRank[C[[All,1 ;; k]],tol],{k,1,K}]];


(* SCIENTIFIC PRESENTATION: MULTIPLICITY / GRAM / CARTAN *)
M3SU3I3YCoordinates[W_List]:=N[{ (W[[1]] - W[[2]])/2, (W[[1]] + W[[2]] - 2 W[[3]])/3
}];

M3ScientificResultPanel[title_,subtitle_,body_]:=Panel[ Column[{
  Style[title,12.4,Bold,FontFamily->"Times",FontColor->RGBColor[.12,.25,.42]],
  If[subtitle==="",Nothing,Style[subtitle,8.9,GrayLevel[.38],FontFamily->"Times"]], body },Alignment->Center,Spacings->.34],
  Background->White, FrameMargins->{{10,10},{8,8}}, FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.7]]
];

(* SU(3) weights in physical coordinates I3,Y, with GT multiplicity. *)
M3GTWeightDiagramArticle[data_Association,selectedN_:None]:=Module[
  {rows,pts,mults,maxm,xr,yr,selected,selectedPt,prims,overlay,radius,color,textSize},

rows=data["Rows"]; If[rows==={},Return[Panel["There are no weights in the irrep."]]];

pts=M3SU3I3YCoordinates /@ Lookup[rows,"W"]; mults=Lookup[rows,"MultiplicityGT"]; maxm=Max[1,Max[mults]];
  xr=MinMax[pts[[All,1]]] + {-0.9,0.9}; yr=MinMax[pts[[All,2]]] + {-0.9,0.9}; radius[m_]:=.070 + .095 Sqrt[N[m/maxm]];
  color[m_]:=Blend[ { RGBColor[.94,.965,.99], RGBColor[.73,.83,.91], RGBColor[.42,.61,.78] }, Rescale[N[m],{1,maxm},{.05,1.}] ];
  textSize[m_]:=9.6 + .55 Sqrt[N[m/maxm]];

prims=Flatten @ MapThread[ Function[{row,pt,m}, { Tooltip[ {
(* very faint outer halo to separate the node from the background *)
FaceForm[Directive[White,Opacity[.98]]], EdgeForm[None], Disk[pt,1.13 radius[m]],

(* main node *)
FaceForm[Directive[color[m],Opacity[.96]]], EdgeForm[Directive[RGBColor[.43,.52,.61],AbsoluteThickness[.75]]], Disk[pt,radius[m]] },
  Column[{ Style["SU(3) weight",Bold], Row[{"n = ",row["n"]}], Row[{"W = ",row["W"]}], Row[{"m_GT(W) = ",m}] },Spacings->.18] ],

(* The number is drawn last, with no lines or contours over it. *)
Text[ Style[ ToString[m], textSize[m], Bold, FontFamily->"Helvetica", FontColor->RGBColor[.08,.11,.14] ], pt ] } ], {rows,pts,mults}
  ];

selected=If[ selectedN===None, Missing["NotSelected"], SelectFirst[rows,# ["n"]===selectedN&,Missing["NotFound"]] ]; selectedPt=If[
  AssociationQ[selected], M3SU3I3YCoordinates[selected["W"]], Missing["NotSelected"] ];

  (* Double ring for the selected weight, placed outside the node so as not to cover the number. *)
overlay=If[ ListQ[selectedPt], { Directive[White,AbsoluteThickness[4.2]],Circle[selectedPt,.252],
  Directive[RGBColor[.80,.15,.08],AbsoluteThickness[2.65]],Circle[selectedPt,.252] }, {} ];

M3ScientificResultPanel[ "Weight structure and theoretical multiplicity",
  "The integer centered at each weight is m_GT(W); node size and intensity encode the same multiplicity.", Graphics[
  Join[prims,overlay], Frame->True, Axes->False, AspectRatio->1, PlotRange->{xr,yr}, PlotRangePadding->Scaled[.025],
  GridLines -> None, FrameLabel->{
  Style[TraditionalForm[Subscript["I",3]],11.8,Bold,FontFamily->"Times"], Style[TraditionalForm["Y"],11.8,Bold,FontFamily->"Times"]
  }, FrameStyle->Directive[GrayLevel[.18],AbsoluteThickness[1.05]],
  FrameTicksStyle->Directive[GrayLevel[.16],9.4,FontFamily->"Times"], Background->White, ImagePadding->{{52,18},{44,20}},
  ImageSize->{500,380}, BaseStyle->{FontFamily->"Times",FontSize->10.2}, PerformanceGoal->"Quality" ] ]
];

(* normalized singular spectrum of C. *)
M3NormalizedSingularSpectrumArticle[C_?MatrixQ,tol_:10^-8]:=Module[ {sv,norm,vals,K,floor,plotVals,ymin,rank,activePts,inactivePts},

sv=N[SingularValueList[C]];
  If[sv==={} || Max[sv]<=0,Return[M3ScientificResultPanel["Singular spectrum of C","","There are no nonzero singular values."]]];

norm=First[sv]; vals=sv/norm; K=Length[vals]; floor=Min[10^-14,tol/100]; plotVals=Max[#,floor]& /@ vals; ymin=Min[floor,tol/10];
  rank=Count[vals,x_/;x>tol]; activePts=Select[Transpose[{Range[K],plotVals}],#[[2]]>tol&];
  inactivePts=Select[Transpose[{Range[K],plotVals}],#[[2]]<=tol&];

M3ScientificResultPanel[ "Normalized singular spectrum of C", "Modes above rankTol contribute to the numerical rank of C.", Show[
  ListPlot[ Transpose[{Range[K],plotVals}],Joined->True, PlotStyle->Directive[RGBColor[.15,.34,.56],AbsoluteThickness[2.35]],
  Frame->True,Axes->False,ScalingFunctions->{None,"Log10"}, PlotRange->{{.65,K+.35},{ymin,1.35}}, FrameLabel->{
  Style["singular index  k",10.4,FontFamily->"Times"],
  Style[Row[{Subscript["\[Sigma]","k"],"(C) / ",Subscript["\[Sigma]","1"],"(C)"}],10.4,FontFamily->"Times"] },
  FrameStyle->Directive[GrayLevel[.20],AbsoluteThickness[1.0]], FrameTicksStyle->Directive[GrayLevel[.18],9.0,FontFamily->"Times"],
  GridLines->{Automatic,Automatic}, GridLinesStyle->Directive[GrayLevel[.94],AbsoluteThickness[.42]], Epilog->{
  {Directive[RGBColor[.82,.23,.12],Dashed,AbsoluteThickness[1.55]],Line[{{.55,tol},{K+.45,tol}}]}, Inset[
  Framed[Style[Row[{"rank C = ",rank,"   |   rankTol = ",ScientificForm[tol,1]}],8.5,FontFamily->"Times",FontColor->GrayLevel[.22]],
  Background->White,FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.65]],RoundingRadius->4,FrameMargins->{{6,6},{3,3}}],
  Scaled[{.70,.18}] ] },
  Background->White,ImagePadding->{{56,14},{42,18}},ImageSize->440,BaseStyle->{FontFamily->"Times",FontSize->10.2} ],
  ListPlot[activePts,PlotStyle->Directive[RGBColor[.10,.39,.66]],PlotMarkers->{Graphics[{EdgeForm[White],Disk[]}],8}],
  If[inactivePts==={},Graphics[{}],ListPlot[inactivePts,PlotStyle->Directive[RGBColor[.80,.22,.12]],PlotMarkers->{Graphics[{EdgeForm[White],Disk[]}],8}]]
  ] ]
];
(* convergence of rank(C_K)/m_GT. *)
M3RankConvergenceArticle[C_?MatrixQ,mGT_Integer?NonNegative,tol_:10^-8]:=Module[ {ranks,K,ratio,ymax,firstComplete,pts},

ranks=M3RunningRanks[C,tol]; K=Length[ranks];
  If[K==0 || mGT<=0,Return[M3ScientificResultPanel["Bank convergence","","There are not enough data for convergence."]]];

ratio=N[ranks/mGT]; ymax=Max[1.08,1.05 Max[ratio]]; firstComplete=SelectFirst[Range[K],ratio[[#]]>=1-10^-12&,Missing["NotReached"]];
  pts=Transpose[{Range[K],ratio}];

M3ScientificResultPanel[
 "Frame-bank convergence",
 "Complete generation is reached when the normalized rank reaches 1.",
 Show[
  ListLinePlot[
   pts,
   PlotStyle -> Directive[RGBColor[.08, .45, .30], AbsoluteThickness[2.55]],
   PlotMarkers -> {Graphics[{EdgeForm[White], FaceForm[RGBColor[.08, .45, .30]], Disk[]}], 7.5},
   Frame -> True,
   Axes -> False,
   PlotRange -> {{.75, Max[1.25, K + .25]}, {0, ymax}},
   FrameLabel -> {
  Style[
   Row[{
     "cumulative number of frames ",
     Style["K", Italic]
   }],
   10.4,
   FontFamily -> "Times"
  ],

  Style[
 Row[{
   Style["rank", FontSlant -> "Plain"],
   "(",
   Subscript[
    Style["C", Italic],
    Style["K", Italic]
   ],
   ") / ",
   Subscript[
    Style["m", Italic],
    Style["GT", FontSlant -> "Plain"]
   ]
 }],
 10.4,
 FontFamily -> "Times"
]
},
   FrameStyle -> Directive[GrayLevel[.20], AbsoluteThickness[1.0]],
   FrameTicksStyle -> Directive[GrayLevel[.18], 9.0, FontFamily -> "Times"],
   GridLines -> None,
   Epilog -> Join[
     {
       {
        Directive[
 GrayLevel[.45],
 Dashed,
 AbsoluteThickness[1.1]
],
        Line[{{.65, 1}, {K + .35, 1}}]
       }
     },
     If[
      IntegerQ[firstComplete],
      {
       Directive[RGBColor[.12, .32, .56], AbsoluteThickness[1.2], Dashing[{.015, .01}]],
       Line[{{firstComplete, 0}, {firstComplete, 1}}],
       Inset[
 Style[
  Row[{
    Superscript[Style["K", Italic], "*"],
    " = ",
    firstComplete
  }],
  9,
  FontFamily -> "Times",
  FontColor -> RGBColor[.12,.32,.56]
 ],
 {firstComplete + .18, .12},
 {Left, Bottom}
]
      },
      {}
     ]
   ],
   Background -> White,
   ImagePadding -> {{56, 14}, {42, 18}},
   ImageSize -> {500,380},
   BaseStyle -> {FontFamily -> "Times", FontSize -> 10.2}
  ]
 ]
]
];
M3MultiplicityDeficitMapArticle[data_Association,selectedN_:None]:=Module[ { rows,pts,selected,selectedPt,xr,yr,delta,maxDelta,
  zeroRows,posRows,missingRows,zeroPrims,posPrims,missingPrims, selectedOverlay,legend,plot,zeroColor,posColor,edgeColor,zeroR,posR
  },

rows=data["Rows"]; If[rows==={},Return[Panel["There are no weights in the irrep."]]];

pts=M3SU3I3YCoordinates /@ Lookup[rows,"W"]; xr=MinMax[pts[[All,1]]] + {-0.9,0.9}; yr=MinMax[pts[[All,2]]] + {-0.9,0.9};

delta[row_]:=Module[{r=Lookup[row,"RankC",Missing["NoRank"]]}, If[IntegerQ[r],Max[0,row["MultiplicityGT"]-r],Missing["NoRank"]] ];

maxDelta=Max[1,Max[Replace[delta /@ rows,_Missing->0,{1}]]]; zeroColor=RGBColor[.73,.84,.93]; edgeColor=RGBColor[.28,.45,.61];
  zeroR=.070; posR[d_]:=.105 + .030 Sqrt[N[d/maxDelta]]; posColor[d_]:=Blend[
  {RGBColor[.97,.68,.35],RGBColor[.89,.30,.14],RGBColor[.65,.07,.06]}, Rescale[N[d],{1,maxDelta},{.18,1.}] ];

zeroRows=Pick[Transpose[{rows,pts}],Map[TrueQ[delta[#]===0]&,rows]];
  posRows=Pick[Transpose[{rows,pts}],Map[IntegerQ[delta[#]]&&delta[#]>0&,rows]];
  missingRows=Pick[Transpose[{rows,pts}],Map[MissingQ[delta[#]]&,rows]];

zeroPrims=Flatten @ Map[ Function[rp, With[{row=rp[[1]],pt=rp[[2]]}, { Tooltip[ { FaceForm[Directive[zeroColor,Opacity[.96]]],
  EdgeForm[Directive[edgeColor,AbsoluteThickness[.72]]], Disk[pt,zeroR], FaceForm[None],
  EdgeForm[Directive[White,Opacity[.52],AbsoluteThickness[.55]]], Circle[pt,.73 zeroR] }, Column[{
  Style["Complete generation",Bold,FontColor->RGBColor[.12,.36,.56]], Row[{"n = ",row["n"]}],Row[{"m_GT = ",row["MultiplicityGT"]}],
  Row[{"rank C = ",row["RankC"]}],"Delta = 0" },Spacings->.16] ] } ] ],zeroRows ];

posPrims=Flatten @ Map[ Function[rp, With[{row=rp[[1]],pt=rp[[2]],d=delta[rp[[1]]]}, { Tooltip[ {
  FaceForm[Directive[posColor[d],Opacity[.98]]], EdgeForm[Directive[RGBColor[.55,.07,.05],AbsoluteThickness[1.15]]],
  Disk[pt,posR[d]], FaceForm[None], EdgeForm[Directive[White,Opacity[.45],AbsoluteThickness[.7]]], Circle[pt,.76 posR[d]] },
  Column[{ Style["Generation deficit",Bold,FontColor->RGBColor[.67,.09,.06]],
  Row[{"n = ",row["n"]}],Row[{"m_GT = ",row["MultiplicityGT"]}], Row[{"rank C = ",row["RankC"]}],Row[{"Delta = ",d}]
  },Spacings->.16] ], Text[Style[ToString[d],8.8,Bold,FontFamily->"Times",FontColor->White],pt] } ] ],posRows ];

missingPrims=Flatten @ Map[ Function[rp, With[{row=rp[[1]],pt=rp[[2]]}, { Tooltip[
  {FaceForm[White],EdgeForm[Directive[GrayLevel[.52],AbsoluteThickness[.95]]],Circle[pt,.062]},
  Column[{Style["Unresolved rank",Bold,FontColor->GrayLevel[.38]],Row[{"n = ",row["n"]}]}] ] } ] ],missingRows ];

selected=If[selectedN===None,Missing["NotSelected"],SelectFirst[rows,# ["n"]===selectedN&,Missing["NotFound"]]];
  selectedPt=If[AssociationQ[selected],M3SU3I3YCoordinates[selected["W"]],Missing["NotSelected"]];
  selectedOverlay=If[ListQ[selectedPt],{ Directive[White,AbsoluteThickness[4.2]],Circle[selectedPt,.27],
  Directive[RGBColor[.08,.28,.53],AbsoluteThickness[2.55]],Circle[selectedPt,.27] },{}];

legend=Framed[ Grid[{{
  Graphics[{FaceForm[zeroColor],EdgeForm[Directive[edgeColor,AbsoluteThickness[.75]]],Disk[{.5,.5},.24]},PlotRange->{{0,1},{0,1}},ImageSize->{22,22},Background->White],
  Style["Delta = 0",9.0,FontFamily->"Times"],Spacer[12],
  Graphics[{FaceForm[RGBColor[.86,.22,.11]],EdgeForm[Directive[RGBColor[.55,.07,.05],AbsoluteThickness[.85]]],Disk[{.5,.5},.27]},PlotRange->{{0,1},{0,1}},ImageSize->{22,22},Background->White],
  Style["Delta > 0",9.0,FontFamily->"Times"],Spacer[12],
  Graphics[{FaceForm[White],EdgeForm[Directive[GrayLevel[.52],AbsoluteThickness[.9]]],Circle[{.5,.5},.23]},PlotRange->{{0,1},{0,1}},ImageSize->{22,22},Background->White],
  Style["unresolved rank",9.0,FontFamily->"Times"] }},Alignment->Center,Spacings->{.50,.2}], Background->RGBColor[.995,.997,.999],
  FrameStyle->Directive[GrayLevel[.85],AbsoluteThickness[.65]],RoundingRadius->5,FrameMargins->{{9,9},{5,5}} ];

plot=Graphics[ Join[zeroPrims,missingPrims,posPrims,selectedOverlay], Frame->True,Axes->False,AspectRatio->1,
  PlotRange->{xr,yr},PlotRangePadding->Scaled[.025],
  GridLines->Automatic,GridLinesStyle->Directive[GrayLevel[.94],AbsoluteThickness[.42]], FrameLabel->{
  Style[TraditionalForm[Subscript["I",3]],11.5,Bold,FontFamily->"Times"], Style[TraditionalForm["Y"],11.5,Bold,FontFamily->"Times"]
  }, FrameStyle->Directive[GrayLevel[.20],AbsoluteThickness[1.0]],
  FrameTicksStyle->Directive[GrayLevel[.18],9.1,FontFamily->"Times"],
  Background->White,ImagePadding->{{50,16},{42,18}},ImageSize->440,BaseStyle->{FontFamily->"Times",FontSize->10.2} ];

M3ScientificResultPanel[ "Generation deficit on the weight diagram",
  "Delta(W)=m_GT(W)-rank C(W). Red nodes identify exclusively the weights with positive deficit.",
  Column[{plot,legend},Alignment->Center,Spacings->.34] ]
];
M3MultiplicitySummaryArticle[finite_Association,iso_Association,multDiag_Association,gramDiag_Association,records_List]:=Module[
  {mGT,rankC,deficit,jtext,status}, mGT=multDiag["MultiplicityGT"]; rankC=gramDiag["RankC"];
  deficit=If[IntegerQ[rankC],mGT-rankC,Missing["NoRank"]]; jtext=Row[Riffle[ToString /@ iso["JValues"],", "]];
  status=If[IntegerQ[rankC] && rankC==mGT, Style["complete generation",Darker[Green],Bold],
  Style["generation not certified",Darker[Orange],Bold] ];

Panel[ Grid[ { {Style["weight n",Bold],finite["nN"],Style["W",Bold],finite["W"],Style["JValues",Bold],jtext},
  {Style["m_GT",Bold],mGT,Style["rank C = rank G",Bold],rankC,Style["Delta",Bold],deficit},
  {Style["frames used",Bold],Length[records],Style["status",Bold],status,SpanFromLeft,SpanFromLeft} }, Alignment->Left,
  Spacings->{1.0,.5}, Dividers->{False,{False,True,True,False}}, Background->{1->{RGBColor[.96,.975,.99]}},
  BaseStyle->{FontFamily->"Times",FontSize->9.5} ], Background->White, FrameMargins->{{9,9},{6,6}},
  FrameStyle->Directive[GrayLevel[.83],AbsoluteThickness[.7]] ]
];

(* Appendix F / Fig. 7 \[LongDash] multiplicity and rank figures *)
M3CartanCoordinates3[W_List]:=N[{W[[1]]-W[[2]],(W[[1]]+W[[2]]-2 W[[3]])/Sqrt[3]}];

(* Enumeration of ALL GT weights and multiplicities in a single pass. It is much faster than recounting patterns for each composition. *)
M3AllWeightMultiplicitiesGT[lamN_List]:=Module[{l1,l2,l3,reaped,weights,keys,counts,unique},
  If[Length[lamN]=!=3||!VectorQ[lamN,IntegerQ]||!(lamN[[1]]>=lamN[[2]]>=lamN[[3]]>=0),
  Return[Failure["InvalidPartition",<|"Message"->"Lambda must be a nonincreasing integer partition."|>]] ]; {l1,l2,l3}=lamN;
  reaped=Reap[ Do[ Do[ Do[ Sow[{nu,mu1+mu2-nu,Total[lamN]-mu1-mu2}], {nu,mu2,mu1}], {mu2,l3,l2}], {mu1,l2,l1}] ][[2]];
  weights=If[reaped==={},{},First[reaped]]; If[weights==={},Return[{}]]; keys=M3WtKey/@weights; counts=Counts[keys];
  unique=DeleteDuplicates[weights]; SortBy[ Table[<|"n"->w,"MultiplicityGT"->Lookup[counts,M3WtKey[w],0]|>,{w,unique}],
  {#["n"][[3]],#["n"][[2]],#["n"][[1]]}&]
];

(* Empirical rank for ONE weight. It first uses a moderate bank and increases sampling only if the theoretical multiplicity has not yet been reached. *)
M3WeightRankSample[lamN_List,nN_List,mGT_Integer,maxFrames_Integer,gridN_Integer,includeConjugate_,rankTol_]:=Module[
  {Nrep,finite,iso,try,records,coords,C,diag,rank=Missing["NoRank"],used=0,status="no frames",localFrames,localGrid},
  Nrep=Total[lamN]; finite=M3FiniteDataFromIntegerData[lamN,nN];
  If[FailureQ[finite],Return[<|"Rank"->rank,"FramesUsed"->0,"Status"->"invalid finite data"|>]]; iso=M3IsospinData[finite];
  If[FailureQ[iso]||iso["Multiplicity"]==0,Return[<|"Rank"->0,"FramesUsed"->0,"Status"->"zero multiplicity"|>]];

try[fcount_,gcount_]:=Module[{rr,cc,cx,wd,ex,mat,gd},
  rr=M3CandidateFrames[N[lamN/Nrep],N[nN/Nrep],fcount,gcount,includeConjugate,rankTol];
  If[!ListQ[rr]||rr==={},Return[{Missing["NoRank"],0,"no frames"}]]; cc=DeleteCases[ Table[
  wd=Quiet[Check[M3WignerCoordinates[finite,iso,r0["V"],rankTol],$Failed],{Power::indet,Infinity::indet}];
  If[AssociationQ[wd],wd["Coordinates"],
If[mGT<=4, ex=Quiet[Check[M3ExactOrbitCoordinates[finite,iso,r0["V"],rankTol],$Failed],{Power::indet,Infinity::indet}];
  If[AssociationQ[ex],ex["Coordinates"],Nothing], Nothing ] ], {r0,rr}],Nothing];
  If[cc==={},Return[{Missing["NoRank"],Length[rr],"no coordinates"}]]; mat=M3CoordinateMatrix[cc];
  If[!MatrixQ[mat],Return[{Missing["NoRank"],Length[rr],"invalid C"}]]; gd=M3GramDiagnostics[mat,rankTol];
  If[!AssociationQ[gd],Return[{Missing["NoRank"],Length[rr],"invalid Gram"}]];
  {gd["RankC"],Length[cc],If[gd["RankC"]==mGT,"complete","incomplete"]} ];

localFrames=Min[maxFrames,Max[4,2 mGT]]; localGrid=Min[gridN,Max[10,8+2 Min[mGT,4]]]; {rank,used,status}=try[localFrames,localGrid];

If[(!IntegerQ[rank]||rank<mGT)&&(localFrames<maxFrames||localGrid<gridN), {rank,used,status}=try[maxFrames,gridN] ];

<|"Rank"->rank,"FramesUsed"->used,"Status"->status|>
];

M3WeightMapCacheKey[lamN_List,maxFrames_,gridN_,includeConjugate_,rankTol_,computeFrameRanks_]:=
  ToString[InputForm[{lamN,maxFrames,gridN,includeConjugate,rankTol,computeFrameRanks}]];

M3WeightMultiplicityData[lamN_List,maxFrames_Integer?Positive,gridN_Integer?Positive,includeConjugate_,rankTol_?NumericQ,computeFrameRanks_]:=Module[
  {key,cached,baseRows,Nrep,rows,finite,iso,rankData,rank,maxMult},
  key=M3WeightMapCacheKey[lamN,maxFrames,gridN,includeConjugate,rankTol,computeFrameRanks];
  cached=Lookup[$M3WeightMapCache,key,Missing["NotCached"]]; If[AssociationQ[cached],Return[cached]];

baseRows=M3AllWeightMultiplicitiesGT[lamN]; If[FailureQ[baseRows],Return[baseRows]]; Nrep=Total[lamN];

rows=Table[ finite=M3FiniteDataFromIntegerData[lamN,row0["n"]]; If[FailureQ[finite],
  Join[row0,<|"W"->{0,0,0},"Coords"->{0.,0.},"MultiplicityIsospin"->Missing["NoIso"],
  "RankC"->Missing["NoRank"],"FullRankInSample"->Missing["NoRank"],"FramesUsed"->0,"Status"->"finite failure"|>],
  iso=M3IsospinData[finite]; rankData=If[TrueQ[computeFrameRanks],
  M3WeightRankSample[lamN,row0["n"],row0["MultiplicityGT"],maxFrames,gridN,includeConjugate,rankTol],
  <|"Rank"->Missing["NotComputed"],"FramesUsed"->0,"Status"->"not computed"|>]; rank=rankData["Rank"];
  Join[row0,<|"W"->finite["W"],"Coords"->M3CartanCoordinates3[finite["W"]],
  "MultiplicityIsospin"->If[AssociationQ[iso],iso["Multiplicity"],Missing["NoIso"]], "RankC"->rank,
  "FullRankInSample"->If[IntegerQ[rank],rank==row0["MultiplicityGT"],Missing["NoRank"]],
  "FramesUsed"->rankData["FramesUsed"],"Status"->rankData["Status"]|>] ], {row0,baseRows}];

maxMult=Max[1,Max[Lookup[rows,"MultiplicityGT",{1}]]];
  cached=<|"Lambda"->lamN,"N"->Nrep,"Rows"->rows,"MaximumMultiplicity"->maxMult, "FrameRanksComputed"->TrueQ[computeFrameRanks],
  "MismatchCount"->Count[Lookup[rows,"FullRankInSample",{}],False], "MissingRankCount"->Count[Lookup[rows,"RankC",{}],_Missing]|>;
  $M3WeightMapCache=Join[$M3WeightMapCache,<|key->cached|>]; cached
];
M3WeightMultiplicityData[lamN_List]:=M3WeightMultiplicityData[lamN,12,24,True,10^-8,False];
M3WeightMultiplicityData[lamN_List,maxFrames_Integer?Positive]:=M3WeightMultiplicityData[lamN,maxFrames,24,True,10^-8,False];
M3WeightMultiplicityData[lamN_List,maxFrames_Integer?Positive,gridN_Integer?Positive]:=M3WeightMultiplicityData[lamN,maxFrames,gridN,True,10^-8,False];
M3WeightMultiplicityData[lamN_List,maxFrames_Integer?Positive,gridN_Integer?Positive,includeConjugate_]:=M3WeightMultiplicityData[lamN,maxFrames,gridN,includeConjugate,10^-8,False];
M3WeightMultiplicityData[lamN_List,maxFrames_Integer?Positive,gridN_Integer?Positive,includeConjugate_,rankTol_?NumericQ]:=M3WeightMultiplicityData[lamN,maxFrames,gridN,includeConjugate,rankTol,False];


(* GELFAND-TSETLIN PATTERN TABLE *)
M3GTTable[gt_Association] := Module[{patterns},
  patterns = gt["Patterns"];
  Grid[
    Prepend[
      Table[{i,patterns[[i]]["Mu1"],patterns[[i]]["Mu2"],patterns[[i]]["Nu1"],patterns[[i]]["j"],patterns[[i]]["Weight"]},
        {i,Length[patterns]}],
      Style[#,Bold]& /@ {"#","mu1","mu2","nu1","j","weight"}],
    Frame -> All,Alignment -> Center,Background -> {1 -> {GrayLevel[.93]}},
    BaseStyle -> {FontFamily -> "Times",FontSize -> 9.5}]
];


M3FailurePanel[result_,title_String:"Calculation error"] := Module[{message,payload},
  payload=If[FailureQ[result]&&Length[result]>=2,result[[2]],Missing["NoPayload"]];
  message=Which[
    FailureQ[result]&&AssociationQ[payload],Lookup[payload,"Message",ToString[result,InputForm]],
    FailureQ[result],ToString[result,InputForm],
    True,ToString[result,InputForm]
  ];
  Panel[
    Column[{
      Style[title,11,Bold,RGBColor[.65,.10,.10],FontFamily->"Times"],
      Style[message,9.5,GrayLevel[.30],FontFamily->"Times"]
    },Spacings->.6],
    Background->RGBColor[1.,.975,.975],FrameMargins->12
  ]
];


M3FrameTable[records_List,independentIndices_List]:=Module[{rows}, rows=Table[ { i, Lookup[records[[i]],"Branch",""],
NumberForm[N[Lookup[records[[i]],"Point",{0,0}]],{5,3}], Lookup[records[[i]],"CoordinateSource","-"],
If[MemberQ[independentIndices,i],Style["independent",Darker[Green],Bold],Style["redundant",Darker[Orange]]] },
{i,Length[records]} ]; Grid[ Prepend[rows,Style[#,Bold]&/@{"frame","branch","Q point","coordinates","rank contribution"}],
Frame->All,Alignment->Left,Background->{1->{GrayLevel[.94]}}, BaseStyle->{FontFamily->"Times",FontSize->9.3} ] ];

(* GRAPHICAL INTERFACE *)
ClearAll[M3WignerFormulaPanel]; M3WignerFormulaPanel[]:=Framed[ Row[{
Style["Orbital amplitude by channel",10.2,Bold,FontFamily->"Times",FontColor->RGBColor[.14,.27,.45]], Spacer[22], Style[
TraditionalForm@HoldForm[ Subscript[A,j][\[Beta]]== Sqrt[(2 j+1) Subscript[\[ScriptCapitalN],W,j]]
Subsuperscript[d,m,\[CapitalDelta]j]^j[\[Beta]] ], 15.0,FontFamily->"Times",FontColor->GrayLevel[.12] ] },Alignment->Center],
Background->RGBColor[.975,.986,.998], FrameStyle->Directive[RGBColor[.72,.80,.89],AbsoluteThickness[.75]],
FrameMargins->{{16,16},{9,9}},RoundingRadius->6 ];

ClearAll[M3WignerSummaryArticle]; M3WignerSummaryArticle[finite_Association,iso_Association]:=Module[{jvals,jtext,dark},
jvals=iso["JValues"]; jtext=If[jvals==={},"-",Row[Riffle[TraditionalForm/@jvals,",  "]]]; dark=RGBColor[.12,.27,.45]; Framed[
Grid[{{ Column[{Style["Projected weight",8.4,GrayLevel[.42],FontFamily->"Times"],
Style[Row[{"W = ",finite["W"]}],11.2,Bold,FontFamily->"Times",FontColor->dark]},Spacings->.10],
Column[{Style["Magnetic projection",8.4,GrayLevel[.42],FontFamily->"Times"],
Style[Row[{"m = ",iso["m"]}],11.2,Bold,FontFamily->"Times",FontColor->dark]},Spacings->.10],
Column[{Style["Shift",8.4,GrayLevel[.42],FontFamily->"Times"],
Style[Row[{"\[CapitalDelta]j = ",iso["DeltaJ"]}],11.2,Bold,FontFamily->"Times",FontColor->dark]},Spacings->.10],
Column[{Style["Allowed channels",8.4,GrayLevel[.42],FontFamily->"Times"],
Style[Row[{Subscript["J","W"]," = {",jtext,"}"}],10.6,Bold,FontFamily->"Times",FontColor->dark]},Spacings->.10],
Column[{Style["Multiplicity",8.4,GrayLevel[.42],FontFamily->"Times"],
Style[Row[{Subscript["m",\[CapitalLambda]],"(W) = ",Length[jvals]}],11.4,Bold,FontFamily->"Times",FontColor->RGBColor[.06,.38,.59]]},Spacings->.10]
}},Alignment->Center,Spacings->{1.45,.2}], Background->RGBColor[.965,.980,.995],
FrameStyle->Directive[RGBColor[.80,.85,.91],AbsoluteThickness[.75]], FrameMargins->{{14,14},{8,8}},RoundingRadius->6 ] ];

ClearAll[M3WignerRotation2]; M3WignerRotation2[beta_]:={{Cos[beta/2],-Sin[beta/2]},{Sin[beta/2],Cos[beta/2]}};

ClearAll[M3WignerAngularResponseFigure]; M3WignerAngularResponseFigure[finite_Association,iso_Association]:=Module[
{m,deltaJ,jvals,beta,exprs,k,colors,legend,maxVal,samples,betaTicks}, m=iso["m"];deltaJ=iso["DeltaJ"];jvals=iso["JValues"];
If[jvals==={},Return[Panel["There are no allowed isospin channels for this weight."]]];

k=Length[jvals]; colors=If[ k==1, {RGBColor[.14,.34,.60]}, Table[ColorData["DarkRainbow"][Rescale[i,{1,k},{.10,.90}]],{i,k}]
];

exprs=Table[ With[{jj=j},Abs[M3WignerDPolynomial[jj,m,deltaJ,M3WignerRotation2[beta]]]^2], {j,jvals} ];

samples=Flatten[Table[N[exprs/.beta->b],{b,Subdivide[0,Pi,80]}]]; maxVal=Quiet@Check[Max[Select[samples,NumericQ]],1.];
If[!NumericQ[maxVal]||maxVal<=0,maxVal=1.];

betaTicks={ {0,Style["0",9.3,FontFamily->"Times"]}, {Pi/4,Style[TraditionalForm[\[Pi]/4],9.3,FontFamily->"Times"]},
{Pi/2,Style[TraditionalForm[\[Pi]/2],9.3,FontFamily->"Times"]},
{3 Pi/4,Style[TraditionalForm[3 \[Pi]/4],9.3,FontFamily->"Times"]},
{Pi,Style[TraditionalForm[\[Pi]],9.3,FontFamily->"Times"]} };

legend=LineLegend[ colors, (Style[Row[{"j = ",#}],8.8,FontFamily->"Times"]&/@jvals), LegendLayout->"Column",
LegendMarkerSize->14, LabelStyle->Directive[FontFamily->"Times",FontSize->8.8], LegendFunction->Identity, Spacings->.25 ];

M3ScientificResultPanel[ "Wigner-D angular response",
"Angular intensity |d^j_{m,\[CapitalDelta]j}(\[Beta])|^2 for each allowed isospin channel.", Plot[
Evaluate[exprs],{beta,0,Pi}, PlotStyle->(Directive[#,AbsoluteThickness[2.35]]&/@colors), PlotRange->{0,1.06 maxVal},
PlotRangePadding->{{Scaled[.01],Scaled[.01]},{0,Scaled[.035]}}, Frame->True,Axes->False, FrameLabel->{ { Style[
TraditionalForm[ Superscript[ Row[{"|",Subsuperscript["d",Row[{"m",",","\[CapitalDelta]j"}],"j"],"(\[Beta])","|"}], 2 ] ],
12,FontFamily->"Times" ], None }, { Style[TraditionalForm[\[Beta]],13,Italic,FontFamily->"Times"], None } },
FrameTicks->{{Automatic,None},{betaTicks,None}}, FrameStyle->Directive[GrayLevel[.16],AbsoluteThickness[1.05]],
FrameTicksStyle->Directive[GrayLevel[.16],9.2,FontFamily->"Times"], GridLines->None,
Epilog->{Inset[legend,Scaled[{.965,.955}],{Right,Top}]}, Background->White, ImagePadding->{{58,18},{52,20}}, ImageSize->535,
PlotPoints->85, MaxRecursion->1, PerformanceGoal->"Quality", BaseStyle->{FontFamily->"Times",FontSize->10.2} ] ] ];

(* Four representative orbital-coordinate distributions built directly from records. *)
ClearAll[M3OrbitalCoordinatesFramesData]; M3OrbitalCoordinatesFramesData[records_List,iso_Association,tol_:10^-9]:=Module[
{jvals,valid,plusInterior,representativePlus,nearCandidates,interior,nearBoundary,
plusRecords,minusRecords,pairCandidates,pair,plusRec,minusRec, probability,pInterior,pNear,pPlus,pMinus,
branchOf,pointDistance,qDistance,sameQ,pairScore},

jvals=Lookup[iso,"JValues",{}]; If[!ListQ[jvals]||jvals==={},
Return[Failure["NoJValues",<|"Message"->"No allowed isospin channels are available."|>]] ];

branchOf[rec_]:=ToString[Lookup[rec,"Branch",""]];

valid=Select[ records, AssociationQ[#]&&
KeyExistsQ[#,"Coordinates"]&&VectorQ[N[# ["Coordinates"]],NumericQ]&&Norm[N[# ["Coordinates"]]]>tol&&
Length[# ["Coordinates"]]===Length[jvals]&& KeyExistsQ[#,"Point"]&&ListQ[# ["Point"]]&&Length[# ["Point"]]===2&&
KeyExistsQ[#,"Q"]&&MatrixQ[N[# ["Q"]]]&& KeyExistsQ[#,"Heron"]&&NumericQ[N[# ["Heron"]]]&& MemberQ[{"+","-"},branchOf[#]]& ];

If[valid==={},
Return[Failure["NoValidOrbitalRecords",<|"Message"->"No stored frame contains valid coordinates, Point, Q, Heron, and branch information."|>]]
];

plusRecords=Select[valid,branchOf[#]=="+"&]; minusRecords=Select[valid,branchOf[#]=="-"&];

If[plusRecords==={}||minusRecords==={},
Return[Failure["NoConjugateBranchPair",<|"Message"->"The stored records do not contain both conjugate branches with valid orbital coordinates."|>]]
];

pointDistance[rp_,rm_]:=Norm[N[rp["Point"]]-N[rm["Point"]]]; qDistance[rp_,rm_]:=Max[Abs[Flatten[N[rp["Q"]]-N[rm["Q"]]]]];
sameQ[r1_,r2_]:=TrueQ[pointDistance[r1,r2]<=1000 tol||qDistance[r1,r2]<=1000 tol];

(* First identify one stored conjugate pair at exactly the same Q. *)
pairCandidates=Flatten[ Table[ If[sameQ[rp,rm],{{rp,rm}},Nothing], {rp,plusRecords},{rm,minusRecords} ], 2 ];

If[pairCandidates==={},
Return[Failure["NoConjugateBranchPair",<|"Message"->"No stored Branch + / Branch - pair was found at the same bistochastic point Q."|>]]
];

pairScore[pr_]:=With[ {rp=pr[[1]],rm=pr[[2]]}, pointDistance[rp,rm]+qDistance[rp,rm]+10^-6 Abs[N[rp["Heron"]]-N[rm["Heron"]]]
]; pair=First@MinimalBy[pairCandidates,pairScore]; plusRec=pair[[1]]; minusRec=pair[[2]];

(* Use Branch + for the two representative non-pair cases. *)
plusInterior=Select[ valid, branchOf[#]=="+"&&N[# ["Heron"]]>10 tol& ]; If[plusInterior==={},
Return[Failure["NoPositiveHeronPlusFrames",<|"Message"->"No admissible Branch + frame with strictly positive Heron value is available."|>]]
];

(* The generic interior point must be a different Q from the conjugate pair. *)
representativePlus=Select[plusInterior,!sameQ[#,plusRec]&]; If[representativePlus==={},
Return[Failure["NoDistinctInteriorFrame",<|"Message"->"No interior Branch + frame at a Q distinct from the selected conjugate pair is available."|>]]
]; interior=First@MaximalBy[representativePlus,N[# ["Heron"]]&];

(* The near-boundary case is also chosen at a third, distinct Q whenever possible. *)
nearCandidates=Select[representativePlus,!sameQ[#,interior]&]; If[nearCandidates==={},
Return[Failure["NoDistinctNearBoundaryFrame",<|"Message"->"A third distinct positive-Heron Branch + point is required for the near-boundary panel."|>]]
]; nearBoundary=First@MinimalBy[nearCandidates,N[# ["Heron"]]&];

probability[rec_]:=Module[{c,p}, c=N[rec["Coordinates"]/Norm[rec["Coordinates"]]]; p=Chop[Abs[c]^2,tol];
If[Total[p]>0,p=N[p/Total[p]]]; p ];

pInterior=probability[interior]; pNear=probability[nearBoundary]; pPlus=probability[plusRec]; pMinus=probability[minusRec];

If[!And@@(Length[#]===Length[jvals]&/@{pInterior,pNear,pPlus,pMinus}),
Return[Failure["OrbitalDimensionMismatch",<|"Message"->"The stored coordinate dimension does not match iso[\"JValues\"]."|>]]
];

<| "JValues"->jvals, "InteriorRecord"->interior, "NearBoundaryRecord"->nearBoundary, "PlusRecord"->plusRec,
"MinusRecord"->minusRec, "InteriorProbability"->pInterior, "NearBoundaryProbability"->pNear, "PlusProbability"->pPlus,
"MinusProbability"->pMinus |> ];

ClearAll[M3OrbitalCoordinatesFramesFigure];
M3OrbitalCoordinatesFramesFigure[records_List,iso_Association,tol_:10^-9]:=Module[
{data,jvals,recs,probs,headers,labels,sameFlags,colors,yMax,jTicks, pointText,heronText,metadata,panel,panels},

data=M3OrbitalCoordinatesFramesData[records,iso,tol];
If[FailureQ[data],Return[M3FailurePanel[data,"The representative orbital-coordinate figure could not be constructed"]]];

jvals=data["JValues"]; recs={data["InteriorRecord"],data["NearBoundaryRecord"],data["PlusRecord"],data["MinusRecord"]};
probs={data["InteriorProbability"],data["NearBoundaryProbability"],data["PlusProbability"],data["MinusProbability"]};
headers={"interior","near boundary",TraditionalForm[Row[{"branch ",Subscript["V","+"]}]],TraditionalForm[Row[{"branch ",Subscript["V","-"]}]]};
labels={"(a)","(b)","(c)","(d)"}; sameFlags={False,False,True,True};
colors={RGBColor[.14,.34,.60],RGBColor[.22,.48,.38],RGBColor[.14,.34,.60],RGBColor[.86,.43,.10]};
yMax=Max[.12,Min[1.05,1.10 Max[Flatten[probs]]]];

jTicks=Table[ {i,Style[TraditionalForm[jvals[[i]]],9.2,FontFamily->"Times"]}, {i,Length[jvals]} ];

pointText[rec_,same_]:=Row[{ If[TrueQ[same],"same Q = ","Q = "],
"(",NumberForm[N[rec["Point"][[1]]],{4,3}],", ",NumberForm[N[rec["Point"][[2]]],{4,3}],")" }];
heronText[rec_]:=Row[{TraditionalForm[Subscript["H","Heron"]]," = ",NumberForm[N[rec["Heron"]],4]}];
metadata[rec_,same_]:=Style[Row[{pointText[rec,same],Spacer[12],heronText[rec]}],8.8,GrayLevel[.28],FontFamily->"Times"];

panel[i_]:=Column[ {
Row[{Style[labels[[i]],10.2,Bold,FontFamily->"Times"],Spacer[7],Style[headers[[i]],10.2,Bold,FontFamily->"Times"]}],
metadata[recs[[i]],sameFlags[[i]]], ListPlot[ Transpose[{Range[Length[jvals]],probs[[i]]}], Joined->False,
PlotMarkers->{Graphics[{EdgeForm[None],Disk[]}],7}, PlotStyle->Directive[colors[[i]],PointSize[.012]], Filling->Axis,
FillingStyle->Directive[colors[[i]],AbsoluteThickness[1.45]], PlotRange->{{.55,Length[jvals]+.45},{0,yMax}},
PlotRangePadding->{{0,0},{0,Scaled[.025]}}, Frame->True,Axes->False, FrameTicks->{{Automatic,None},{jTicks,None}},
FrameLabel->{ Style[TraditionalForm[j],10.5,FontFamily->"Times"],
Style[TraditionalForm[Superscript[Row[{"|",Subscript["c","j"],"(V)|"}],2]],10.5,FontFamily->"Times"] },
FrameStyle->Directive[GrayLevel[.16],AbsoluteThickness[.95]],
FrameTicksStyle->Directive[GrayLevel[.16],8.8,FontFamily->"Times"], Background->White, ImagePadding->{{48,12},{42,12}},
ImageSize->360, BaseStyle->{FontFamily->"Times",FontSize->10.2} ] }, Alignment->Center, Spacings->{.18,.18} ];

panels=panel/@Range[4]; Grid[ {{panels[[1]],panels[[2]]},{panels[[3]],panels[[4]]}}, Alignment->{Center,Top},
Spacings->{.65,.60}, Background->White ] ];

ClearAll[M3WignerOrbitalProbabilityHeatmap]; M3WignerOrbitalProbabilityHeatmap[finite_Association,iso_Association]:=Module[
{p,q,W,m,deltaJ,jvals,k,nBeta,betaGrid,nWeights,probMatrix,raw,den,cf,xTicks,yTicks,heat,legend,separators},
p=finite["p"];q=finite["q"];W=finite["W"];m=iso["m"];deltaJ=iso["DeltaJ"];jvals=iso["JValues"];
If[jvals==={},Return[Panel["There are no allowed isospin channels for this weight."]]];
k=Length[jvals];nBeta=180;betaGrid=Subdivide[0.,N[Pi],nBeta]; nWeights=N[M3NWeight[p,q,W,#]&/@jvals];
probMatrix=Transpose@Table[
raw=Table[Sqrt[N[(2 jvals[[a]]+1) nWeights[[a]]]] M3WignerDPolynomial[jvals[[a]],m,deltaJ,M3WignerRotation2[b]],{a,k}];
den=Total[Abs[N[raw]]^2]; If[NumericQ[den]&&den>10^-14,Clip[N[Abs[raw]^2/den],{0,1}],ConstantArray[0.,k]], {b,betaGrid} ];
cf=Function[z,Blend[{RGBColor[.985,.990,.997],RGBColor[.84,.91,.96],RGBColor[.50,.72,.85],RGBColor[.20,.48,.70],RGBColor[.08,.28,.52],RGBColor[.94,.58,.12]},Clip[z,{0,1}]]];
xTicks=Table[{1+nBeta frac,Switch[frac,0,Style["0",9.2,FontFamily->"Times"],1/4,Style[TraditionalForm[\[Pi]/4],9.2,FontFamily->"Times"],
1/2,Style[TraditionalForm[\[Pi]/2],9.2,FontFamily->"Times"],3/4,Style[TraditionalForm[3 \[Pi]/4],9.2,FontFamily->"Times"],
1,Style[TraditionalForm[\[Pi]],9.2,FontFamily->"Times"]]},{frac,{0,1/4,1/2,3/4,1}}];
yTicks=Table[{a,Style[ToString[jvals[[k-a+1]]],9.5,Bold,FontFamily->"Times"]},{a,k}];
separators=If[k<=1,{},Table[{Directive[White,Opacity[.72],AbsoluteThickness[.75]],Line[{{.5,a+.5},{nBeta+1.5,a+.5}}]},{a,1,k-1}]];
legend=BarLegend[{cf,{0,1}},LegendLabel->Style[TraditionalForm[Subscript[P,j][\[Beta]]],9.6,Bold,FontFamily->"Times"],
LabelStyle->Directive[FontFamily->"Times",FontSize->9.0],LegendMarkerSize->{14,185}]; heat=ArrayPlot[
Reverse[probMatrix],ColorFunctionScaling->False,ColorFunction->cf,PlotRange->{0,1},AspectRatio->.72,
Frame->True,Axes->False,FrameTicks->{{yTicks,None},{xTicks,None}},
FrameLabel->{{Style["isospin channel  j",11.4,Bold,FontFamily->"Times"],None},{Style[TraditionalForm[\[Beta]],12.4,Bold,FontFamily->"Times"],None}},
FrameStyle->Directive[GrayLevel[.16],AbsoluteThickness[1.05]],FrameTicksStyle->Directive[GrayLevel[.16],9.1,FontFamily->"Times"],
Epilog->separators,PlotLegends->Placed[legend,Right],Background->White,
ImagePadding->{{60,62},{50,20}},ImageSize->535,BaseStyle->{FontFamily->"Times",FontSize->10.2} ]; M3ScientificResultPanel[
"Orbital probability by isospin channel",
"Map of P_j(\[Beta]): the most intense regions show which channel dominates at each angle.", heat ] ];


M3SelectedPhysicsData[lam_List,n_List,finite_Association,iso_Association,rec_Association,tol_:10^-9] := Module[
  {V,coords,rg,sigmaData,sigma,sigmaDiag,pOrb,pPhys},
  V=rec["V"];
  coords=Lookup[rec,"Coordinates",Missing["NoCoordinates"]];
  If[!VectorQ[coords,NumericQ]||Norm[N[coords]]<=tol,
    Return[Failure["InvalidCoordinates",<|
      "Message"->"The selected frame does not contain valid numerical orbital coordinates."
    |>]]
  ];
  coords=N[coords/Norm[coords]];

  rg=M3RhoGammaFromFrame[lam,n,V,tol];
  If[FailureQ[rg]||!AssociationQ[rg],Return[rg]];
  sigmaData=M3SigmaFromMatrix[finite,iso,rg["Gamma"],tol];
  If[FailureQ[sigmaData]||!AssociationQ[sigmaData],Return[sigmaData]];
  sigma=sigmaData["Sigma"];
  sigmaDiag=M3SigmaDiagnostics[sigma,tol];
  If[FailureQ[sigmaDiag]||!AssociationQ[sigmaDiag],Return[sigmaDiag]];

  pOrb=Chop[Abs[coords]^2,tol];
  If[Total[pOrb]>0,pOrb=N[pOrb/Total[pOrb]]];
  pPhys=Chop[Re[sigmaDiag["PhysicalJProbabilities"]],tol];
  If[Total[pPhys]>0,pPhys=N[pPhys/Total[pPhys]]];

  <|"RhoGamma"->rg,"Sigma"->sigma,"SigmaDiagnostics"->sigmaDiag,
    "POrb"->pOrb,"PPhys"->pPhys|>
];

(* Reference panels \[LongDash] graphical counterparts of SU3Interactive and MultiphotonProbabilities, built entirely with the M3 backend. *)

If[!AssociationQ[$M3ReferenceTorusCache],$M3ReferenceTorusCache=<||>];
If[!AssociationQ[$M3ReferenceSectorCache],$M3ReferenceSectorCache=<||>];
If[!AssociationQ[$M3ReferenceGlobalDetectionCache],$M3ReferenceGlobalDetectionCache=<||>];

(* Deterministic sampling of the relative torus T^2. The global phase is fixed to zero because it cancels in t_phi rho t_phi^\[Dagger]. *)
M3ReferencePhasePairs[n_Integer?Positive] := Module[{a,b},
  a=(Sqrt[5]-1)/2; b=Sqrt[2]-1;
  Table[2 Pi {FractionalPart[k a],FractionalPart[k b]},{k,1,n}]
];

M3ReferenceOrderedPolygon[pts_List] := Module[{u,c},
  u=DeleteDuplicates[Chop[N[pts],10^-12]];
  If[Length[u]<=2,Return[u]];
  c=Mean[u];
  SortBy[u,Arg[(#[[1]]-c[[1]])+I (#[[2]]-c[[2]])]&]
];

(* Light barycentric grid of the normalized simplex n1+n2+n3=1. *)
M3ReferenceSimplexGrid[Nrep_Integer?Positive] := Module[{step,idx,ks},
  step=Max[1,Ceiling[Nrep/15]];
  idx=DeleteDuplicates[Join[Range[0,Nrep,step],{Nrep}]];
  ks=N[idx/Nrep];
  {
    Directive[GrayLevel[.88],AbsoluteThickness[.36],Opacity[.80]],
    Table[Line[M3CartanCoordinates3 /@ {{c,0,1-c},{c,1-c,0}}],{c,ks}],
    Table[Line[M3CartanCoordinates3 /@ {{0,c,1-c},{1-c,c,0}}],{c,ks}],
    Table[Line[M3CartanCoordinates3 /@ {{0,1-c,c},{1-c,0,c}}],{c,ks}]
  }
];

M3ReferenceMatrixGrid[A_?MatrixQ] := Module[{m},
  m=Map[If[NumericQ[#],If[Abs[N[#]]<10^-10,0,NumberForm[N[#],{5,3}]],#]&,
    Chop[N[A],10^-10],{2}];
  Grid[m,Frame->All,Alignment->Center,
    Dividers->Directive[GrayLevel[.84],AbsoluteThickness[.55]],
    Background->White,ItemSize->All,BaseStyle->{FontFamily->"Times",FontSize->8.7}]
];

(* Data for the panel equivalent to MultiphotonParameterSetter: nu(phi)=diag[U t_phi rho t_phi^\[Dagger] U^\[Dagger]]. *)
M3ReferenceTorusMapData[
  lam_List,n_List,finite_Association,rec_Association,U_?MatrixQ,
  sampleCount_Integer?Positive,tol_:10^-9] := Module[
  {key,cached,rg,rho,gamma,a,target,parts,typical,phases,mapPoints,mapCoords,
   t,nu,Nrep,rows,weightPoints,weightCoords,weylPts,weylCoords},

  key=Hash[{finite["LambdaN"],finite["nN"],Chop[N[rec["V"]],tol],Chop[N[U],tol],sampleCount}];
  cached=Lookup[$M3ReferenceTorusCache,key,Missing["NotCached"]];
  If[AssociationQ[cached],Return[cached]];

  rg=M3RhoGammaFromFrame[lam,n,rec["V"],tol];
  If[FailureQ[rg]||!AssociationQ[rg],Return[rg]];
  rho=rg["Rho"]; gamma=rg["Gamma"];
  a=Reverse@Sort[Chop[Re@Eigenvalues[N[rho]],tol]];
  Nrep=finite["N"];
  target=Nrep a;
  parts=M3Partitions3[Nrep];
  typical=First@MinimalBy[parts,Norm[N[#]-target]&];

  phases=M3ReferencePhasePairs[sampleCount];
  mapPoints=Table[
    t=DiagonalMatrix[{1,Exp[I ph[[1]]],Exp[I ph[[2]]]}];
    nu=Chop[Re@Diagonal[N[U] . t . N[rho] . ConjugateTranspose[t] . ConjugateTranspose[N[U]]],tol];
    If[Total[nu]!=0,N[nu/Total[nu]],N[nu]],
    {ph,phases}
  ];
  mapCoords=M3CartanCoordinates3 /@ mapPoints;

  rows=M3AllWeightMultiplicitiesGT[typical];
  If[FailureQ[rows],Return[rows]];
  weightPoints=N[Lookup[rows,"n",{}]/Nrep];
  weightCoords=M3CartanCoordinates3 /@ weightPoints;

  weylPts=DeleteDuplicates[N[Permutations[typical]/Nrep]];
  weylCoords=M3ReferenceOrderedPolygon[M3CartanCoordinates3 /@ weylPts];

  cached=<|
    "Rho"->rho,"Gamma"->gamma,"AccessibilitySpectrum"->a,
    "TypicalPartition"->typical,"TargetPartition"->target,
    "SelectedWeight"->finite["nN"],"N"->Nrep,
    "MapPoints"->mapPoints,"MapCoordinates"->mapCoords,
    "WeightPoints"->weightPoints,"WeightCoordinates"->weightCoords,
    "WeylCoordinates"->weylCoords,"SampleCount"->Length[mapPoints],
    "U"->U,"LambdaBar"->lam,"NBar"->n
  |>;
  $M3ReferenceTorusCache=Join[$M3ReferenceTorusCache,<|key->cached|>];
  cached
];


M3ReferenceTorusDashboard[data_Association]:=Module[ {
Nrep,gamma,Umat,a,typical,nN,mapCoords,weightCoords,weylCoords,simplex,
selectedCoord,mapGraphic,mapGraphicArticle,legend,legendArticle,paramTable,specErr,rankGamma,dark,teal,mag,
gridGray,pointBlue,torusBlack,articleXr,articleYr,articleAspect, gammaCard,uCard,coordCard,equationCard,headerCard },

Nrep=data["N"]; gamma=data["Gamma"]; Umat=data["U"]; a=data["AccessibilitySpectrum"]; typical=data["TypicalPartition"];
nN=data["SelectedWeight"]; mapCoords=data["MapCoordinates"]; weightCoords=data["WeightCoordinates"];
weylCoords=data["WeylCoordinates"];

simplex=M3ReferenceOrderedPolygon[ M3CartanCoordinates3/@{{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}} ];

selectedCoord=M3CartanCoordinates3[N[nN/Nrep]]; specErr=Norm[Sort[N[a]]-Sort[N[data["LambdaBar"]]]];
rankGamma=MatrixRank[N[gamma],Tolerance->10^-9];

dark=RGBColor[.10,.24,.41]; teal=RGBColor[.08,.56,.76]; mag=RGBColor[.86,.12,.10]; gridGray=GrayLevel[.86];
pointBlue=RGBColor[.18,.47,.74]; torusBlack=GrayLevel[.06];

gammaCard=Framed[ Column[ { Style[ TraditionalForm@HoldForm[\[CapitalGamma]], 11.5,Bold,FontFamily->"Times",FontColor->dark
], M3ReferenceMatrixGrid[gamma] }, Alignment->Center, Spacings->.28 ], Background->RGBColor[.988,.993,.998],
FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.65]], FrameMargins->{{9,9},{8,8}}, RoundingRadius->6 ];

uCard=Framed[ Column[ { Style[ TraditionalForm@HoldForm[\[ScriptCapitalU]], 11.5,Bold,FontFamily->"Times",FontColor->dark ],
M3ReferenceMatrixGrid[Umat] }, Alignment->Center, Spacings->.28 ], Background->RGBColor[.988,.993,.998],
FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.65]], FrameMargins->{{9,9},{8,8}}, RoundingRadius->6 ];

paramTable=Grid[ { { Style["N",Bold],Nrep, Style["rank \[CapitalGamma]",Bold],rankGamma }, {
Style["accessibility spectrum  a",Bold], NumberForm[N[a],{5,3}], SpanFromLeft,SpanFromLeft }, {
Style["typical partition",Bold],typical, Style["selected weight",Bold],nN } }, Frame->All, Alignment->Left,
Spacings->{.8,.43}, Background->{1->{RGBColor[.97,.982,.995]}}, Dividers->Directive[GrayLevel[.82],AbsoluteThickness[.5]],
BaseStyle->{FontFamily->"Times",FontSize->9.35} ];

headerCard=M3ScientificResultPanel[ "Toroidal-map parameters",
"Fix \[Rho], \[CapitalGamma], and \[ScriptCapitalU], and scan the relative phase torus.", Column[ { Grid[
{{gammaCard,uCard}}, Alignment->Top, Spacings->{1.2,0} ], paramTable }, Alignment->Center, Spacings->.55 ] ];

mapGraphic=Graphics[ { { FaceForm[RGBColor[.975,.978,.985]], EdgeForm[Directive[GrayLevel[.76],AbsoluteThickness[.95]]],
Polygon[simplex] },

M3ReferenceSimplexGrid[Nrep],

If[ Length[weylCoords]>=3, { FaceForm[Directive[teal,Opacity[.26]]],
EdgeForm[Directive[Darker[teal,.32],AbsoluteThickness[1.6]]], Polygon[weylCoords] }, {} ],

{ Directive[pointBlue,PointSize[.0048],Opacity[.44]], Point[weightCoords] },

{ Directive[torusBlack,PointSize[.0029],Opacity[.72]], Point[mapCoords] },

{ Directive[mag,PointSize[.018],Opacity[.95]], Point[selectedCoord] } },

Frame->True, Axes->False, AspectRatio->1,

PlotRange->{{-1.08,1.08},{-1.22,.68}}, PlotRangePadding->Scaled[.012],

FrameLabel->{ Style[ TraditionalForm@HoldForm[Subscript[X,C]], 11.3,FontFamily->"Times" ], Style[
TraditionalForm@HoldForm[Subscript[Y,C]], 11.3,FontFamily->"Times" ] },

FrameStyle->Directive[GrayLevel[.20],AbsoluteThickness[.95]],
FrameTicksStyle->Directive[GrayLevel[.22],9.0,FontFamily->"Times"],

Background->White, ImageSize->500, ImagePadding->{{54,18},{48,18}}, BaseStyle->{FontFamily->"Times",FontSize->10.2} ];

articleXr=MinMax[simplex[[All,1]]]; articleYr=MinMax[simplex[[All,2]]];
articleXr=articleXr+{-0.055,0.055} (articleXr[[2]]-articleXr[[1]]);
articleYr=articleYr+{-0.055,0.055} (articleYr[[2]]-articleYr[[1]]);
articleAspect=(articleYr[[2]]-articleYr[[1]])/(articleXr[[2]]-articleXr[[1]]);

mapGraphicArticle=Graphics[ { { FaceForm[Directive[RGBColor[.985,.986,.992],Opacity[.82]]],
EdgeForm[Directive[GrayLevel[.50],AbsoluteThickness[.72],Opacity[.72]]], Polygon[simplex] },
{ Directive[GrayLevel[.80],AbsoluteThickness[.28],Opacity[.30]], M3ReferenceSimplexGrid[Nrep][[2;;]] },

If[ Length[weylCoords]>=3, { FaceForm[Directive[teal,Opacity[.15]]],
EdgeForm[Directive[Darker[teal,.28],AbsoluteThickness[1.05],Opacity[.82]]], Polygon[weylCoords] }, {} ],

{ Directive[pointBlue,PointSize[.0045],Opacity[.38]], Point[weightCoords] },

{ Directive[torusBlack,PointSize[.0030],Opacity[.73]], Point[mapCoords] },

},

Frame->False, Axes->False, PlotRange->{articleXr,articleYr}, PlotRangePadding->0, PlotRangeClipping->False,
AspectRatio->articleAspect, Background->White, ImageSize->500, ImagePadding->6, ImageMargins->0,
BaseStyle->{FontFamily->"Times",FontSize->10.2} ];

legendArticle=Grid[ {{ Style["\[FilledSquare]",11.5,torusBlack], Style["torus image",8.8,FontFamily->"Times"], Spacer[10],
Style["\[FilledSquare]",11.5,pointBlue], Style["output weights",8.8,FontFamily->"Times"] }}, Alignment->Center,
Spacings->{.28,.10}, Background->None, Frame->None ];

coordCard=Framed[ Row[ { TraditionalForm@ HoldForm[ Subscript[X,C]==Subscript[n,1]-Subscript[n,2] ],

Spacer[26],

TraditionalForm@ HoldForm[ Subscript[Y,C]== (Subscript[n,1]+Subscript[n,2]-2 Subscript[n,3])/Sqrt[3] ] }, Alignment->Center
], Background->RGBColor[.982,.988,.998], FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.6]],
FrameMargins->{{13,13},{7,7}}, RoundingRadius->5, BaseStyle->{FontFamily->"Times",FontSize->11.0,FontColor->dark} ];

legend=Grid[ {{ Style["\[FilledSquare]",12,torusBlack], Style["torus image",9.0,FontFamily->"Times"],

Spacer[12],

Style["\[FilledSquare]",12,pointBlue], Style["spectral region / typical irrep",9.0,FontFamily->"Times"],

Spacer[12],

Style["\[FilledCircle]",11.5,mag], Style["input weight",9.0,FontFamily->"Times"] }}, Alignment->Center, Spacings->{.35,.15}
];

equationCard=Framed[ Style[ TraditionalForm@ HoldForm[ Subscript[\[Nu],\[Phi]]== Diagonal[ \[ScriptCapitalU] .
Subscript[t,\[Phi]] . \[Rho] . Superscript[Subscript[t,\[Phi]],\[Dagger]] . Superscript[\[ScriptCapitalU],\[Dagger]] ] ],
12.4, FontFamily->"Times", FontColor->dark ], Background->RGBColor[.978,.986,.998],
FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.62]], FrameMargins->{{14,14},{7,7}}, RoundingRadius->5 ];

Column[ { Style[ "Toroidal map and spectrum", 14.2,Bold,FontFamily->"Times",FontColor->dark ],

Grid[ {{ Column[ { headerCard, equationCard }, Alignment->Center, Spacings->.42 ],

Column[ { M3ScientificResultPanel[ "Toroidal map to output occupations", "", Column[ { mapGraphic, coordCard, legend },
Alignment->Center, Spacings->.34 ] ],

Column[ { Style[ "Article-style toroidal map", 12.2,Bold,FontFamily->"Times",FontColor->dark ], mapGraphicArticle,
legendArticle }, Alignment->Center, Spacings->.24 ] }, Alignment->Center, Spacings->.70 ] }}, Alignment->{Left,Top},
Spacings->{1.05,0}, BaseStyle->{FontFamily->"Times",FontSize->10.2} ] }, Alignment->Center, Spacings->.58 ] ];

M3ReferenceCompatibleSectors[nN_List]:=

Select[M3Partitions3[Total[nN]],TrueQ[M3GTMultiplicity[#,nN]>0]&];

M3ReferenceSectorData[ lamN_List,nN_List,Gamma_?MatrixQ,U_?MatrixQ,tol_:10^-9]:=Module[
{key,cached,sector,prob,finite,iso,sigma,transitions,detection,joint},

key=Hash[{lamN,nN,Chop[N[Gamma],tol],Chop[N[U],tol]}]; cached=Lookup[$M3ReferenceSectorCache,key,Missing["NotCached"]];
If[AssociationQ[cached],Return[cached]];

sector=M3SectorBlockDataCached[lamN,nN,Gamma,tol]; If[FailureQ[sector]||!AssociationQ[sector],Return[sector]];
prob=Lookup[sector,"Probability",0.]; finite=Lookup[sector,"FiniteData",Missing["NoFiniteData"]];
iso=Lookup[sector,"IsospinData",Missing["NoIsospinData"]];

If[!(AssociationQ[finite]&&AssociationQ[iso])||!KeyExistsQ[sector,"SigmaData"], Return[Failure["ZeroOrUnavailableSector",<|
"Message"->"The selected sector does not contain a normalizable physical block for this Gamma.",
"Lambda"->lamN,"Probability"->prob|>]] ];

sigma=sector["SigmaData"]["Sigma"]; transitions=M3AllTransitions[finite,iso,U,tol];
If[FailureQ[transitions],Return[transitions]]; detection=M3PhysicalOutputDistribution[transitions,sigma,tol];
If[FailureQ[detection],Return[detection]];

joint=Map[ Join[#,<|"ConditionalProbability"->#["Probability"],"JointProbability"->N[prob #["Probability"]]|>]&,
detection["Results"] ];

(* Only the information actually used by the interface is kept in cache. The transition blocks, sigma, and the Detection structure can be large and are discarded once the output probabilities have been constructed. *)
cached=<| "Lambda"->lamN,"nIn"->nN,"SectorProbability"->prob, "Multiplicity"->iso["Multiplicity"],"Rows"->joint,
"ConditionalTotal"->detection["TotalProbability"], "JointTotal"->N[prob detection["TotalProbability"]] |>;
$M3ReferenceSectorCache=Join[$M3ReferenceSectorCache,<|key->cached|>]; cached ];

M3ReferenceWeightDiagram[lamN_List,nN_List]:=Module[ { Nrep,rows,pts,mults,maxm,simplex,weyl,selected,highest,nodeSize,prims,
inputColor,highestColor,legend }, Nrep=Total[lamN]; rows=M3AllWeightMultiplicitiesGT[lamN];
If[FailureQ[rows]||rows==={},Return[Panel["The sector diagram could not be constructed."]]];

pts=M3CartanCoordinates3/@N[Lookup[rows,"n"]/Nrep]; mults=Lookup[rows,"MultiplicityGT"]; maxm=Max[1,Max[mults]];

simplex=M3ReferenceOrderedPolygon[M3CartanCoordinates3/@{{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}}];
weyl=M3ReferenceOrderedPolygon[M3CartanCoordinates3/@N[DeleteDuplicates[Permutations[lamN]]/Nrep]];
selected=M3CartanCoordinates3[N[nN/Nrep]]; highest=M3CartanCoordinates3[N[lamN/Nrep]];

inputColor=RGBColor[.82,.11,.10]; highestColor=RGBColor[.93,.46,.10];

nodeSize[m_]:=.0035 + .0038 Sqrt[N[m/maxm]];

prims=MapThread[ { Directive[GrayLevel[.12],PointSize[nodeSize[#2]]], Point[#1] }&, {pts,mults} ];

legend=Framed[ Row[{ Style["\[FilledCircle]",12,GrayLevel[.10]],Style[" irrep weights   ",8.8,FontFamily->"Times"],
Style["\[FilledCircle]",13,inputColor],Style[" input weight n/N   ",8.8,FontFamily->"Times"],
Style["\[FilledDiamond]",12,highestColor],Style[" highest weight \[CapitalLambda]/N",8.8,FontFamily->"Times"] },Spacer[4]],
Background->Directive[White,Opacity[.88]], FrameStyle->Directive[GrayLevel[.82],AbsoluteThickness[.7]], RoundingRadius->5,
FrameMargins->{{6,6},{4,4}} ];

Legended[ Graphics[ {
{FaceForm[RGBColor[.962,.965,.973]],EdgeForm[Directive[GrayLevel[.70],AbsoluteThickness[.9]]],Polygon[simplex]},
M3ReferenceSimplexGrid[Nrep],
{FaceForm[Directive[RGBColor[.28,.80,.82],Opacity[.30]]],EdgeForm[Directive[RGBColor[.06,.50,.56],AbsoluteThickness[1.25]]],Polygon[weyl]},
prims, { highestColor, EdgeForm[Directive[Darker[highestColor,.15],AbsoluteThickness[.6]]], GeometricTransformation[
Polygon[{{0,.028},{.022,0},{0,-.028},{-.022,0}}], TranslationTransform[highest] ] },
{Directive[inputColor],PointSize[.017],Point[selected]} }, Frame->True, Axes->False, AspectRatio->1,
PlotRange->{{-1.12,1.12},{-1.22,.70}}, PlotRangePadding->Scaled[.014],
FrameStyle->Directive[GrayLevel[.20],AbsoluteThickness[1.0]],
FrameTicksStyle->Directive[GrayLevel[.22],9.0,FontFamily->"Times"],
FrameTicks->{{Automatic,Automatic},{Automatic,Automatic}}, FrameLabel->{
Style[TraditionalForm[(Subscript["n",1]-Subscript["n",2])/N],10.8,FontFamily->"Times"],
Style[TraditionalForm[(Subscript["n",1]+Subscript["n",2]-2 Subscript["n",3])/(Sqrt[3] N)],10.8,FontFamily->"Times"] },
Background->White, ImageSize->350, ImagePadding->{{54,18},{42,22}}, BaseStyle->{FontFamily->"Times",FontSize->10.2} ],
Placed[legend,Below] ] ];


(* Appendix G / Figs. 11\[Dash]13 \[LongDash] 3D output-probability representation with a common 2D/3D layout: 
plot on the left and aligned color bar on the right. *)
$M3PaperMainImageSize={350,285};
$M3PaperLegendPaneSize={70,285};
$M3PaperLegendMarkerSize={11,155};
$M3PaperPlotLegendSpacing=.15;
$M3PaperCanvasImageSize={450,300};

ClearAll[M3ReferencePaperCompose];
M3ReferencePaperCompose[main_,legend_] := Graphics[
  {
    Inset[main,{180,150},Center,$M3PaperMainImageSize],
    Inset[Pane[legend,$M3PaperLegendPaneSize,Alignment->Center],
      {405,150},Center,$M3PaperLegendPaneSize]
  },
  PlotRange->{{0,450},{0,300}},PlotRangePadding->0,PlotRangeClipping->True,
  AspectRatio->300/450,ImageSize->$M3PaperCanvasImageSize,
  ImagePadding->0,ImageMargins->0,Background->White
];

ClearAll[M3ReferenceHexagon2D,M3ReferenceHexPrism];

M3ReferenceHexagon2D[{x_?NumericQ,y_?NumericQ},r_?NumericQ] :=
  Table[{x+r Cos[Pi/6+k Pi/3],y+r Sin[Pi/6+k Pi/3]},{k,0,5}];

M3ReferenceHexPrism[{x_?NumericQ,y_?NumericQ},r_?NumericQ,h_?NumericQ,col_] := Module[
  {base,top,sideCol,topCol,edges},
  base=M3ReferenceHexagon2D[{x,y},r];
  top=Append[#,h]& /@ base;
  sideCol=Blend[{col,GrayLevel[.34]},.12];
  topCol=col;
  edges=Directive[GrayLevel[.26],Opacity[.26],AbsoluteThickness[.24]];

  {
    {
      FaceForm[Directive[sideCol,Opacity[.96]]],EdgeForm[edges],
      Table[
        Polygon[{Append[base[[i]],0.],Append[base[[1+Mod[i,6]]],0.],
          top[[1+Mod[i,6]]],top[[i]]}],
        {i,1,6}
      ]
    },
    {
      FaceForm[Directive[topCol,Opacity[1.]]],EdgeForm[edges],Polygon[top]
    }
  }
];


M3ReferenceSectorProbabilityPlot3D[data_Association,mode_String:"Joint"] := Module[
  {rows,Nrep,coords,vals,valid,vmax,zmax,probLabel,probName,subtitle,sectorPoly,selectedCoord,
   dominantIndex,dominantCoord,dominantVal,xr,yr,dx,dy,hexR,distances,colorF,teal,mag,
   basePlane,bars,dominantGraphic,selectedGraphic,xRange,yRange,xAxisY,yAxisX,xTickVals,
   yTickVals,xTickSize,yTickSize,axesGraphic,xLabelPos,yLabelPos,displayXr,displayYr,axisZ,
   plot,legendTicks,barLegend},

  rows=Lookup[data,"Rows",{}];
  If[rows==={},Return[M3ScientificResultPanel["3D distribution over output weights","","There are no output weights for this sector."]]];

  Nrep=Total[data["Lambda"]];
  coords=M3CartanCoordinates3 /@ (N[Lookup[rows,"nOut"]]/Nrep);
  vals=If[mode==="Conditional",Lookup[rows,"ConditionalProbability"],Lookup[rows,"JointProbability"]];
  vals=Map[Function[z,If[NumericQ[z],Max[0.,Chop[Re@N[z],10^-13]],Indeterminate]],vals];

  valid=Select[Range[Min[Length[coords],Length[vals]]],
    Function[k,VectorQ[coords[[k]],NumericQ]&&NumericQ[vals[[k]]]]];
  If[valid==={},Return[M3ScientificResultPanel["3D distribution over output weights","","There are no numerical probabilities to plot."]]];

  rows=rows[[valid]]; coords=coords[[valid]]; vals=vals[[valid]];
  vmax=Max[10^-14,Max[vals]]; zmax=1.08 vmax;

  probLabel=If[mode==="Conditional",
    TraditionalForm@HoldForm[Subscript[P,\[CapitalLambda]][Superscript[n,"\[Prime]"]]],
    TraditionalForm@HoldForm[P[\[CapitalLambda],Superscript[n,"\[Prime]"]]]
  ];
  probName=If[mode==="Conditional","P_Lambda(n')","P(Lambda,n')"];
  subtitle=If[mode==="Conditional","Probability conditioned on the selected sector.",
    "Joint probability by sector over the physical output weights."];

  sectorPoly=M3ReferenceOrderedPolygon[
    M3CartanCoordinates3 /@ N[DeleteDuplicates[Permutations[data["Lambda"]]]/Nrep]];
  selectedCoord=M3CartanCoordinates3[N[data["nIn"]/Nrep]];

  xr=MinMax[coords[[All,1]]]; yr=MinMax[coords[[All,2]]];
  If[Length[sectorPoly]>=3,
    xr={Min[xr[[1]],Min[sectorPoly[[All,1]]]],Max[xr[[2]],Max[sectorPoly[[All,1]]]]};
    yr={Min[yr[[1]],Min[sectorPoly[[All,2]]]],Max[yr[[2]],Max[sectorPoly[[All,2]]]]};
  ];

  dx=Max[.040,.030 Max[10^-6,xr[[2]]-xr[[1]]]];
  dy=Max[.040,.030 Max[10^-6,yr[[2]]-yr[[1]]]];
  xr=xr+{-dx,dx}; yr=yr+{-dy,dy};

  distances=DeleteCases[
    Flatten@Table[If[j>i,Norm[coords[[i]]-coords[[j]]],Nothing],{i,Length[coords]},{j,Length[coords]}],
    x_/;!NumericQ[x]||x<=10^-10
  ];
  hexR=If[distances==={},.038,Clip[.35 Min[distances],{.020,.049}]];

  teal=RGBColor[.09,.55,.65]; mag=RGBColor[.86,.10,.10];

  colorF[v_?NumericQ]:=Blend[
    {RGBColor[.16,.34,.67],RGBColor[.12,.66,.86],RGBColor[.78,.88,.42],
     RGBColor[.97,.62,.16],RGBColor[.82,.12,.09]},
    Clip[v/vmax,{0.,1.}]
  ];

  basePlane=Graphics3D[{
    If[Length[sectorPoly]>=3,
      {FaceForm[Directive[RGBColor[.86,.95,.96],Opacity[.10]]],
       EdgeForm[Directive[GrayLevel[.68],AbsoluteThickness[.65],Opacity[.60]]],
       Polygon[(Append[#,0.]&) /@ sectorPoly]},
      {}
    ]
  },Lighting->"Neutral"];

  bars=Graphics3D[
    MapThread[
      Function[{row,p,v},
        Tooltip[
          M3ReferenceHexPrism[p,hexR,v,colorF[v]],
          Column[{Style["Output weight",Bold,FontFamily->"Times"],
            Row[{"n' = ",row["nOut"]}],Row[{probName," = ",ScientificForm[v,4]}]},Spacings->.15]
        ]
      ],
      {rows,coords,vals}
    ],
    Lighting->{{"Ambient",GrayLevel[.80]},
      {"Directional",GrayLevel[.96],{2.3,-2.8,4.4}},
      {"Directional",GrayLevel[.72],{-2.0,1.6,3.1}}}
  ];

  dominantIndex=First@Ordering[vals,-1];
  dominantCoord=coords[[dominantIndex]]; dominantVal=vals[[dominantIndex]];

  dominantGraphic=Graphics3D[
    {Directive[RGBColor[.83,.12,.08],PointSize[.012]],
     Point[{dominantCoord[[1]],dominantCoord[[2]],dominantVal}]},
    Lighting->"Neutral"
  ];

  selectedGraphic=Graphics3D[
    {Directive[RGBColor[.88,.05,.76],PointSize[.013],Opacity[.92]],
     Point[Append[selectedCoord,.018 vmax]]},
    Lighting->"Neutral"
  ];

  xRange=xr[[2]]-xr[[1]]; yRange=yr[[2]]-yr[[1]];
  xAxisY=yr[[1]]-.055 yRange; yAxisX=xr[[2]]+.055 xRange;
  axisZ=-.010 zmax;

  displayXr={xr[[1]]-.020 xRange,yAxisX+.085 xRange};
  displayYr={xAxisY-.085 yRange,yr[[2]]+.035 yRange};

  xTickVals=Select[{-1.,-.5,0.,.5,1.},xr[[1]]-10^-10<=#<=xr[[2]]+10^-10&];
  yTickVals=Select[{-1.,-.5,0.,.5,1.},yr[[1]]-10^-10<=#<=yr[[2]]+10^-10&];
  xTickSize=.012 yRange; yTickSize=.012 xRange;

  xLabelPos={Mean[xr],xAxisY-.075 yRange,axisZ};
  yLabelPos={yAxisX+.067 xRange,Mean[yr],axisZ};

  axesGraphic=Graphics3D[{
    {Directive[GrayLevel[.34],AbsoluteThickness[.95]],
     Line[{{xr[[1]],xAxisY,axisZ},{xr[[2]],xAxisY,axisZ}}],
     Line[{{yAxisX,yr[[1]],axisZ},{yAxisX,yr[[2]],axisZ}}]},
    Table[
      {Directive[GrayLevel[.38],AbsoluteThickness[.70]],
       Line[{{x,xAxisY-xTickSize,axisZ},{x,xAxisY+xTickSize,axisZ}}],
       Text[Style[NumberForm[x,{2,1}],8.7,GrayLevel[.20],FontFamily->"Times"],
         {x,xAxisY-.045 yRange,axisZ}]},
      {x,xTickVals}
    ],
    Table[
      {Directive[GrayLevel[.38],AbsoluteThickness[.70]],
       Line[{{yAxisX-yTickSize,y,axisZ},{yAxisX+yTickSize,y,axisZ}}],
       Text[Style[NumberForm[y,{2,1}],8.7,GrayLevel[.20],FontFamily->"Times"],
         {yAxisX+.040 xRange,y,axisZ}]},
      {y,yTickVals}
    ],
    Text[Style[TraditionalForm[Subscript[X,C]],11.3,GrayLevel[.10],FontFamily->"Times"],xLabelPos],
    Text[Style[TraditionalForm[Subscript[Y,C]],11.3,GrayLevel[.10],FontFamily->"Times"],yLabelPos]
  },Lighting->"Neutral"];

  plot=Show[
    basePlane,bars,dominantGraphic,selectedGraphic,axesGraphic,
    Boxed->False,Axes->False,
    PlotRange->{displayXr,displayYr,{axisZ-.003 zmax,zmax}},
    PlotRangeClipping->False,PlotRangePadding->Scaled[.005],
    BoxRatios->{1.,.88,.44},ViewPoint->{1.42,-2.10,1.72},
    ViewVertical->{0,0,1},ViewAngle->.28,SphericalRegion->True,
    Background->White,ImageSize->$M3PaperMainImageSize,
    ImagePadding->{{16,16},{14,14}},
    BaseStyle->{FontFamily->"Times",FontSize->10.2}
  ];

  legendTicks=Table[
    {q,If[q==0,Style["0",8.3,FontFamily->"Times"],
      Style[ScientificForm[q,2],8.3,FontFamily->"Times"]]},
    {q,{0,vmax/4,vmax/2,3 vmax/4,vmax}}
  ];

  barLegend=BarLegend[
    {Function[q,colorF[q]],{0,vmax}},
    Ticks->legendTicks,
    LegendLabel->Style[probLabel,9.5,FontFamily->"Times"],
    LabelStyle->Directive[FontFamily->"Times",FontSize->8.3,GrayLevel[.14]],
    LegendMarkerSize->$M3PaperLegendMarkerSize
  ];

  M3ReferencePaperCompose[plot,barLegend]
];


(* Flat 2D companion to the 3D probability plots: outputs are shown on the normalized three-port occupation simplex, with probability 
encoded by a sequential color map. *)

M3ReferencePaperProbabilitySimplex[
  coords_List,vals_List,Nrep_Integer?Positive,selectedCoord_List,
  probLabel_,tooltipRows_List,supportPoly_:Automatic,dominantCoord_:Automatic
] := Module[
  {good,v,vmax,simplex,poly,xr,yr,xSpan,ySpan,plotAspect,distances,hexR,colorScale,
   probColor,gridStep,gridIdx,gridKs,gridPrims,dominant,legendTicks,barLegend,mainGraphic},

  good=Select[Range[Min[Length[coords],Length[vals],Length[tooltipRows]]],
    Function[k,VectorQ[coords[[k]],NumericQ]&&NumericQ[vals[[k]]]]];

  If[good==={},
    Return[Graphics[
      Text[Style["There are no numerical probabilities to plot.",11,FontFamily->"Times"]],
      Background->White
    ]]
  ];

  v=N[vals[[good]]];
  vmax=Max[10^-14,Max[v]];

  simplex=M3ReferenceOrderedPolygon[
    M3CartanCoordinates3 /@ {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}}
  ];

  poly=Which[
    supportPoly===Automatic,simplex,
    ListQ[supportPoly]&&Length[supportPoly]>=3,M3ReferenceOrderedPolygon[N[supportPoly]],
    True,simplex
  ];

  xr=MinMax[simplex[[All,1]]]; yr=MinMax[simplex[[All,2]]];
  xSpan=xr[[2]]-xr[[1]]; ySpan=yr[[2]]-yr[[1]];
  xr=xr+{-0.035 xSpan,0.035 xSpan};
  yr=yr+{-0.035 ySpan,0.035 ySpan};
  plotAspect=(yr[[2]]-yr[[1]])/(xr[[2]]-xr[[1]]);

  distances=DeleteCases[
    Flatten@Table[
      If[j>i,Norm[coords[[good[[i]]]]-coords[[good[[j]]]]],Nothing],
      {i,Length[good]},{j,Length[good]}
    ],
    x_/;!NumericQ[x]||x<=10^-10
  ];
  hexR=If[distances==={},.035,Clip[.42 Min[distances],{.014,.050}]];

  (* Linear probability color scale. *)
  colorScale[t_?NumericQ]:=Blend[
    {RGBColor[.16,.34,.67],RGBColor[.12,.66,.86],RGBColor[.78,.88,.42],
     RGBColor[.97,.62,.16],RGBColor[.82,.12,.09]},
    Clip[t,{0.,1.}]
  ];
  probColor[p_?NumericQ]:=colorScale[Clip[p/vmax,{0.,1.}]];

  (* Very light barycentric grid: enough to show the simplex structure but
     deliberately weaker than the probability markers. *)
  gridStep=Max[1,Ceiling[Nrep/14]];
  gridIdx=DeleteDuplicates[Join[Range[0,Nrep,gridStep],{Nrep}]];
  gridKs=N[gridIdx/Nrep];
  gridPrims={
    Directive[RGBColor[.68,.67,.74],AbsoluteThickness[.32],Opacity[.28]],
    Table[Line[M3CartanCoordinates3 /@ {{c,0,1-c},{c,1-c,0}}],{c,gridKs}],
    Table[Line[M3CartanCoordinates3 /@ {{0,c,1-c},{1-c,c,0}}],{c,gridKs}],
    Table[Line[M3CartanCoordinates3 /@ {{0,1-c,c},{1-c,0,c}}],{c,gridKs}]
  };

  dominant=If[
    dominantCoord===Automatic,
    coords[[good[[First@Ordering[v,-1]]]]],
    dominantCoord
  ];

  mainGraphic=Graphics[
    {
      (* Full normalized occupation simplex. *)
      {
        FaceForm[Directive[RGBColor[.955,.952,.975],Opacity[.72]]],
        EdgeForm[Directive[GrayLevel[.48],AbsoluteThickness[.72],Opacity[.70]]],
        Polygon[simplex]
      },

      gridPrims,

      (* Selected irrep support / permutahedron, when different from the full
         simplex.  This gives the same "region on a triangle" reading as the
         reference-paper panels. *)
      If[
        Length[poly]===Length[simplex]&&
          Norm[Sort[Flatten[N[poly]]]-Sort[Flatten[N[simplex]]]]<10^-10,
        {},
        {
          FaceForm[Directive[RGBColor[.55,.82,.88],Opacity[.10]]],
          EdgeForm[Directive[RGBColor[.20,.45,.58],AbsoluteThickness[.75],Opacity[.62]]],
          Polygon[poly]
        }
      ],

      (* Probability heat map on the discrete weight lattice. *)
      Table[
        With[
          {k=good[[ii]],p=N[vals[[good[[ii]]]]],pt=N[coords[[good[[ii]]]]]},
          Tooltip[
            {
              EdgeForm[Directive[GrayLevel[.24],AbsoluteThickness[.32],Opacity[.56]]],
              FaceForm[Directive[probColor[p],Opacity[.985]]],
              Polygon[M3ReferenceHexagon2D[pt,hexR]]
            },
            Column[
              {
                Style["Output weight",Bold,FontFamily->"Times"],
                Row[{"n' = ",tooltipRows[[k]]}],
                Row[{"P = ",ScientificForm[p,4]}]
              },
              Spacings->.12
            ]
          ]
        ],
        {ii,Length[good]}
      ],

      (* Input occupation and most probable output. *)
      {Directive[RGBColor[.50,.08,.55],PointSize[.012]],Point[selectedCoord]},
      {
        FaceForm[None],
        EdgeForm[Directive[GrayLevel[.08],AbsoluteThickness[1.05],Opacity[.92]]],
        Polygon[M3ReferenceHexagon2D[dominant,1.12 hexR]]
      }
    },

    Frame->False,Axes->False,
    PlotRange->{xr,yr},PlotRangePadding->0,
    AspectRatio->plotAspect,Background->White,
    ImageSize->$M3PaperMainImageSize,
    ImagePadding->{{8,8},{8,8}},
    ImageMargins->0,PlotRangeClipping->False,
    BaseStyle->{FontFamily->"Times",FontSize->10}
  ];

  legendTicks=Table[
    {q,If[q==0,Style["0",8.3,FontFamily->"Times"],
      Style[ScientificForm[q,2],8.3,FontFamily->"Times"]]},
    {q,{0,vmax/4,vmax/2,3 vmax/4,vmax}}
  ];
  barLegend=BarLegend[
    {Function[q,probColor[q]],{0,vmax}},
    Ticks->legendTicks,
    LegendLabel->Style[probLabel,9.5,FontFamily->"Times"],
    LabelStyle->Directive[FontFamily->"Times",FontSize->8.3,GrayLevel[.14]],
    LegendMarkerSize->$M3PaperLegendMarkerSize
  ];

  (* Use the same fixed two-column publication layout as the 3D plots. *)
  M3ReferencePaperCompose[mainGraphic,barLegend]
];


M3ReferenceSectorProbabilityMap[data_Association,mode_String:"Joint"] := Module[
  {rows,Nrep,coords,vals,valid,probLabel,sectorPoly,selectedCoord,dominantCoord,dominantIndex},

  rows=Lookup[data,"Rows",{}];
  If[rows==={},Return[Graphics[Text[Style["There are no output weights for this sector.",11,FontFamily->"Times"]],Background->White]]];

  Nrep=Total[data["Lambda"]];
  coords=M3CartanCoordinates3 /@ (N[Lookup[rows,"nOut"]]/Nrep);
  vals=If[mode==="Conditional",Lookup[rows,"ConditionalProbability"],Lookup[rows,"JointProbability"]];
  vals=Map[Function[z,If[NumericQ[z],Max[0.,Chop[Re@N[z],10^-13]],Indeterminate]],vals];

  valid=Select[Range[Min[Length[coords],Length[vals]]],
    Function[k,VectorQ[coords[[k]],NumericQ]&&NumericQ[vals[[k]]]]];
  If[valid==={},Return[Graphics[Text[Style["There are no numerical probabilities to plot.",11,FontFamily->"Times"]],Background->White]]];

  probLabel=If[mode==="Conditional",
    TraditionalForm@HoldForm[Subscript[P,\[CapitalLambda]][Superscript[n,"\[Prime]"]]],
    TraditionalForm@HoldForm[P[\[CapitalLambda],Superscript[n,"\[Prime]"]]]];

  sectorPoly=M3ReferenceOrderedPolygon[
    M3CartanCoordinates3 /@ N[DeleteDuplicates[Permutations[data["Lambda"]]]/Nrep]];
  selectedCoord=M3CartanCoordinates3[N[data["nIn"]/Nrep]];
  dominantIndex=valid[[First@Ordering[N[vals[[valid]]],-1]]];
  dominantCoord=coords[[dominantIndex]];

  M3ReferencePaperProbabilitySimplex[
    coords,vals,Nrep,selectedCoord,probLabel,Lookup[rows,"nOut",{}],sectorPoly,dominantCoord]
];


(* Appendix G.6.3 / Figs. 11\[Dash]13 \[LongDash] total physical distribution from the exact sum over all sectors. *)

M3ReferenceGlobalDetectionData[
  nN_List,Gamma_?MatrixQ,U_?MatrixQ,tol_:10^-9] := Module[
  {key,cached,data},

  key=Hash[{nN,Chop[N[Gamma],tol],Chop[N[U],tol],tol}];
  cached=Lookup[$M3ReferenceGlobalDetectionCache,key,Missing["NotCached"]];
  If[AssociationQ[cached]||FailureQ[cached],Return[cached]];

  data=M3FullSchurWeylDetection[nN,Gamma,U,tol];
  $M3ReferenceGlobalDetectionCache=
    Join[$M3ReferenceGlobalDetectionCache,<|key->data|>];
  data
];


M3ReferenceGlobalProbabilityPlot3D[data_Association,nN_List] := Module[
  {rows,Nrep,coords,vals,valid,vmax,zmax,selectedCoord,dominantIndex,dominantCoord,dominantVal,
   simplex,xr,yr,dx,dy,hexR,distances,colorF,bars,basePlane,dominantGraphic,selectedGraphic,
   xRange,yRange,xAxisY,yAxisX,xTickVals,yTickVals,xTickSize,yTickSize,axesGraphic,
   xLabelPos,yLabelPos,displayXr,displayYr,axisZ,plot,legendTicks,barLegend},

  rows=Lookup[data,"Results",{}];
  If[rows==={},Return[M3ScientificResultPanel["Total physical distribution 3D","","There are no global probabilities to plot."]]];

  Nrep=Total[nN];
  coords=M3CartanCoordinates3 /@ (N[Lookup[rows,"nOut"]]/Nrep);
  vals=Lookup[rows,"Probability",{}];
  vals=Map[Function[z,If[NumericQ[z],Max[0.,Chop[Re@N[z],10^-13]],Indeterminate]],vals];

  valid=Select[Range[Min[Length[coords],Length[vals]]],
    Function[k,VectorQ[coords[[k]],NumericQ]&&NumericQ[vals[[k]]]]];
  If[valid==={},Return[M3ScientificResultPanel["Total physical distribution 3D","","There are no numerical probabilities to plot."]]];

  rows=rows[[valid]]; coords=coords[[valid]]; vals=vals[[valid]];
  vmax=Max[10^-14,Max[vals]]; zmax=1.08 vmax;

  simplex=M3ReferenceOrderedPolygon[M3CartanCoordinates3 /@ {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}}];
  selectedCoord=M3CartanCoordinates3[N[nN/Nrep]];

  xr=MinMax[Join[coords[[All,1]],simplex[[All,1]]]];
  yr=MinMax[Join[coords[[All,2]],simplex[[All,2]]]];
  dx=Max[.040,.030 Max[10^-6,xr[[2]]-xr[[1]]]];
  dy=Max[.040,.030 Max[10^-6,yr[[2]]-yr[[1]]]];
  xr=xr+{-dx,dx}; yr=yr+{-dy,dy};

  distances=DeleteCases[
    Flatten@Table[If[j>i,Norm[coords[[i]]-coords[[j]]],Nothing],{i,Length[coords]},{j,Length[coords]}],
    x_/;!NumericQ[x]||x<=10^-10
  ];
  hexR=If[distances==={},.038,Clip[.35 Min[distances],{.016,.049}]];

  colorF[v_?NumericQ]:=Blend[
    {RGBColor[.16,.34,.67],RGBColor[.12,.66,.86],RGBColor[.78,.88,.42],
     RGBColor[.97,.62,.16],RGBColor[.82,.12,.09]},
    Clip[v/vmax,{0.,1.}]
  ];

  basePlane=Graphics3D[
    {FaceForm[Directive[GrayLevel[.975],Opacity[.16]]],
     EdgeForm[Directive[GrayLevel[.72],AbsoluteThickness[.65],Opacity[.60]]],
     Polygon[(Append[#,0.]&) /@ simplex]},
    Lighting->"Neutral"
  ];

  bars=Graphics3D[
    MapThread[
      Function[{row,p,v},
        Tooltip[
          M3ReferenceHexPrism[p,hexR,v,colorF[v]],
          Column[{Style["Output weight",Bold,FontFamily->"Times"],
            Row[{"n' = ",row["nOut"]}],Row[{"P(n') = ",ScientificForm[v,4]}]},Spacings->.15]
        ]
      ],
      {rows,coords,vals}
    ],
    Lighting->{{"Ambient",GrayLevel[.72]},{"Directional",White,{2.3,-2.8,4.4}},
      {"Directional",RGBColor[.92,.97,1.],{-2.0,1.6,3.1}}}
  ];

  dominantIndex=First@Ordering[vals,-1];
  dominantCoord=coords[[dominantIndex]]; dominantVal=vals[[dominantIndex]];

  dominantGraphic=Graphics3D[
    {Directive[RGBColor[.83,.12,.08],PointSize[.012]],
     Point[{dominantCoord[[1]],dominantCoord[[2]],dominantVal}]},
    Lighting->"Neutral"
  ];

  selectedGraphic=Graphics3D[
    {Directive[RGBColor[.88,.05,.76],PointSize[.013],Opacity[.92]],
     Point[Append[selectedCoord,.018 vmax]]},
    Lighting->"Neutral"
  ];

  xRange=xr[[2]]-xr[[1]]; yRange=yr[[2]]-yr[[1]];
  xAxisY=yr[[1]]-.055 yRange; yAxisX=xr[[2]]+.055 xRange; axisZ=-.010 zmax;
  displayXr={xr[[1]]-.020 xRange,yAxisX+.085 xRange};
  displayYr={xAxisY-.085 yRange,yr[[2]]+.035 yRange};

  xTickVals=Select[{-1.,-.5,0.,.5,1.},xr[[1]]-10^-10<=#<=xr[[2]]+10^-10&];
  yTickVals=Select[{-1.,-.5,0.,.5,1.},yr[[1]]-10^-10<=#<=yr[[2]]+10^-10&];
  xTickSize=.012 yRange; yTickSize=.012 xRange;
  xLabelPos={Mean[xr],xAxisY-.075 yRange,axisZ};
  yLabelPos={yAxisX+.067 xRange,Mean[yr],axisZ};

  axesGraphic=Graphics3D[
    {
      {Directive[GrayLevel[.34],AbsoluteThickness[.95]],
       Line[{{xr[[1]],xAxisY,axisZ},{xr[[2]],xAxisY,axisZ}}],
       Line[{{yAxisX,yr[[1]],axisZ},{yAxisX,yr[[2]],axisZ}}]},
      Table[
        {Directive[GrayLevel[.42],AbsoluteThickness[.9]],
         Line[{{x,xAxisY-xTickSize,axisZ},{x,xAxisY+xTickSize,axisZ}}],
         Text[Style[NumberForm[x,{2,1}],8.7,GrayLevel[.20],FontFamily->"Times"],
           {x,xAxisY-.045 yRange,axisZ}]},
        {x,xTickVals}
      ],
      Table[
        {Directive[GrayLevel[.42],AbsoluteThickness[.9]],
         Line[{{yAxisX-yTickSize,y,axisZ},{yAxisX+yTickSize,y,axisZ}}],
         Text[Style[NumberForm[y,{2,1}],8.7,GrayLevel[.20],FontFamily->"Times"],
           {yAxisX+.040 xRange,y,axisZ}]},
        {y,yTickVals}
      ],
      Text[Style[TraditionalForm[Subscript[X,C]],11.3,GrayLevel[.10],FontFamily->"Times"],xLabelPos],
      Text[Style[TraditionalForm[Subscript[Y,C]],11.3,GrayLevel[.10],FontFamily->"Times"],yLabelPos]
    },
    Lighting->"Neutral"
  ];

  plot=Show[
    basePlane,bars,dominantGraphic,selectedGraphic,axesGraphic,
    Boxed->False,Axes->False,
    PlotRange->{displayXr,displayYr,{axisZ-.003 zmax,zmax}},
    PlotRangeClipping->False,PlotRangePadding->Scaled[.005],
    BoxRatios->{1.,.88,.44},ViewPoint->{1.42,-2.10,1.72},
    ViewVertical->{0,0,1},ViewAngle->.28,SphericalRegion->True,
    Background->White,ImageSize->$M3PaperMainImageSize,
    ImagePadding->{{16,16},{14,14}},BaseStyle->{FontFamily->"Times",FontSize->10.2}
  ];

  legendTicks=Table[
    {q,If[q==0,Style["0",8.3,FontFamily->"Times"],
      Style[ScientificForm[q,2],8.3,FontFamily->"Times"]]},
    {q,{0,vmax/4,vmax/2,3 vmax/4,vmax}}
  ];

  barLegend=BarLegend[
    {Function[q,colorF[q]],{0,vmax}},
    Ticks->legendTicks,
    LegendLabel->Style[TraditionalForm@HoldForm[P[Superscript[n,"\[Prime]"]]],9.5,FontFamily->"Times"],
    LabelStyle->Directive[FontFamily->"Times",FontSize->8.3,GrayLevel[.14]],
    LegendMarkerSize->$M3PaperLegendMarkerSize
  ];

  M3ReferencePaperCompose[plot,barLegend]
];


M3ReferenceGlobalProbabilityMap[data_Association,nN_List] := Module[
  {rows,Nrep,coords,vals,valid,selectedCoord,dominantIndex,dominantCoord,probLabel},

  rows=Lookup[data,"Results",{}];
  If[rows==={},Return[Graphics[Text[Style["There are no global probabilities to plot.",11,FontFamily->"Times"]],Background->White]]];

  Nrep=Total[nN];
  coords=M3CartanCoordinates3 /@ (N[Lookup[rows,"nOut"]]/Nrep);
  vals=Lookup[rows,"Probability",{}];
  vals=Map[Function[z,If[NumericQ[z],Max[0.,Chop[Re@N[z],10^-13]],Indeterminate]],vals];

  valid=Select[Range[Min[Length[coords],Length[vals]]],
    Function[k,VectorQ[coords[[k]],NumericQ]&&NumericQ[vals[[k]]]]];
  If[valid==={},Return[Graphics[Text[Style["There are no numerical probabilities to plot.",11,FontFamily->"Times"]],Background->White]]];

  selectedCoord=M3CartanCoordinates3[N[nN/Nrep]];
  dominantIndex=valid[[First@Ordering[N[vals[[valid]]],-1]]];
  dominantCoord=coords[[dominantIndex]];
  probLabel=TraditionalForm@HoldForm[P[Superscript[n,"\[Prime]"]]];

  M3ReferencePaperProbabilitySimplex[
    coords,vals,Nrep,selectedCoord,probLabel,Lookup[rows,"nOut",{}],Automatic,dominantCoord]
];


M3ReferenceGlobalDetectionExplorer[nN_List,Gamma_?MatrixQ,U_?MatrixQ,tol_:10^-9] := With[
  {nIn=nN,gam=N[Gamma],umat=N[U],rtol=tol},
  DynamicModule[{globalData=Missing["NotComputed"],status="Not computed"},
    Column[
      {
        Panel[
          Column[
            {
              Style["Total physical distribution: sum over all sectors",11.7,Bold,FontFamily->"Times"],
              Row[{
                Button[
                  "Compute total distribution",
                  status="Computing all sectors...";
                  globalData=M3ReferenceGlobalDetectionData[nIn,gam,umat,rtol];
                  status=If[AssociationQ[globalData],"Computation completed","The computation produced an error"],
                  Method->"Queued",ImageSize->190
                ],
                Spacer[12],
                Dynamic[Style[status,9.3,FontFamily->"Times",GrayLevel[.24]]]
              }]
            },
            Spacings->.45
          ],
          Background->White,FrameMargins->{{10,10},{8,8}},
          FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.65]]
        ],

        Dynamic[
          Which[
            MissingQ[globalData],
            Panel[
              Style["Press 'Compute total distribution' to sum the contributions P(Lambda,n') from all compatible sectors.",
                9.2,FontFamily->"Times",GrayLevel[.30]],
              Background->RGBColor[.985,.99,1.],FrameMargins->10,
              FrameStyle->Directive[GrayLevel[.88],AbsoluteThickness[.55]]
            ],

            FailureQ[globalData]||!AssociationQ[globalData],
            M3FailurePanel[globalData,"The total physical distribution could not be constructed"],

            True,
            Column[
              {
                M3ScientificResultPanel[
                  "Global-sum check",
                  "The total distribution must be normalized after summing all sectors.",
                  Grid[
                    {
                      {Style["compatible sectors",Bold],Length[Lookup[globalData,"SectorData",{}]],
                       Style["Sum P(n')",Bold],NumberForm[Lookup[globalData,"TotalProbability",Indeterminate],{10,8}]},
                      {Style["normalization error",Bold],ScientificForm[Lookup[globalData,"NormalizationError",Indeterminate],3],
                       SpanFromLeft,SpanFromLeft}
                    },
                    Alignment->{{Left,Right,Left,Right}},Spacings->{1.0,.55},
                    BaseStyle->{FontFamily->"Times",FontSize->9.3}
                  ]
                ],
                Grid[{{M3ReferenceGlobalProbabilityPlot3D[globalData,nIn],
                  M3ReferenceGlobalProbabilityMap[globalData,nIn]}},
                  Alignment->{Center,Top},Spacings->{.55,.45},ItemSize->All]
              },
              Alignment->Center,Spacings->.55
            ]
          ],
          TrackedSymbols:>{globalData,status}
        ]
      },
      Alignment->Center,Spacings->.55
    ]
  ]
];


M3ReferenceIRRExplorer[
  nN_List,Gamma_?MatrixQ,U_?MatrixQ,referenceLambda_List,tol_:10^-9] := Module[
  {sectors,init},
  sectors=M3ReferenceCompatibleSectors[nN];
  If[sectors==={},Return[Panel["There are no sectors compatible with the input weight."]]];
  init=First@FirstPosition[sectors,referenceLambda,{1}];

  With[{sectorList=sectors,nIn=nN,gam=N[Gamma],umat=N[U],start=init,rtol=tol},
    DynamicModule[{sectorIndex=start,mode="Joint"},
      Column[{
        Panel[
          Dynamic@Row[{
            Style["Sector No.",Bold,FontFamily->"Times"],Spacer[8],
            Slider[Dynamic[sectorIndex],{1,Length[sectorList],1},ImageSize->165],Spacer[7],
            Style[Round[sectorIndex],10,Bold,FontFamily->"Times"],Spacer[20],
            Style["Highest weight  Lambda = ",Bold,FontFamily->"Times"],
            Style[sectorList[[Clip[Round[sectorIndex],{1,Length[sectorList]}]]],10.4,Bold,FontFamily->"Times"],
            Spacer[18],
            SetterBar[Dynamic[mode],{"Joint"->"P(Lambda,n')","Conditional"->"P_Lambda(n')"}]
          }],
          Background->White,FrameMargins->{{10,10},{7,7}},
          FrameStyle->Directive[GrayLevel[.84],AbsoluteThickness[.65]]
        ],

        Dynamic[
          Module[{idx,lamN,data,leftPanel,rightPanel,bottomPanel},
            idx=Clip[Round[sectorIndex],{1,Length[sectorList]}];
            lamN=sectorList[[idx]];
            data=M3ReferenceSectorData[lamN,nIn,gam,umat,rtol];
            If[FailureQ[data]||!AssociationQ[data],
              M3FailurePanel[data,"The selected sector could not be evaluated"],
              leftPanel=M3ScientificResultPanel[
                "Sector diagram, input weight, and highest weight",
                "The cyan region is the irrep \[CapitalLambda]",
                M3ReferenceWeightDiagram[lamN,nIn]
              ];
              rightPanel=M3ReferenceSectorProbabilityPlot3D[data,mode];
              bottomPanel=M3ReferenceSectorProbabilityMap[data,mode];
              Grid[
                {
                  {leftPanel,rightPanel},
                  {Item[bottomPanel,Alignment->Center],SpanFromLeft}
                },
                Alignment->{Center,Top},
                Spacings->{.55,.45},
                ItemSize->All
              ]
            ]
          ],
          TrackedSymbols:>{sectorIndex,mode}
        ],

        M3ReferenceGlobalDetectionExplorer[nIn,gam,umat,rtol]
      },Alignment->Center,Spacings->.55]
    ]
  ]
];


(* Orthonormal physical-frame basis \[LongDash] minimize E = 1/2 ||G - I_d||_F^2 over admissible frames. *)
(* diagnostics for the criterion G = I_d*)
M3QBMetrics[coordinates_List,tol_:10^-8]:=Module[ {coords,C,G,d,off,energy,frob},

If[coordinates==={}, Return[Failure["EmptyCoordinates",<|"Message"->"There are no coordinates to evaluate the orthonormal-basis question."|>]] ];

coords=Table[ If[!VectorQ[N[c],NumericQ]||Norm[N[c]]<=tol,
Return[Failure["InvalidCoordinateVector",<|"Message"->"A zero or nonnumeric coordinate vector was found."|>]] ]; N[c/Norm[c]], {c,coordinates} ];

d=Length[coords]; C=Transpose[coords]; G=Chop[ConjugateTranspose[C] . C,10^-13]; off=G-DiagonalMatrix[Diagonal[G]];

energy=N[Total[Abs[Flatten[off]]^2]/2]; frob=N[Sqrt[Total[Abs[Flatten[G-IdentityMatrix[d]]]^2]]];

<| "Coordinates"->coords, "Gram"->G, "Energy"->energy, "FrobeniusToIdentity"->frob |> ];

M3QBEnergyCoordinates[coordinates_List]:=Module[{coords,C,G,off}, If[coordinates==={},Return[Infinity]];
If[!And@@(VectorQ[N[#],NumericQ]&&Norm[N[#]]>10^-14&/@coordinates),Return[Infinity]]; coords=N[#/Norm[#]&/@coordinates]; C=Transpose[coords];
G=ConjugateTranspose[C] . C; off=G-DiagonalMatrix[Diagonal[G]]; N[Total[Abs[Flatten[off]]^2]/2] ];

(* Removes projectively duplicated rays: |<c_a,c_b>|^2 ~= 1. *)
M3QBUniqueRecords[records_List,tol_:10^-10]:=Module[{valid,out={},c}, valid=Select[ records, Function[r, AssociationQ[r]&&KeyExistsQ[r,"Coordinates"]&&
VectorQ[N[r["Coordinates"]],NumericQ]&&Norm[N[r["Coordinates"]]]>tol ] ];

Do[ c=N[rec["Coordinates"]/Norm[rec["Coordinates"]]]; If[ !AnyTrue[out,Function[r,Abs[Conjugate[r["Coordinates"]] . c]^2>=1-tol]],
out=Append[out,Join[rec,<|"Coordinates"->c|>]] ], {rec,valid} ]; out ];

M3QBSubsetEnergy[G_?MatrixQ,inds_List]:=Module[{sub,off}, sub=G[[inds,inds]]; off=sub-DiagonalMatrix[Diagonal[sub]]; N[Total[Abs[Flatten[off]]^2]/2] ];
M3QBSeedSubsets[ records_List,d_Integer?Positive,nSeeds_Integer:8,exhaustiveLimit_Integer:120000 ]:=Module[
{K,coords,C,G,score,nComb,allSubs,seeds,greedyFromStart,starts, randomSubs,nRandom,candidates,improveOne,current,remaining,proposals,best,
changed,seedValue},

K=Length[records]; If[K<d, Return[Failure["InsufficientFrames",<|
"Message"->"The current bank contains fewer distinct physical rays than the multiplicity. Increase 'maximum frames' and/or 'fiber resolution'.",
"Frames"->K,"Multiplicity"->d|>]] ];

coords=N[#["Coordinates"]&/@records]; C=Transpose[coords]; G=Chop[ConjugateTranspose[C] . C,10^-13]; score[inds_List]:=M3QBSubsetEnergy[G,inds];
nComb=Binomial[K,d];

If[nComb<=exhaustiveLimit, allSubs=Subsets[Range[K],{d}]; seeds=TakeSmallestBy[allSubs,score,Min[nSeeds,Length[allSubs]]];,

greedyFromStart[s_Integer]:=Module[{sel={s},rem,pick}, rem=DeleteCases[Range[K],s]; While[Length[sel]<d,
pick=First@MinimalBy[rem,Function[j,Total[Abs[G[[j,sel]]]^2]]]; sel=Append[sel,pick]; rem=DeleteCases[rem,pick] ]; Sort[sel] ];

starts=DeleteDuplicates[greedyFromStart/@Range[K]]; nRandom=Min[800,Max[150,25 K]]; seedValue=Hash[{K,d,Round[10^10 Total[Abs[Flatten[G]]]]}];
randomSubs=BlockRandom[ SeedRandom[seedValue]; Table[Sort@RandomSample[Range[K],d],{nRandom}] ];

candidates=DeleteDuplicates@Join[starts,randomSubs]; candidates=TakeSmallestBy[candidates,score,Min[Max[6 nSeeds,40],Length[candidates]]];

improveOne[s0_List]:=Module[{cur=s0,rem,props,cand,iter=0}, changed=True; While[TrueQ[changed]&&iter<4, iter++; changed=False;
rem=Complement[Range[K],cur]; props=DeleteDuplicates@Flatten[ Table[Sort[ReplacePart[cur,pos->j]],{pos,Length[cur]},{j,rem}],1 ]; If[props=!={},
cand=First@MinimalBy[props,score]; If[score[cand]<score[cur]-10^-13, cur=cand; changed=True ] ] ]; cur ];

seeds=DeleteDuplicates[improveOne/@candidates]; seeds=TakeSmallestBy[seeds,score,Min[nSeeds,Length[seeds]]] ];

<| "Seeds"->seeds, "SeedEnergies"->(score/@seeds), "GramBank"->G, "NumberOfCombinations"->nComb, "Method"->If[nComb<=exhaustiveLimit,
"exact combinatorics + global cloud + annealing + polishing", "greedy/random/swaps + global cloud + annealing + polishing" ] |> ];

(* physical frame at a continuous coordinate (x,y) *)
M3QBRecordAtPoint[ lam_List,n_List,finite_Association,iso_Association, point_List,branch_,tol_:10^-9 ]:=Module[{x,y,Q,H,b,V,wd,c,fres,uerr},
If[Length[point]=!=2,Return[$Failed]]; {x,y}=N[point]; If[!NumericQ[x]||!NumericQ[y]||x<0||x>1||y<0||y>1,Return[$Failed]];

Q=M3QSection[lam,n,x,y]; If[FailureQ[Q]||!MatrixQ[Q],Return[$Failed]]; If[!M3BistochasticSectionQ[Q,100 tol],Return[$Failed]]; H=N[M3HeronPolynomialQ[Q]];
If[!NumericQ[H]||H< -100 tol,Return[$Failed]];

b=Which[ branch===-1,-1, branch===1,1, StringContainsQ[ToString[branch],"-"],-1, True,1 ];

V=Quiet[Check[M3ReconstructUnitary3[Q,b,tol],$Failed]]; If[!MatrixQ[V],Return[$Failed]];

fres=N[Norm[Q . N[lam]-N[n]]]; uerr=N[Sqrt[Total[Abs[Flatten[ConjugateTranspose[V] . V-IdentityMatrix[3]]]^2]]];
If[!NumericQ[fres]||!NumericQ[uerr]||fres>10^-6||uerr>10^-6,Return[$Failed]];

wd=Quiet[Check[M3WignerCoordinates[finite,iso,V,tol],$Failed],{Power::indet,Infinity::indet}]; If[!AssociationQ[wd],Return[$Failed]];
c=N[wd["Coordinates"]]; If[!VectorQ[c,NumericQ]||Norm[c]<=tol,Return[$Failed]]; c=N[c/Norm[c]];

<| "Point"->{x,y}, "Branch"->If[b==1,"+","-"], "Q"->Q, "V"->V, "Heron"->H, "Coordinates"->c, "FiberResidual"->fres, "UnitarityError"->uerr,
"CoordinateSource"->"continuous Wigner-D" |> ];
M3QBQuasiPoints[count_Integer?Positive]:=Module[{a,b}, a=(Sqrt[5]-1)/2; b=Sqrt[2]-1; Table[ {FractionalPart[(k+1/3) a],FractionalPart[(k+1/7) b]},
{k,1,count} ] ];

M3QBBuildCloud[ lam_List,n_List,finite_Association,iso_Association,baseRecords_List, pointCount_Integer?Positive,tol_:10^-9
]:=Module[{pts,raw,attempted=0,r}, pts=M3QBQuasiPoints[pointCount]; raw=Reap[ Do[ Do[ attempted++; r=M3QBRecordAtPoint[lam,n,finite,iso,p,b,tol];
If[AssociationQ[r],Sow[r]], {b,{1,-1}} ], {p,pts} ] ]; raw=If[Length[raw[[2]]]>=1,First[raw[[2]]],{}]; <| "Records"->Join[baseRecords,raw],
"AttemptedEvaluations"->attempted, "ValidContinuousRecords"->Length[raw] |> ];
M3QBCloudDescent[seedRecords_List,cloudRecords_List,maxSweeps_Integer:8]:=Module[ {recs,coords,cloud,cloudCoords,d,energy,startEnergy,accepted=0,sweep,
contrib,order,i,others,oldLocal,scores,k,newLocal,candC,duplicateQ,changed},

recs=seedRecords; coords=N[#["Coordinates"]/Norm[#["Coordinates"]]&/@recs]; cloud=Select[cloudRecords,AssociationQ[#]&&KeyExistsQ[#,"Coordinates"]&];
cloudCoords=N[#["Coordinates"]/Norm[#["Coordinates"]]&/@cloud]; d=Length[coords]; energy=M3QBEnergyCoordinates[coords]; startEnergy=energy;

If[d<=1, Return[<|"Records"->recs,"Coordinates"->coords,"StartEnergy"->energy,"FinalEnergy"->energy,"AcceptedMoves"->0,"Sweeps"->0|>] ];

sweep=0; changed=True; While[TrueQ[changed]&&sweep<maxSweeps, sweep++; changed=False; contrib=Table[
Total@Table[If[j==i,0.,Abs[Conjugate[coords[[i]]] . coords[[j]]]^2],{j,d}], {i,d} ]; order=Reverse@Ordering[contrib];

Do[ others=Delete[coords,i]; oldLocal=Total[Abs[Conjugate[coords[[i]]] . #]^2&/@others];

scores=Table[ candC=cloudCoords[[kk]]; duplicateQ=AnyTrue[others,Abs[Conjugate[candC] . #]^2>1-10^-9&];
If[duplicateQ,Infinity,Total[Abs[Conjugate[candC] . #]^2&/@others]], {kk,Length[cloudCoords]} ];

k=First@Ordering[scores,1]; newLocal=scores[[k]]; If[NumericQ[newLocal]&&newLocal<oldLocal-10^-12, recs[[i]]=cloud[[k]]; coords[[i]]=cloudCoords[[k]];
energy=N[energy-oldLocal+newLocal]; accepted++; changed=True ], {i,order} ] ];

<| "Records"->recs, "Coordinates"->coords, "StartEnergy"->startEnergy, "FinalEnergy"->M3QBEnergyCoordinates[coords], "AcceptedMoves"->accepted,
"Sweeps"->sweep |> ];
M3QBAnnealSeed[ lam_List,n_List,finite_Association,iso_Association,seedRecords_List, steps_Integer:1400,tol_:10^-9,randomSeed_Integer:12345 ]:=Module[
{recs,coords,bestRecs,bestCoords,d,energy,bestEnergy,cache=<||>,evals=0, accepted=0,key,evalPoint,t,frac,step,temp,contrib,order,i,pt,b,newB,newPt,
cand,candC,others,oldLocal,newLocal,newEnergy,delta,acceptQ},

recs=seedRecords; If[!And@@(AssociationQ/@recs),Return[Failure["InvalidAnnealSeed",<|"Message"->"Invalid seed for annealing."|>]]];
If[!And@@(ListQ[Lookup[#,"Point",None]]&/@recs), Return[Failure["NoContinuousChart",<|"Message"->"The seed has no continuous (x,y) coordinates."|>]] ];

coords=N[#["Coordinates"]/Norm[#["Coordinates"]]&/@recs]; d=Length[coords]; energy=M3QBEnergyCoordinates[coords]; bestEnergy=energy; bestRecs=recs;
bestCoords=coords;

key[p_List,bb_Integer]:=ToString[InputForm[{Round[N[p[[1]]],10^-9],Round[N[p[[2]]],10^-9],bb}]]; evalPoint[p_List,bb_Integer]:=Module[{k,r}, k=key[p,bb];
If[KeyExistsQ[cache,k],Return[cache[k]]]; evals++; r=M3QBRecordAtPoint[lam,n,finite,iso,p,bb,tol]; AssociateTo[cache,k->r]; r ];

BlockRandom[ SeedRandom[randomSeed]; Do[ frac=If[steps<=1,1.,N[(t-1)/(steps-1)]]; step=0.045 (0.00045/0.045)^frac; temp=0.020 (2.*10^-5/0.020)^frac;

contrib=Table[ Total@Table[If[j==ii,0.,Abs[Conjugate[coords[[ii]]] . coords[[j]]]^2],{j,d}], {ii,d} ]; order=Reverse@Ordering[contrib];
i=If[RandomReal[]<0.72, RandomChoice[Take[order,UpTo[Min[3,d]]]], RandomInteger[{1,d}] ];

pt=N[recs[[i]]["Point"]]; b=If[StringContainsQ[ToString[recs[[i]]["Branch"]],"-"],-1,1]; newB=If[RandomReal[]<0.14,-b,b]; newPt=If[RandomReal[]<0.06, pt,
(Clip[#,{0.,1.}]&/@(pt+step RandomReal[{-1,1},2])) ];

cand=evalPoint[newPt,newB]; If[AssociationQ[cand], candC=N[cand["Coordinates"]/Norm[cand["Coordinates"]]]; others=Delete[coords,i];
If[!AnyTrue[others,Abs[Conjugate[candC] . #]^2>1-10^-9&], oldLocal=Total[Abs[Conjugate[coords[[i]]] . #]^2&/@others];
newLocal=Total[Abs[Conjugate[candC] . #]^2&/@others]; newEnergy=N[energy-oldLocal+newLocal]; delta=newEnergy-energy;
acceptQ=TrueQ[delta<=0]|| Log[Max[RandomReal[],10.^-300]] < -delta/Max[temp,10^-12]; If[acceptQ, recs[[i]]=cand; coords[[i]]=candC; energy=newEnergy;
accepted++; If[energy<bestEnergy, bestEnergy=energy; bestRecs=recs; bestCoords=coords ] ] ] ], {t,1,steps} ] ];

<| "Records"->bestRecs, "Coordinates"->bestCoords, "FinalEnergy"->bestEnergy, "CoordinateEvaluations"->evals, "AcceptedMoves"->accepted |> ];
M3QBLocalPolish[ lam_List,n_List,finite_Association,iso_Association,seedRecords_List,
stepSchedule_List:{0.030,0.018,0.010,0.006,0.0035,0.0020,0.0010,0.0005,0.00025}, maxSweeps_Integer:2,tol_:10^-9 ]:=Module[
{recs,coords,d,energy,startEnergy,cache=<||>,evals=0,accepted=0, key,evalPoint,dirs,step,sweep,changed,contrib,order,i,pt,b,specs,bestRec,
bestCoord,bestE,others,oldLocal,cand,candC,newLocal,newE,newPt},

recs=seedRecords; If[!And@@(AssociationQ/@recs),Return[Failure["InvalidPolishSeed",<|"Message"->"Invalid seed for polishing."|>]]];
If[!And@@(ListQ[Lookup[#,"Point",None]]&/@recs), Return[Failure["NoContinuousChart",<|"Message"->"The seed has no (x,y) coordinates."|>]] ];

coords=N[#["Coordinates"]/Norm[#["Coordinates"]]&/@recs]; d=Length[coords]; energy=M3QBEnergyCoordinates[coords]; startEnergy=energy;

key[p_List,bb_Integer]:=ToString[InputForm[{Round[N[p[[1]]],10^-10],Round[N[p[[2]]],10^-10],bb}]]; evalPoint[p_List,bb_Integer]:=Module[{k,r}, k=key[p,bb];
If[KeyExistsQ[cache,k],Return[cache[k]]]; evals++; r=M3QBRecordAtPoint[lam,n,finite,iso,p,bb,tol]; AssociateTo[cache,k->r]; r ];

dirs=Table[{Cos[2 Pi k/16.],Sin[2 Pi k/16.]},{k,0,15}];

Do[ sweep=0; changed=True; While[TrueQ[changed]&&sweep<maxSweeps, sweep++; changed=False; contrib=Table[
Total@Table[If[j==ii,0.,Abs[Conjugate[coords[[ii]]] . coords[[j]]]^2],{j,d}], {ii,d} ]; order=Reverse@Ordering[contrib];

Do[ pt=N[recs[[i]]["Point"]]; b=If[StringContainsQ[ToString[recs[[i]]["Branch"]],"-"],-1,1]; others=Delete[coords,i];
oldLocal=Total[Abs[Conjugate[coords[[i]]] . #]^2&/@others]; bestRec=recs[[i]]; bestCoord=coords[[i]]; bestE=energy;

specs=DeleteDuplicates@Join[ Table[{(Clip[#,{0.,1.}]&/@(pt+step dir)),b},{dir,dirs}], Table[{(Clip[#,{0.,1.}]&/@(pt+step dir)),-b},{dir,dirs}], {{pt,-b}}
];

Do[ newPt=spec[[1]]; cand=evalPoint[newPt,spec[[2]]]; If[AssociationQ[cand], candC=N[cand["Coordinates"]/Norm[cand["Coordinates"]]];
If[!AnyTrue[others,Abs[Conjugate[candC] . #]^2>1-10^-9&], newLocal=Total[Abs[Conjugate[candC] . #]^2&/@others]; newE=N[energy-oldLocal+newLocal];
If[newE<bestE-10^-13, bestE=newE; bestRec=cand; bestCoord=candC ] ] ], {spec,specs} ];

If[bestE<energy-10^-13, recs[[i]]=bestRec; coords[[i]]=bestCoord; energy=bestE; accepted++; changed=True ], {i,order} ] ], {step,stepSchedule} ];

<| "Records"->recs, "Coordinates"->coords, "StartEnergy"->startEnergy, "FinalEnergy"->M3QBEnergyCoordinates[coords], "CoordinateEvaluations"->evals,
"AcceptedMoves"->accepted |> ];
M3QBSearchHybrid[ lam_List,n_List,finite_Association,iso_Association,records_List, tol_:10^-8,nSeeds_Integer:8 ]:=Module[
{d,unique,seedData,seeds,bankRecords,bankEnergy,nondegenerate, cloudPointCount,cloudData,cloud,cloudRuns,cloudBest,cloudEnergy,
annealInputs,annealRuns={},ann,annBest, polishInputs,polishRuns={},pol,best,finalRecords,met, finalCandidates,seedHash},

d=iso["Multiplicity"]; If[d<=0, Return[Failure["ZeroMultiplicity",<|"Message"->"The selected weight has zero multiplicity."|>]] ];

unique=M3QBUniqueRecords[records,10^-10]; If[Length[unique]<d, Return[Failure["InsufficientQuestionBBank",<|
"Message"->"The current bank does not contain enough distinct rays. Increase 'maximum frames' or 'fiber resolution'.",
"UniqueFrames"->Length[unique],"Multiplicity"->d|>]] ];

seedData=M3QBSeedSubsets[unique,d,nSeeds]; If[FailureQ[seedData],Return[seedData]]; seeds=seedData["Seeds"]; If[!ListQ[seeds]||seeds==={},
Return[Failure["NoQuestionBSeeds",<|"Message"->"Seeds for the orthonormal-basis question could not be constructed."|>]] ];

bankRecords=unique[[First[seeds]]]; bankEnergy=M3QBEnergyCoordinates[(#["Coordinates"]&/@bankRecords)];

(* d=1 is trivial. *)
If[d==1, met=M3QBMetrics[(#["Coordinates"]&/@bankRecords),Min[tol,10^-8]]; Return[<|"Metrics"->met|>] ];

nondegenerate=Abs[N[lam[[2]]-lam[[3]]]]>100 tol;

If[nondegenerate, cloudPointCount=Max[1200,220 d]; cloudData=M3QBBuildCloud[N[lam],N[n],finite,iso,unique,cloudPointCount,Max[10^-10,tol]];
cloud=cloudData["Records"];

cloudRuns=Table[ M3QBCloudDescent[unique[[seed]],cloud,8], {seed,seeds} ];
cloudRuns=Select[cloudRuns,AssociationQ[#]&&NumericQ[Lookup[#,"FinalEnergy",Infinity]]&];

If[cloudRuns==={}, cloudBest=<|"Records"->bankRecords,"FinalEnergy"->bankEnergy|>; cloudEnergy=bankEnergy;,
cloudBest=First@MinimalBy[cloudRuns,Lookup[#,"FinalEnergy",Infinity]&]; cloudEnergy=cloudBest["FinalEnergy"] ];

annealInputs=TakeSmallestBy[cloudRuns,Lookup[#,"FinalEnergy",Infinity]&,Min[4,Length[cloudRuns]]]; If[annealInputs==={},annealInputs={cloudBest}];

Do[ seedHash=Hash[{Round[10^8 Lookup[cr,"FinalEnergy",1.]],Length[annealRuns],finite["LambdaN"],finite["nN"]}]; ann=Quiet[Check[
M3QBAnnealSeed[N[lam],N[n],finite,iso,cr["Records"],1500,Max[10^-10,tol],seedHash], $Failed ],{Power::indet,Infinity::indet}];
If[AssociationQ[ann],AppendTo[annealRuns,ann]], {cr,annealInputs} ];

If[annealRuns==={}, annBest=<|"Records"->cloudBest["Records"],"FinalEnergy"->cloudEnergy|>;,
annBest=First@MinimalBy[annealRuns,Lookup[#,"FinalEnergy",Infinity]&] ];

polishInputs=DeleteDuplicatesBy[ Join[ TakeSmallestBy[annealRuns,Lookup[#,"FinalEnergy",Infinity]&,Min[2,Length[annealRuns]]],
{annBest,<|"Records"->cloudBest["Records"],"FinalEnergy"->cloudEnergy|>} ], Round[10^10 Lookup[#,"FinalEnergy",Infinity]]& ];

Do[ pol=Quiet[Check[ M3QBLocalPolish[N[lam],N[n],finite,iso,pr["Records"],
{0.030,0.018,0.010,0.006,0.0035,0.0020,0.0010,0.0005,0.00025},2,Max[10^-10,tol]], $Failed ],{Power::indet,Infinity::indet}];
If[AssociationQ[pol],AppendTo[polishRuns,pol]], {pr,polishInputs} ];

finalCandidates=Join[ polishRuns, annealRuns, cloudRuns, {<|"Records"->bankRecords,"FinalEnergy"->bankEnergy|>} ]; finalCandidates=Select[finalCandidates,
AssociationQ[#]&&KeyExistsQ[#,"Records"]&&NumericQ[Lookup[#,"FinalEnergy",Infinity]]& ];
best=First@MinimalBy[finalCandidates,Lookup[#,"FinalEnergy",Infinity]&]; finalRecords=best["Records"];,

(* degenerate spectrum: keep the best result from the bank. *)
cloudEnergy=bankEnergy; finalRecords=bankRecords ];

met=M3QBMetrics[(#["Coordinates"]&/@finalRecords),Min[tol,10^-8]]; If[FailureQ[met],Return[met]]; <|"Metrics"->met|> ];

M3QBCacheKey[ lam_List,n_List,finite_Association,iso_Association,records_List,tol_ ]:=Hash[ { Chop[N[lam],10^-12],Chop[N[n],10^-12],
finite["LambdaN"],finite["nN"],iso["Multiplicity"], Table[ {Lookup[r,"Point",None],Lookup[r,"Branch",None],Chop[N[Lookup[r,"Coordinates",{}]],10^-10]},
{r,records} ], tol,"QB-hybrid-compact-v7" }, "SHA256" ];

M3QBSearchCached[
  lam_List,
  n_List,
  finite_Association,
  iso_Association,
  records_List,
  tol_:10^-8
] := Module[{key, res},

  If[!AssociationQ[$M3QBCache],
    $M3QBCache = <||>
  ];

  key = M3QBCacheKey[
    lam,
    n,
    finite,
    iso,
    records,
    tol
  ];

  If[
    KeyExistsQ[$M3QBCache, key],
    Return[$M3QBCache[key]]
  ];

  res = M3QBSearchHybrid[
    lam,
    n,
    finite,
    iso,
    records,
    tol,
    8
  ];

  AssociateTo[
    $M3QBCache,
    key -> res
  ];

  res
];

M3QBResultPanel[result_Association]:=Module[ {met,G,summary},

met=result["Metrics"]; G=met["Gram"];

summary=Framed[ Grid[{ {Style["final E",Bold],ScientificForm[met["Energy"],6]}, {Style["||G-I||_F",Bold],ScientificForm[met["FrobeniusToIdentity"],6]} },
Frame->All, Alignment->{{Left,Right}}, Spacings->{1.5,.6}, BaseStyle->{FontFamily->"Times",FontSize->10} ], Background->White,
FrameStyle->Directive[GrayLevel[.78],AbsoluteThickness[.8]], FrameMargins->{{12,12},{9,9}}, RoundingRadius->4 ];

Column[{ summary, M3GramHeatmap[G] },Alignment->Center,Spacings->.8] ];

M3QBExplorer[ lam_List,n_List,finite_Association,iso_Association,records_List,rankTol_:10^-8 ]:=Module[{result},
result=M3QBSearchCached[N[lam],N[n],finite,iso,records,rankTol];

Column[{ Style["orthonormal basis of physically admissible frames",15,Bold, FontFamily->"Times",FontColor->RGBColor[.08,.20,.38]], Style[
"Goal: find d=m_Lambda(n) physical frames such that G_ab=K(Va,Vb)=delta_ab.", 10.2,FontFamily->"Times",GrayLevel[.28] ], If[AssociationQ[result],
M3QBResultPanel[result], M3FailurePanel[result,"The search for the orthonormal-basis question could not be completed"] ] },Alignment->Center,Spacings->.65]
];

ClearAll[M3SU3Interface]; M3SU3Interface[]:=Manipulate[ Module[ { lam,n,nmin,Nrep,finite,iso,gt,multDiag,rawRecords,records,coordinates,C,gramDiag,G,
independentIndices,frameIndex,frameChoices,rec,V,Q,orbitCoords, header,tabbar,content,phys,U,transitions,detPhys,detOrb,comparison,
weightMapData,tabs1,tabs2 }, lam=M3LambdaFromIntegerRatios[lamA,lamB,lamC]; n=M3NbarFromIntegerRatios[nA,nB,nC];

If[!TrueQ[M3MajorizedQ3[lam,n]], Return@Panel[ Column[{ Style["Physically inadmissible data",14,Bold,RGBColor[.65,.10,.10],FontFamily->"Times"],
Style[Row[{TraditionalForm[Overscript["n","_"]]," is not majorized by ",TraditionalForm[Overscript[\[Lambda],"_"]],"."}],10,FontFamily->"Times"],
Row[{Style[TraditionalForm[Overscript[\[Lambda],"_"]],Bold]," = ",N[lam],Spacer[18],Style[TraditionalForm[Overscript["n","_"]],Bold]," = ",N[n]}]
},Spacings->.7],Background->RGBColor[1.,.975,.975],FrameMargins->15 ] ];

nmin=M3MinimalCompatibleN[lam,n]; If[FailureQ[nmin]||!IntegerQ[nmin]||nmin<=0, Return[M3FailurePanel[nmin,"N_min could not be determined"]] ];
Nrep=scaleFactor nmin;

finite=M3FiniteRepresentationData[lam,n,Nrep];
If[FailureQ[finite]||!AssociationQ[finite],Return[M3FailurePanel[finite,"The finite representation could not be constructed"]]];

iso=M3IsospinData[finite]; If[FailureQ[iso]||!AssociationQ[iso],Return[M3FailurePanel[iso,"The isospin basis could not be constructed"]]];

gt=M3GTData[finite]; If[FailureQ[gt]||!AssociationQ[gt],Return[M3FailurePanel[gt,"Gelfand-Tsetlin could not be constructed"]]];

multDiag=M3MultiplicityDiagnostics[finite,iso];
If[FailureQ[multDiag]||!AssociationQ[multDiag],Return[M3FailurePanel[multDiag,"Multiplicity comparison failed"]]];

(* frame bank *)
rawRecords=M3CandidateFrames[N[lam],N[n],maxFrames,gridN,includeConjugate,rankTol];
If[FailureQ[rawRecords],Return[M3FailurePanel[rawRecords,"The frame fiber could not be constructed"]]];
If[!ListQ[rawRecords]||rawRecords==={},Return[Panel["No physically admissible frames were found for these parameters."]]];

records=DeleteCases[ Table[ Module[{wd,ex,c,source}, wd=M3WignerCoordinates[finite,iso,r0["V"],rankTol]; Which[ AssociationQ[wd],
c=wd["Coordinates"];source="Wigner-D"; Join[r0,<|"Coordinates"->c,"CoordinateSource"->source|>], True,
ex=M3ExactOrbitCoordinates[finite,iso,r0["V"],rankTol]; If[AssociationQ[ex], c=ex["Coordinates"];source="exact ambient";
Join[r0,<|"Coordinates"->c,"CoordinateSource"->source|>], Nothing ] ] ], {r0,rawRecords} ],Nothing ];

If[records==={},Return[M3FailurePanel[Failure["NoValidCoordinates",<|"Message"->"No frame produced valid orbital coordinates."|>],"The projected states could not be constructed"]]];

coordinates=(#["Coordinates"]&/@records); C=M3CoordinateMatrix[coordinates];
If[FailureQ[C]||!MatrixQ[C],Return[M3FailurePanel[C,"C could not be constructed"]]]; gramDiag=M3GramDiagnostics[C,rankTol];
If[FailureQ[gramDiag]||!AssociationQ[gramDiag],Return[M3FailurePanel[gramDiag,"Gram could not be constructed"]]]; G=gramDiag["Gram"];
independentIndices=gramDiag["IndependentFrameIndices"]; If[selectedFrame<1||selectedFrame>Length[records],selectedFrame=1];
frameIndex=Clip[selectedFrame,{1,Length[records]}]; rec=records[[frameIndex]];V=rec["V"];Q=rec["Q"];
orbitCoords=N[rec["Coordinates"]/Norm[rec["Coordinates"]]];

frameChoices=Table[i->Row[{"frame ",i,"  ",records[[i]]["Branch"]}],{i,Length[records]}]; header=Panel[ Column[{
Style["SU(3): geometry, multiplicity, and probabilities",17,Bold,FontFamily->"Times"], Grid[{ {
Style[TraditionalForm[Overscript[\[Lambda],"_"]],Bold],NumberForm[N[lam],{6,4}], Style[TraditionalForm[Overscript["n","_"]],Bold],NumberForm[N[n],{6,4}],
Style[Subscript["N","min"],Bold],nmin,Style["N",Bold],Nrep }, { Style[TraditionalForm[\[CapitalLambda]],Bold],finite["LambdaN"],
Style["(p,q,r)",Bold],{finite["p"],finite["q"],finite["r"]},
Style["W",Bold],finite["W"],Style[Row[{Subscript["m",\[CapitalLambda]],"(n)"}],Bold],iso["Multiplicity"] }, {
Style["frames",Bold],Length[records],Style["rank C = rank G",Bold],gramDiag["RankC"],
Style["bank status",Bold],If[gramDiag["RankC"]==iso["Multiplicity"],Style["complete",Darker[Green],Bold],Style["incomplete",Darker[Orange],Bold]],SpanFromLeft
} },Frame->All,Alignment->Left, Background->{1->{RGBColor[.955,.97,.99]}}, BaseStyle->{FontFamily->"Times",FontSize->9.8},Spacings->{1.05,.55}], Grid[{{
Style["Frame",Bold],PopupMenu[Dynamic[selectedFrame],frameChoices] }},Alignment->Left] },Spacings->.55], Background->White,FrameMargins->12,
FrameStyle->Directive[GrayLevel[.80],AbsoluteThickness[.7]] ];
tabs1={"Summary","Geometry","Multiplicity / Gram / Cartan","Wigner-D / isospin","Saddle / Hessian","\[Sigma] / probabilities"};
tabs2={"Detection","Toroidal map","IRR decomposition / probabilities","orthonormal basis"}; tabbar=Panel[ Column[{ SetterBar[Dynamic[activeTab],tabs1],
SetterBar[Dynamic[activeTab],tabs2] },Spacings->.28], Background->RGBColor[.965,.975,.99],FrameMargins->7,
FrameStyle->Directive[GrayLevel[.82],AbsoluteThickness[.6]] ];

content=Switch[activeTab,

"Summary", Column[{M3FrameTable[records,independentIndices] },Spacings->.75],

"Geometry", Module[{geometryPlots}, geometryPlots=M3GeometryBasePlots[N[lam],N[n]]; Column[{ Grid[{ {
M3FiberFramesFigure[N[lam],N[n],records,independentIndices,geometryPlots], M3UnitarityTriangleFigure[Q] } },Alignment->{Center,Top},Spacings->{1.0,0}],
M3SelectedFrameSummaryPanel[N[lam],N[n],Q,V], M3BistochasticUnistochasticFigure[N[lam],N[n],records,geometryPlots],
M3ConjugateBranchFiberFigure[N[lam],N[n],records,geometryPlots] },Alignment->Center,Spacings->1.0] ],

"Multiplicity / Gram / Cartan", Module[{}, weightMapData=M3WeightMultiplicityData[finite["LambdaN"],maxFrames,gridN,includeConjugate,rankTol,True];
Column[{ M3MultiplicitySummaryArticle[finite,iso,multDiag,gramDiag,records],

Grid[ { { If[AssociationQ[weightMapData], M3GTWeightDiagramArticle[weightMapData,finite["nN"]],
M3FailurePanel[weightMapData,"The weight diagram could not be constructed"] ], M3NormalizedSingularSpectrumArticle[C,rankTol] }, {
M3RankConvergenceArticle[C,multDiag["MultiplicityGT"],rankTol], If[AssociationQ[weightMapData],
M3MultiplicityDeficitMapArticle[weightMapData,finite["nN"]], M3FailurePanel[weightMapData,"The deficit map could not be constructed"] ] } },
Alignment->Center, Spacings->{.65,.65} ],

OpenerView[ { Style["Fidelity",10.5,Bold,FontFamily->"Times",FontColor->RGBColor[.16,.30,.47]], Column[{
Style["",8.8,Italic,GrayLevel[.38],FontFamily->"Times"], M3GramHeatmap[G], M3GTTable[gt] },Alignment->Center,Spacings->.5] }, False ] },Spacings->.75] ],

"Wigner-D / isospin", Column[{ M3WignerSummaryArticle[finite,iso], M3WignerFormulaPanel[], Grid[{{ M3WignerAngularResponseFigure[finite,iso],
M3WignerOrbitalProbabilityHeatmap[finite,iso] }},Alignment->Top,Spacings->{.75,.25}], M3OrbitalCoordinatesFramesFigure[records,iso,rankTol]
},Spacings->.70],

"Saddle / Hessian", Column[ { M3HessianExactContoursFigure[N[lam],N[n],finite,Q,rankTol], M3HessianExactVsQuadraticFigure[N[lam],N[n],finite,Q,rankTol] },
Spacings->1.0, Alignment->Center ],

"\[Sigma] / probabilities", Module[{sigmaScale}, phys=M3SelectedPhysicsData[N[lam],N[n],finite,iso,rec,rankTol]; If[FailureQ[phys]||!AssociationQ[phys],
M3FailurePanel[phys,"The physical density sigma could not be constructed"],
sigmaScale=Max[10^-12,Max[Abs[Flatten[Re[N[phys["Sigma"]]]]]],Max[Abs[Flatten[Im[N[phys["Sigma"]]]]]]]; Column[{ Style[
"Exact physical density and isospin probabilities", 12.5, Bold, FontFamily->"Times" ], Grid[{{
M3SigmaComponentHeatmap[phys["Sigma"],iso["JValues"],"Re",sigmaScale], M3SigmaComponentHeatmap[phys["Sigma"],iso["JValues"],"Im",sigmaScale]
}},Alignment->Top,Spacings->{.65,.15}], M3IsospinProbabilityFigure[iso["JValues"],phys["POrb"],phys["PPhys"]],
M3SigmaSpectrumFigure[phys["Sigma"],phys["SigmaDiagnostics"]], OpenerView[{ Style["Gram",10.5,Bold,FontFamily->"Times",FontColor->RGBColor[.12,.27,.45]],
Column[{ Grid[{{ M3SigmaHeatmapFigure[phys["Sigma"],iso["JValues"]], M3ScientificMatrixPlot[phys["RhoGamma"]["Gamma"],"|Gamma_ij|",Range[3]]
}},Alignment->Top,Spacings->{.65,.15}] },Spacings->.45] },False] },Spacings->.75] ] ],

"Detection", Module[{}, phys=M3SelectedPhysicsData[N[lam],N[n],finite,iso,rec,rankTol];
If[FailureQ[phys]||!AssociationQ[phys],Return[M3FailurePanel[phys,"sigma could not be constructed"]]]; U=M3Interferometer[interferometer,theta,phase];
transitions=M3AllTransitions[finite,iso,U,rankTol];
If[FailureQ[transitions],Return[M3FailurePanel[transitions,"The transitions could not be constructed"]]];
detPhys=M3PhysicalOutputDistribution[transitions,phys["Sigma"],rankTol]; detOrb=M3OrbitalOutputDistribution[transitions,orbitCoords];
If[FailureQ[detPhys]||FailureQ[detOrb],Return[M3FailurePanel[Failure["DetectionFailure",<|"Message"->"The output distributions could not be computed."|>],"Detection computation failed"]]];
comparison=<|"Physical"->detPhys,"Orbital"->detOrb|>; Column[{ Panel[ Row[{ Style[ "Final detection probability", 12.5,Bold,FontFamily->"Times" ],
Spacer[30], Style[ "The pure orbital-state distribution is compared with that of the physical density for the same interferometer.",
9.5,Bold,RGBColor[.10,.22,.58],FontFamily->"Times" ] }], Background->White ], Panel[
Grid[{{Style["Interferometer",Bold],interferometer,Style["theta",Bold],NumberForm[N[theta],{5,3}],Style["phase",Bold],NumberForm[N[phase],{5,3}]}}],
Background->White,FrameMargins->{{9,9},{6,6}},FrameStyle->Directive[GrayLevel[.86],AbsoluteThickness[.6]] ], M3DetectionComparisonFigure[comparison,18]
},Spacings->.8] ],

"Toroidal map", Module[{refData}, U=M3Interferometer[interferometer,theta,phase];
refData=M3ReferenceTorusMapData[N[lam],N[n],finite,rec,U,referenceSamples,rankTol]; If[FailureQ[refData]||!AssociationQ[refData],
M3FailurePanel[refData,"The reference toroidal map could not be constructed"], M3ReferenceTorusDashboard[refData] ] ],

"IRR decomposition / probabilities", Module[{rgRef}, U=M3Interferometer[interferometer,theta,phase];
rgRef=M3RhoGammaFromFrame[N[lam],N[n],rec["V"],rankTol]; If[FailureQ[rgRef]||!AssociationQ[rgRef],
M3FailurePanel[rgRef,"Gamma could not be reconstructed for the IRR explorer"], M3ReferenceIRRExplorer[ finite["nN"], rgRef["Gamma"], U, finite["LambdaN"],
rankTol ] ] ],

"orthonormal basis", M3QBExplorer[ N[lam],N[n],finite,iso,records,rankTol ],

_,M3CaseSummary[N[lam],N[n],finite,iso,records,gramDiag] ];

Pane[ Column[{header,tabbar,content},Spacings->.65], {1080,800},Scrollbars->{False,True},Alignment->{Left,Top} ] ],

Delimiter, Style[Row[{"Spectrum ",TraditionalForm[\[Lambda]]}],Bold,FontFamily->"Times"],
{{lamA,6,Style[TraditionalForm[Subscript[\[Lambda],1]],11]},1,12,1,Appearance->"Labeled"},
{{lamB,3,Style[TraditionalForm[Subscript[\[Lambda],2]],11]},1,12,1,Appearance->"Labeled"},
{{lamC,1,Style[TraditionalForm[Subscript[\[Lambda],3]],11]},1,12,1,Appearance->"Labeled"},

Delimiter, Style[Row[{"Occupation ",TraditionalForm["n"]}],Bold,FontFamily->"Times"],
{{nA,1,Style[TraditionalForm[Subscript["n",1]],11]},1,12,1,Appearance->"Labeled"},
{{nB,1,Style[TraditionalForm[Subscript["n",2]],11]},1,12,1,Appearance->"Labeled"},
{{nC,1,Style[TraditionalForm[Subscript["n",3]],11]},1,12,1,Appearance->"Labeled"},

Delimiter, Style["Finite representation",Bold,FontFamily->"Times"], {{scaleFactor,1,Row[{"N / ",Subscript["N","min"]}]},1,4,1,Appearance->"Labeled"},

Delimiter, Style["Frame bank",Bold,FontFamily->"Times"], {{maxFrames,12,"maximum frames"},4,40,1,Appearance->"Labeled"},
{{gridN,20,"fiber resolution"},10,60,1,Appearance->"Labeled"}, {{includeConjugate,True,"conjugate branches"},{False,True}},
{{rankTol,10^-8,"rank tolerance"},{10^-6->"10^-6",10^-7->"10^-7",10^-8->"10^-8",10^-9->"10^-9",10^-10->"10^-10"},ControlType->PopupMenu},

Delimiter, Style["Interferometer",Bold,FontFamily->"Times"], {{interferometer,"Fourier F3","Output U"},{"Fourier F3","Identity","Beam splitter 1-2"}},
{{theta,Pi/4,TraditionalForm[\[Theta]]},0,Pi/2,Pi/32,Appearance->"Labeled"}, {{phase,0.,TraditionalForm[\[Phi]]},0,2 Pi,Pi/16,Appearance->"Labeled"},

Delimiter, Style["Reference panels",Bold,FontFamily->"Times"], {{referenceSamples,2500,"toroidal-map samples"},500,6000,500,Appearance->"Labeled"},

{{selectedFrame,1},ControlType->None}, {{activeTab,"Summary"},ControlType->None},

ControlPlacement->Left, ContinuousAction->False, SynchronousUpdating->False, SaveDefinitions->False ];


(* Launch the interactive SU(3) analysis interface. *)
M3SU3Interface[]
