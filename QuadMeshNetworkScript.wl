(* ::Package:: *)

(* ::Input::Initialization:: *)
<<NDSolve`FEM`;


(*Note: the preset[] function is required to run before this to build the basic topology/geometry of the network manifold, providing variables such as edgeList, vertexPolyFacesByLines, etc.*)


(* ::Input::Initialization:: *)
(*setupFEM[globalScale]: This function creates the necessary symoblic variables for Mathematica's use.*)
(*globalScale - determines the mesh's resolution. Smaller scale = higher resolution.*)
setupFEM[globalScale_:1.0]:=Module[{},
\[CapitalDelta]w=0.12566370614359174`*globalScale;
\[CapitalDelta]l=\[CapitalDelta]w;
Do[
Subscript[edgeLength, edgeN]=edgeLengthList[[edgeN]];
Subscript[edgeWidth, edgeN]=edgeWidthList[[edgeN]];
Subscript[width, edgeN]=Round[Subscript[edgeWidth, edgeN]/\[CapitalDelta]w];
Subscript[radius, edgeN]=Subscript[edgeWidth, edgeN]/Perimeter[RegularPolygon[Subscript[width, edgeN]]];
Subscript[length, edgeN]=Round[Subscript[edgeLength, edgeN]/\[CapitalDelta]l];

Subscript[meshRegion, edgeN]=MeshRegion[ToElementMesh[Rectangle[{0,1},{Subscript[width, edgeN],Subscript[length, edgeN]}],"MeshElementType"->QuadElement,MaxCellMeasure->1,"MeshOrder"->1]];
(*sqrtZScale: allow the length ("z" scale) of each tube to be adjusted during optimization.*)
Subscript[meshCoordinates, edgeN]={Subscript[radius, edgeN] Cos[2\[Pi] #[[1]]/Subscript[width, edgeN]],Subscript[radius, edgeN] Sin[2\[Pi] #[[1]]/Subscript[width, edgeN]],(Subscript[sqrtZScale, edgeN])^2*\[CapitalDelta]l*#[[2]]}&/@MeshCoordinates[Subscript[meshRegion, edgeN]][[;;Subscript[length, edgeN]*Subscript[width, edgeN]]];
Subscript[meshQuadCells, edgeN]=MeshCells[Subscript[meshRegion, edgeN],2]/.Table[Subscript[length, edgeN]*Subscript[width, edgeN]+i->i,{i,Subscript[length, edgeN]}];
Subscript[meshRegion, edgeN]=MeshRegion[Subscript[meshCoordinates, edgeN]/.{(Subscript[sqrtZScale, edgeN])^2->1},Subscript[meshQuadCells, edgeN]];
Subscript[meshLineCells, edgeN]=MeshCells[Subscript[meshRegion, edgeN],1];
Subscript[meshLabels, edgeN]=Range[Subscript[length, edgeN]*Subscript[width, edgeN]]+If[edgeN==1,0,Last[Subscript[meshLabels, edgeN-1]]];
Subscript[boundaryMeshLabelsS, edgeN]=Table[1+(w-1) Subscript[length, edgeN],{w,Subscript[width, edgeN]}]+(First[Subscript[meshLabels, edgeN]]-1);
Subscript[boundaryMeshLabelsT, edgeN]=Table[w Subscript[length, edgeN],{w,Subscript[width, edgeN]}]+(First[Subscript[meshLabels, edgeN]]-1);
,{edgeN,Length@edgeList}
];
meshCoordinatesAll=Join@@Table[Subscript[meshCoordinates, edgeN],{edgeN,Length@edgeList}];
];


(* ::Input::Initialization:: *)
(*doAttribute[]: This function creates the necessary geometric attributes of the network, e.g., how each vertex polygon looks like.*)
doAttribute[]:=Module[{},
(*Global attributes.*)
(*Use: vertexPolyConnections[[vertex]]=vertex connected by face {1,2,3,...}. The normal direction of the face is given by the sign: S (+), T(-), or no end (0).*)
vertexPolyConnections=Table[Table[Sign[#]edgeList[[Abs[#],-Sign[#]]]&[vertexPolyConnectionsByEdgeN[[v]][[f]]],{f,Length@vertexPolyConnectionsByEdgeN[[v]]}],{v,Length@vertexList}];
(*Use:.*)
vertexPolyFacesByCoordinates=Table[Table[vertexPolyLines[[v]][[Abs[#],Sign[#]]]&/@vertexPolyFacesByLines[[v]][[f]],{f,Length@vertexPolyFacesByLines[[v]]}],{v,Length@vertexList}];
(*Use: Displacement from S to T.*)
edgeDisplacementList=Table[Last@edgePathList[[edgeN]]-First@edgePathList[[edgeN]],{edgeN,Length@edgeList}];

(*Local attributes.*)
(*Use:.*)
vertexPolyFaceSideLengths=Table[
Module[{edgeNSigned=vertexPolyConnectionsByEdgeN[[v]][[faceN]],faceSideLengths,k=-1},
faceSideLengths=Round[\[CapitalDelta]w^-1 vertexPolyLineLengths[[v]][[Abs[vertexPolyFacesByLines[[v]][[faceN]]]]]];
While[faceSideLengths[[k]]==0,k--];
faceSideLengths[[k]]=faceSideLengths[[k]]+Subscript[width, Abs[edgeNSigned]]-Total[faceSideLengths];
faceSideLengths],{v,Length@vertexList},{faceN,Length@vertexPolyFacesByLines[[v]]}];
(*Use: Path tracking. From edgePathList to edgePathInterpolationList. t=1 to t=Subscript[length, edgeN].*)
edgePathInterpolationList=Table[
Module[{path=edgePathList[[edgeN]],intrinsicTime,leftPad,rightPad},
intrinsicTime={0}~Join~Accumulate[Norm/@Differences[path]];
(*Leave some space for vertices.*)
(*leftPad=Round[\[CapitalDelta]l^-1*Sqrt[(\[Pi](edgeWidthList[[edgeN]]/(2\[Pi]))^2*Length@vertexPolyConnections[[edgeList[[edgeN]][[1]]]])/(4\[Pi])]];
rightPad=Round[\[CapitalDelta]l^-1*Sqrt[(\[Pi](edgeWidthList[[edgeN]]/(2\[Pi]))^2*Length@vertexPolyConnections[[edgeList[[edgeN]][[2]]]])/(4\[Pi])]];*)
leftPad=(Subscript[length, edgeN]-1)/Last[intrinsicTime] Round[0.2Subscript[edgeWidth, edgeN]+Norm[Mean@vertexPolyCoordinates[[edgeList[[edgeN]][[1]]]][[vertexPolyFacesByCoordinates[[#[[1]],#[[2]]]]&@FirstPosition[vertexPolyConnectionsByEdgeN,edgeN]]]-Mean@vertexPolyCoordinates[[edgeList[[edgeN]][[1]]]]]];
rightPad=(Subscript[length, edgeN]-1)/Last[intrinsicTime] Round[0.2Subscript[edgeWidth, edgeN]+Norm[Mean@vertexPolyCoordinates[[edgeList[[edgeN]][[2]]]][[vertexPolyFacesByCoordinates[[#[[1]],#[[2]]]]&@FirstPosition[vertexPolyConnectionsByEdgeN,-edgeN]]]-Mean@vertexPolyCoordinates[[edgeList[[edgeN]][[2]]]]]];
(*Print["pads: ",leftPad," ,",rightPad];*)
intrinsicTime=intrinsicTime/Last[intrinsicTime]*(Subscript[length, edgeN]-1+leftPad+rightPad)+1-leftPad;
Interpolation[Transpose[{Transpose[{intrinsicTime}],path}],InterpolationOrder->1]
],
{edgeN,Length@edgeList}];
(*Use: Rotation list. Accumulated. Starting from IdentityMatrix[3]. t=1 to t=Subscript[length, edgeN]. *)
edgePathInterpolationRotationList=ConstantArray[{},Length@edgeList];
Do[edgePathInterpolationRotationList[[edgeN]]=Module[{normalPre=Normalize[edgePathInterpolationList[[edgeN]][1]-edgePathInterpolationList[[edgeN]][0]],normalPost,rotate=N@IdentityMatrix[3]},{rotate}~Join~Table[normalPost=Normalize[edgePathInterpolationList[[edgeN]][t]-edgePathInterpolationList[[edgeN]][t-1]];
rotate=If[Chop[VectorAngle[normalPre,normalPost]]==0.,N@IdentityMatrix[3],RotationMatrix[{normalPre,normalPost}]] . rotate;
normalPre=normalPost;
rotate,{t,2,Subscript[length,edgeN],1/100 (*higher resolution*)}][[;;;;100]]];,{edgeN,Length@edgeList}];
(*Use: Boundary coordinates.*)
edgeBoundaryCoordinatesFatterSList=ConstantArray[{},{Length@edgeList}];
edgeBoundaryCoordinatesFatterTList=ConstantArray[{},{Length@edgeList}];
edgeBoundaryCoordinatesSList=ConstantArray[{},{Length@edgeList}];
edgeBoundaryCoordinatesTList=ConstantArray[{},{Length@edgeList}];
edgeBoundaryCoordinatesNormalSList=ConstantArray[{},{Length@edgeList}];
edgeBoundaryCoordinatesNormalTList=ConstantArray[{},{Length@edgeList}];
Module[{generateBoundaryCoordinates},
generateBoundaryCoordinates[edgeNSigned_]:=Module[{polyCoordinates,polyEdgeVectors,faceSideLengths,minFaceSideLengths,v,faceN,accumulate},
{v,faceN}=FirstPosition[vertexPolyConnectionsByEdgeN,edgeNSigned];
polyCoordinates=vertexPolyCoordinates[[v]][[vertexPolyFacesByCoordinates[[v]][[faceN]]]];
faceSideLengths=vertexPolyFaceSideLengths[[v]][[faceN]];
(*Since some boundary labels are dangling (thanks to the "Min" in cposGlueAll). we cannot count them. We have to recalculate the faceSideLengths.*)
minFaceSideLengths=Table[
Min@With[{
p1=Position[vertexPolyFacesByLines[[v]],vertexPolyFacesByLines[[v]][[faceN]][[e]]],
p2=Position[vertexPolyFacesByLines[[v]],-vertexPolyFacesByLines[[v]][[faceN]][[e]]]
},
vertexPolyFaceSideLengths[[v]][[#[[1]],#[[2]]]]&/@Join[p1,p2]
],{e,Length@vertexPolyFacesByLines[[v]][[faceN]]}];
If[edgeNSigned<0,polyCoordinates=Reverse@RotateLeft[polyCoordinates,1];faceSideLengths=Reverse@faceSideLengths;minFaceSideLengths=Reverse@minFaceSideLengths];
polyEdgeVectors=RotateLeft[polyCoordinates]-polyCoordinates;
accumulate=Accumulate[faceSideLengths]~Prepend~0;
Print[accumulate];
(*The following is to match the cposGlueAll, making both cposGlueAll and cposCoordinateS/T to be zero during minimization.*)
Table[
Most[If[vertexPolyFacesByLines[[v]][[faceN]][[e]]>0,
ArrayPad[#,{{0,faceSideLengths[[e]]-minFaceSideLengths[[e]]}},"Fixed"],
ArrayPad[#,{{faceSideLengths[[e]]-minFaceSideLengths[[e]],0}},"Fixed"]
]]&@If[minFaceSideLengths[[e]]==0,{polyCoordinates[[e]]},ArrayResample[{polyCoordinates[[e]],polyCoordinates[[e]]+polyEdgeVectors[[e]]},{minFaceSideLengths[[e]]+1 (*To cut the total length into n pieces one needs to resample the array by n+1*),3}]],
{e,Length@polyEdgeVectors}]
];
Do[
(*Generate boundary coordinates from vertex polyhedrons.*)
edgeBoundaryCoordinatesSList[[edgeN]]=generateBoundaryCoordinates[edgeN];
edgeBoundaryCoordinatesTList[[edgeN]]=generateBoundaryCoordinates[-edgeN];
(*Alternatively, copy boundary coordinates from previous iterations.*)
(*edgeBoundaryCoordinatesSList[[edgeN]]={Flatten[initialCoordinatesList,1][[Subscript[boundaryMeshLabelsS, edgeN]]]};
edgeBoundaryCoordinatesTList[[edgeN]]={Flatten[initialCoordinatesList,1][[Subscript[boundaryMeshLabelsT, edgeN]]]};*)
(*Alternatively, self-define boundary coordinates.*)
(*edgeBoundaryCoordinatesSList[[1]]=TakeList[(Subscript[meshCoordinates, 1][[Table[1+(w-1) Subscript[length, 1],{w,Subscript[width, 1]}]]]/.Subscript[sqrtZScale, 1]->1),vertexPolyFaceSideLengths[[1]][[1]]];
edgeBoundaryCoordinatesTList[[1]]=TakeList[{#[[1]],#[[2]],#[[3]]/geodesicScale}&/@(Subscript[meshCoordinates, 1][[Table[w Subscript[length, 1],{w,Subscript[width, 1]}]]]/.Subscript[sqrtZScale, 1]->1),vertexPolyFaceSideLengths[[2]][[1]]];*)

(*Normal directions.*)
edgeBoundaryCoordinatesNormalSList[[edgeN]]=If[Length@vertexPolyConnections[[edgeList[[edgeN]][[1]]]]<=3,Normalize@edgeDisplacementList[[edgeN]],Normalize[(#[[1]]-Mean[#])\[Cross](#[[2]]-Mean[#])]&@Flatten[edgeBoundaryCoordinatesSList[[edgeN]],1]];
edgeBoundaryCoordinatesNormalTList[[edgeN]]=If[Length@vertexPolyConnections[[edgeList[[edgeN]][[2]]]]<=3,Normalize@edgeDisplacementList[[edgeN]],Normalize[(#[[1]]-Mean[#])\[Cross](#[[2]]-Mean[#])]&@Flatten[edgeBoundaryCoordinatesTList[[edgeN]],1]];

(*Small expansions added to boundary coordinates to avoid degeneracy.*)
edgeBoundaryCoordinatesFatterSList[[edgeN]]=edgeBoundaryCoordinatesSList[[edgeN]]+TakeList[0.1Subscript[edgeWidth, edgeN]Map[Normalize[#\[Cross]edgeBoundaryCoordinatesNormalSList[[edgeN]]]&,(#-RotateRight[#,1])&[Flatten[edgeBoundaryCoordinatesSList[[edgeN]],1]]],Length/@edgeBoundaryCoordinatesSList[[edgeN]]];
edgeBoundaryCoordinatesFatterTList[[edgeN]]=edgeBoundaryCoordinatesTList[[edgeN]]+TakeList[0.1Subscript[edgeWidth, edgeN]Map[Normalize[#\[Cross]edgeBoundaryCoordinatesNormalTList[[edgeN]]]&,(#-RotateRight[#,1])&[Flatten[edgeBoundaryCoordinatesTList[[edgeN]],1]]],Length/@edgeBoundaryCoordinatesTList[[edgeN]]];

,{edgeN,Length@edgeList}];
];

(*Use: Periodic boundary condition (PBC). We use the fact that the PBC was reflected by Mean@vertexPolyCoordinates[[edgeList[[edgeN]][[2]]]] \[NotEqual] Last@edgePathList[[edgeN]] in the preset.*)
edgePBCLabelsList=Table[Flatten[Table[Thread[{Subscript[x, label],Subscript[y, label],Subscript[z, label]}->{Subscript[x, label],Subscript[y, label],Subscript[z, label]}+Mean@vertexPolyCoordinates[[edgeList[[edgeN]][[2]]]]-Last@edgePathList[[edgeN]]],{label,Subscript[boundaryMeshLabelsT, edgeN]~Join~(Subscript[boundaryMeshLabelsT, edgeN]-1)}]],{edgeN,Length@edgeList}];
];


(* ::Input::Initialization:: *)
(*doConstraint[]: Set up the cost functions and constraints.*)
doConstraint[]:=Module[{},
(*Constraint set 1. Intrinsic constraints.*)
(*1a. Isometry parallelogram.*)
Do[
{Subscript[cisoA, edgeN],Subscript[cisoB, edgeN],Subscript[cisoC, edgeN]}=Transpose@ParallelTable[Simplify@{
(1+(Subscript[d\[Lambda], ToString@f])^2)Total[(meshCoordinatesAll[[f[[1]]]]-meshCoordinatesAll[[f[[3]]]])^2]-Total[{Subscript[x, f[[1]]]-Subscript[x, f[[3]]],Subscript[y, f[[1]]]-Subscript[y, f[[3]]],Subscript[z, f[[1]]]-Subscript[z, f[[3]]]}^2],
(1+(Subscript[d\[Lambda], ToString@f])^2)Total[(meshCoordinatesAll[[f[[2]]]]-meshCoordinatesAll[[f[[4]]]])^2]-Total[{Subscript[x, f[[2]]]-Subscript[x, f[[4]]],Subscript[y, f[[2]]]-Subscript[y, f[[4]]],Subscript[z, f[[2]]]-Subscript[z, f[[4]]]}^2],
(1+(Subscript[d\[Lambda], ToString@f])^2)(meshCoordinatesAll[[f[[1]]]]-meshCoordinatesAll[[f[[3]]]]) . (meshCoordinatesAll[[f[[2]]]]-meshCoordinatesAll[[f[[4]]]])-{Subscript[x, f[[1]]]-Subscript[x, f[[3]]],Subscript[y, f[[1]]]-Subscript[y, f[[3]]],Subscript[z, f[[1]]]-Subscript[z, f[[3]]]} . {Subscript[x, f[[2]]]-Subscript[x, f[[4]]],Subscript[y, f[[2]]]-Subscript[y, f[[4]]],Subscript[z, f[[2]]]-Subscript[z, f[[4]]]}
},{f,(First/@Subscript[meshQuadCells, edgeN])/.(label_Integer:>label+(First[Subscript[meshLabels, edgeN]]-1))}];,
{edgeN,Length@edgeList}
];

Do[Subscript[constraintsAllciso, edgeN]=Total[(Subscript[cisoA, edgeN])^2]+Total[(Subscript[cisoB, edgeN])^2]+Total[(Subscript[cisoC, edgeN])^2],{edgeN,Length@edgeList}];

(*1b. Conformal-weight scalar factors.*)
Do[
Subscript[cconf\[Lambda], edgeN]=Table[(Subscript[d\[Lambda], ToString@f])^2,{f,(First/@Subscript[meshQuadCells, edgeN])/.(label_Integer:>label+(First[Subscript[meshLabels, edgeN]]-1))}];,
{edgeN,Length@edgeList}
];

(*1c. (Soft) constraint for gluing: the boundary coordinates are glued.*)
(*Note: it is possible that the lengths of the boundary coordinates do not match each other (because of the use of Round[]). We simply let the mismatched linger around without being glued, by using Min[] for the lengths.*)
boundaryConstraintsG[meshLabels1_,meshLabels2_]:=Module[
{meshCoordinates1=Table[{Subscript[x, meshLabels1[[i]]],Subscript[y, meshLabels1[[i]]],Subscript[z, meshLabels1[[i]]]},{i,Min[Length@meshLabels1,Length@meshLabels2]}],meshCoordinates2=Table[{Subscript[x, meshLabels2[[i]]],Subscript[y, meshLabels2[[i]]],Subscript[z, meshLabels2[[i]]]},{i,Min[Length@meshLabels1,Length@meshLabels2]}]},
Flatten[meshCoordinates1-meshCoordinates2]
];
findLabelPartitions[faceSideLengths_,twist_:0]:=Module[{accumulate=Accumulate[faceSideLengths]~Prepend~0,range},
range=RotateLeft[Range[accumulate[[-1]]]~Join~Range[accumulate[[-1]]]~Join~Range[accumulate[[-1]]],twist];
Table[range[[(accumulate[[-1]])+accumulate[[r]]+1;;(accumulate[[-1]]+1)+accumulate[[r+1]]]],{r,Length[accumulate]-1}]];
(*Check.*)(*findLabelPartitions[{11,17,22},10]*)
findBoundaryMeshLabelsPerFaceSide[v_,faceN_,faceSideNSigned_]:=Module[{edgeNSigned,boundaryMeshLabels,labelPartitions,faceSideLengths,source,target,k=-1},
edgeNSigned=vertexPolyConnectionsByEdgeN[[v]][[faceN]];
If[edgeNSigned==0,Return[Nothing];];
boundaryMeshLabels=If[edgeNSigned>0,Subscript[boundaryMeshLabelsS, edgeNSigned],Reverse[RotateLeft[Subscript[boundaryMeshLabelsT, -edgeNSigned],1]]];
labelPartitions=findLabelPartitions[vertexPolyFaceSideLengths[[v]][[faceN]],If[edgeNSigned>0,(*S*)0,(*T*)Round[Subscript[width, -edgeNSigned] edgeTwistList[[-edgeNSigned]]/(2\[Pi])]]];
(*Symmetric twist.*)(*labelPartitions=findLabelPartitions[vertexPolyFaceSideLengths[[v]][[faceN]],If[edgeNSigned>0,(*S*)-Round[Subscript[width, edgeNSigned](edgeTwistList[[edgeNSigned]]/(2\[Pi]))/2],(*T*)Round[Subscript[width, -edgeNSigned](edgeTwistList[[-edgeNSigned]]/(2\[Pi]))/2]]];*)
boundaryMeshLabels[[If[faceSideNSigned>0,labelPartitions[[faceSideNSigned]],Reverse[labelPartitions[[-faceSideNSigned]]]]]]
];
cposGlueAll=Flatten@Table[
Module[{faceNList,faceSideNSignedList},
Table[
{faceNList,faceSideNSignedList}=With[{
p1=Position[vertexPolyFacesByLines[[v]],polyLineN],
p2=Position[vertexPolyFacesByLines[[v]],-polyLineN]
},
{p1[[;;,1]]~Join~p2[[;;,1]],p1[[;;,2]]~Join~-p2[[;;,2]]}
];
(*Two and only two polyLines will be glued!*)
If[Length[faceNList]<2,Nothing,boundaryConstraintsG[findBoundaryMeshLabelsPerFaceSide[v,faceNList[[1]],faceSideNSignedList[[1]]],findBoundaryMeshLabelsPerFaceSide[v,faceNList[[2]],faceSideNSignedList[[2]]]]]
,{polyLineN,Length@vertexPolyLines[[v]]}]
],
{v,Length@vertexList}];

(*1d. Fairness + boundary fairness.*)
Module[{triples},
Do[
triples=With[{points=Partition[Range[Subscript[length, edgeN]*Subscript[width, edgeN]],Subscript[length, edgeN]]},
Flatten[Partition[#,3,1]&/@points,1]~Join~Flatten[(Partition[#,3,1,1]&/@Transpose[points])[[2;;-2]],1]
(*Why [[2;;-2]]? Dispose two boundaries.*)
]+(First[Subscript[meshLabels, edgeN]]-1);
(*Extrinsic-only.*)
(*cfair=Table[With[{\[CapitalDelta]=normCoordinates[[triple[[1]]]]-2normCoordinates[[triple[[2]]]]+normCoordinates[[triple[[3]]]],n={Subscript[nx, triple[[2]]],Subscript[ny, triple[[2]]],Subscript[nz, triple[[2]]]}/.findMinimum[[2]]},\[CapitalDelta]-(n.\[CapitalDelta])n],{triple,triples}]*)
(*Scale-relevant.*)
(*Subscript[cfair, edgeN]=Flatten@Table[{Subscript[x, triple[[1]]],Subscript[y, triple[[1]]],Subscript[z, triple[[1]]]}-2{Subscript[x, triple[[2]]],Subscript[y, triple[[2]]],Subscript[z, triple[[2]]]}+{Subscript[x, triple[[3]]],Subscript[y, triple[[3]]],Subscript[z, triple[[3]]]},{triple,triples}];*)
(*Scale-free.*)
Subscript[cfair, edgeN]=Flatten@Table[({Subscript[x, triple[[1]]],Subscript[y, triple[[1]]],Subscript[z, triple[[1]]]}-2{Subscript[x, triple[[2]]],Subscript[y, triple[[2]]],Subscript[z, triple[[2]]]}+{Subscript[x, triple[[3]]],Subscript[y, triple[[3]]],Subscript[z, triple[[3]]]})/(\[Sqrt](Total[({Subscript[x, triple[[1]]],Subscript[y, triple[[1]]],Subscript[z, triple[[1]]]}-{Subscript[x, triple[[2]]],Subscript[y, triple[[2]]],Subscript[z, triple[[2]]]})^2]+Total[({Subscript[x, triple[[2]]],Subscript[y, triple[[2]]],Subscript[z, triple[[2]]]}-{Subscript[x, triple[[3]]],Subscript[y, triple[[3]]],Subscript[z, triple[[3]]]})^2])),{triple,triples}];,
{edgeN,Length@edgeList}
];

triples={
If[MemberQ[Flatten[Table[Subscript[boundaryMeshLabelsS, edgeN],{edgeN,Length@edgeList}]],#[[1]]],#[[1]]+1,#[[1]]-1],
#[[1]],
If[MemberQ[Flatten[Table[Subscript[boundaryMeshLabelsS, edgeN],{edgeN,Length@edgeList}]],#[[2]]],#[[2]]+1,#[[2]]-1]
}&/@Abs[Transpose[{cposGlueAll[[;;;;3,1]],cposGlueAll[[;;;;3,2]]}]/.Subscript[(c_), label_]->label];
cfairb=Flatten@Table[({Subscript[x, triple[[1]]],Subscript[y, triple[[1]]],Subscript[z, triple[[1]]]}-2{Subscript[x, triple[[2]]],Subscript[y, triple[[2]]],Subscript[z, triple[[2]]]}+{Subscript[x, triple[[3]]],Subscript[y, triple[[3]]],Subscript[z, triple[[3]]]})/(\[Sqrt](Total[({Subscript[x, triple[[1]]],Subscript[y, triple[[1]]],Subscript[z, triple[[1]]]}-{Subscript[x, triple[[2]]],Subscript[y, triple[[2]]],Subscript[z, triple[[2]]]})^2]+Total[({Subscript[x, triple[[2]]],Subscript[y, triple[[2]]],Subscript[z, triple[[2]]]}-{Subscript[x, triple[[3]]],Subscript[y, triple[[3]]],Subscript[z, triple[[3]]]})^2])),{triple,triples}];
];
(*PBC.*)
cposGlueAll=cposGlueAll/.Association[Flatten@edgePBCLabelsList];
cfairb=cfairb/.Association[Flatten@edgePBCLabelsList];

(*Constraint set 2. Extrinsic constraints.*)
(*2a. Dirichlet/Neumann boundary conditions.*)
boundaryConstraintsD[boundaryCoordinates_,meshLabels_]:=With[
{meshCoordinates=Table[{Subscript[x, meshLabels[[i]]],Subscript[y, meshLabels[[i]]],Subscript[z, meshLabels[[i]]]},{i,Length@meshLabels}]},
Flatten[boundaryCoordinates-meshCoordinates]
];
boundaryConstraintsN[boundaryCoordinates_,meshLabels_]:=Module[
{meshCoordinates=Table[{Subscript[x, meshLabels[[i]]],Subscript[y, meshLabels[[i]]],Subscript[z, meshLabels[[i]]]},{i,Length@meshLabels}],meshVectors,boundaryVectors},
(*The total lengths and angles are fixed.*)
boundaryVectors=RotateLeft[boundaryCoordinates]-boundaryCoordinates;
meshVectors=RotateLeft[meshCoordinates]-meshCoordinates;
(*Individual length + absolute innerproduct.*)
{Table[Total[meshVectors[[i]]^2]-Total[boundaryVectors[[i]]^2],{i,Length@meshLabels}],
Table[meshVectors[[i]] . boundaryVectors[[i]]-Sqrt[Total[meshVectors[[i]]^2]Total[boundaryVectors[[i]]^2]],{i,Length@meshLabels}]}
];
Do[
Subscript[cposCoordinateS, edgeN]=boundaryConstraintsD[Flatten[edgeBoundaryCoordinatesSList[[edgeN]],1],Subscript[boundaryMeshLabelsS, edgeN]];
{Subscript[cposLengthS, edgeN],Subscript[cposAngleS, edgeN]}=boundaryConstraintsN[Flatten[edgeBoundaryCoordinatesSList[[edgeN]],1],Subscript[boundaryMeshLabelsS, edgeN]];
Subscript[cposCoordinateT, edgeN]=boundaryConstraintsD[RotateLeft[Flatten[edgeBoundaryCoordinatesTList[[edgeN]],1],Round[Subscript[width, edgeN] edgeTwistList[[edgeN]]/(2\[Pi])]],Subscript[boundaryMeshLabelsT, edgeN]];
{Subscript[cposLengthT, edgeN],Subscript[cposAngleT, edgeN]}=boundaryConstraintsN[RotateLeft[Flatten[edgeBoundaryCoordinatesTList[[edgeN]],1],Round[Subscript[width, edgeN] edgeTwistList[[edgeN]]/(2\[Pi])]],Subscript[boundaryMeshLabelsT, edgeN]];,
{edgeN,Length@edgeList}
];

(*2b. Follow-the-path constraints.*)
(*The mean of all mesh labels along a path is set to be equal to the mean of the different steps along the path.*)Do[Subscript[cposPathMean, edgeN]=Flatten@Table[Sum[{Subscript[x, meshLabel],Subscript[y, meshLabel],Subscript[z, meshLabel]},{meshLabel,(Table[t+(w-1) Subscript[length, edgeN],{w,Subscript[width, edgeN]}]+(First[Subscript[meshLabels, edgeN]]-1))}]/Subscript[width, edgeN]-edgePathInterpolationList[[edgeN]][t],{t,Subscript[length, edgeN]}];,
{edgeN,Length@edgeList}];

(*2c. Terminal constraints (weak). The mean of all boundary mesh coordinates at the terminal is set to be equal to the vertex's center, i.e., the mean of the corresponding vertexPolyCoordinates. This is only for vertices that are terminals.*)
cposTerminalMean=Flatten@Table[
With[{boundaryMeshLabels=If[pos[[2]]==1,Subscript[boundaryMeshLabelsS, pos[[1]]],Subscript[boundaryMeshLabelsT, pos[[1]]]]},
Sum[{Subscript[x, boundaryMeshLabel],Subscript[y, boundaryMeshLabel],Subscript[z, boundaryMeshLabel]},{boundaryMeshLabel,boundaryMeshLabels}]/Subscript[width, pos[[1]]]-Mean[vertexPolyCoordinates[[edgeList[[pos[[1]],pos[[2]]]]]]]
],
{pos,Flatten[Position[edgeList,#]&/@Pick[vertexList,vertexPolyIsTerminals,True],1]}];

(*2d. Terminal constraints (strong). All boundary mesh labels are pushed away from the vertex's center.*)
cposTerminalCircle=Flatten@Table[
With[{boundaryMeshLabels=If[pos[[2]]==1,Subscript[boundaryMeshLabelsS, pos[[1]]],Subscript[boundaryMeshLabelsT, pos[[1]]]]},
Sum[(Total[({Subscript[x, boundaryMeshLabel],Subscript[y, boundaryMeshLabel],Subscript[z, boundaryMeshLabel]}-Mean[vertexPolyCoordinates[[edgeList[[pos[[1]],pos[[2]]]]]]])^2]-(Subscript[radius, pos[[1]]])^2)^2,{boundaryMeshLabel,boundaryMeshLabels}]
],
{pos,Flatten[Position[edgeList,#]&/@Pick[vertexList,vertexPolyIsTerminals,True],1]}];

];


(* ::Input::Initialization:: *)
(*doInitialize[edged\[Lambda]List,edgeSqrtZScaleList]: Build the initial network.*)
(*edged\[Lambda]List - initial values of the onformal-weight scalar factor for each tube.*)
(*edgeSqrtZScaleList - initial values of the "z" scale for each tube.*)
doInitialize[edged\[Lambda]List_,edgeSqrtZScaleList_]:=Module[{},
parametersList=ConstantArray[{},Length@edgeList];
Do[
parametersList[[edgeN]]=Flatten@Table[{Subscript[x,i],Subscript[y,i],Subscript[z,i]},{i,Subscript[meshLabels, edgeN]}];
parametersList[[edgeN]]=parametersList[[edgeN]]~Join~Table[Subscript[d\[Lambda], ToString@f],{f,(First/@Subscript[meshQuadCells, edgeN])/.(label_Integer:>label+(First[Subscript[meshLabels, edgeN]]-1))}];
parametersList[[edgeN]]=parametersList[[edgeN]]~Join~{Subscript[sqrtZScale, edgeN]};,
{edgeN,Length@edgeList}
];

(*Coordinates initialization.*)
initialCoordinatesList=ConstantArray[{},Length@edgeList];
Do[
(*Normal direction extension: along edge path. Refined.*)
initialCoordinatesList[[edgeN]]=Module[
{coordinatesS=Flatten[edgeBoundaryCoordinatesFatterSList[[edgeN]],1],coordinatesT=RotateLeft[Flatten[edgeBoundaryCoordinatesFatterTList[[edgeN]],1],Round[Subscript[width, edgeN] edgeTwistList[[edgeN]]/(2\[Pi])]],
coordinatesSFake,coordinatesTFake,coordinatesFake,
normalPre=Normalize[edgePathInterpolationList[[edgeN]][1]-edgePathInterpolationList[[edgeN]][0]],
normalPost=Normalize[edgePathInterpolationList[[edgeN]][Subscript[length, edgeN]]-edgePathInterpolationList[[edgeN]][Subscript[length, edgeN]-1]],
accumulate,rotate},
accumulate=Accumulate[{edgePathInterpolationList[[edgeN]][1]}~Join~Table[Norm[edgePathInterpolationList[[edgeN]][t]-edgePathInterpolationList[[edgeN]][t-1]]*normalPre,{t,2,Subscript[length, edgeN]}]];
rotate=If[Chop[VectorAngle[edgeBoundaryCoordinatesNormalSList[[edgeN]],normalPre]]==0.,N@IdentityMatrix[3],RotationMatrix[{edgeBoundaryCoordinatesNormalSList[[edgeN]],normalPre}]];
coordinatesSFake=Table[rotate . (coordinatesS[[i]]-Mean[coordinatesS])+accumulate[[1]],{i,Length@coordinatesS}];
rotate=If[Chop[VectorAngle[edgeBoundaryCoordinatesNormalTList[[edgeN]],normalPost]]==0.,N@IdentityMatrix[3],RotationMatrix[{edgeBoundaryCoordinatesNormalTList[[edgeN]],normalPost}]];
rotate=(edgePathInterpolationRotationList[[edgeN]][[Subscript[length, edgeN]]])\[ConjugateTranspose] . rotate;
coordinatesTFake=Table[rotate . (coordinatesT[[i]]-Mean[coordinatesT])+accumulate[[Subscript[length, edgeN]]],{i,Length@coordinatesT}];
coordinatesFake=Interpolation[{{{1},coordinatesSFake},{{Subscript[length, edgeN]},coordinatesTFake}},InterpolationOrder->1]/@Range[Subscript[length, edgeN]];
Flatten[Transpose[Table[
Table[edgePathInterpolationRotationList[[edgeN]][[t]] . (coordinatesFake[[t]][[i]]-Mean[coordinatesFake[[t]]])+edgePathInterpolationList[[edgeN]][t],{i,Length@coordinatesFake[[t]]}],
{t,Subscript[length, edgeN]}]],1]
];,
{edgeN,Length@edgeList}];


findMinimumList=Table[{10^9,Thread[parametersList[[edgeN]]->(Flatten@initialCoordinatesList[[edgeN]]~Join~edged\[Lambda]List[[edgeN]]~Join~edgeSqrtZScaleList[[edgeN]])]},{edgeN,Length@edgeList}];

globalFindMinimumFactor=(10^3 *Sum[Subscript[edgeLength, edgeN] Subscript[edgeWidth, edgeN],{edgeN,Length@edgeList}]^-1);
parametersAll=Flatten@parametersList;
findMinimumAll={globalFindMinimumFactor*10^9,Flatten@findMinimumList[[;;,2]]};

];


(* ::Input::Initialization:: *)
(*doChangeMode[conf\[Lambda],zScale,stillKeepd\[Lambda],stillKeepZScaleForEdges]: Some extra modifications.*)
doChangeMode[conf\[Lambda]_:True,zScale_:True,stillKeepd\[Lambda]_:False,stillKeepZScaleForEdges_:{}]:=Module[{},
(*No conf\[Lambda]: disallows d\[Lambda], for precise isometric immersion purposes.*)
If[conf\[Lambda]==False,
wconf\[Lambda]=\[Infinity];
With[{d\[Lambda]ReplaceAssocList=Table[If[stillKeepd\[Lambda],AssociationThread[Keys[#]->Sqrt[Chop[Values[#]^2]+Keys[#]^2]],#]&@Association[findMinimumList[[edgeN,2,3Subscript[width, edgeN]*Subscript[length, edgeN]+1;;-2]]],{edgeN,Length@edgeList}]},
Do[Subscript[constraintsAllciso, edgeN]=Chop[Subscript[constraintsAllciso, edgeN]/.d\[Lambda]ReplaceAssocList[[edgeN]]],{edgeN,Length@edgeList}];
];
With[{sqrt\[Lambda]ScaleList=Table[Values@KeySort@Merge[(Mod[Min@ToExpression[#[[1]][[2]]]-findMinimumList[[edgeN,2,1,1,2]]-1,Subscript[length, edgeN]]+1)->Chop[Sqrt[1+#[[2]]^2]]&/@findMinimumList[[edgeN,2,3Subscript[width, edgeN]*Subscript[length, edgeN]+1;;-2]],Mean],{edgeN,Length@edgeList}]},
Do[Subscript[cposPathMean, edgeN]=Flatten@Table[Sum[{Subscript[x, meshLabel],Subscript[y, meshLabel],Subscript[z, meshLabel]},{meshLabel,(Table[t+(w-1) Subscript[length, edgeN],{w,Subscript[width, edgeN]}]+(First[Subscript[meshLabels, edgeN]]-1))}]/Subscript[width, edgeN]-edgePathInterpolationList[[edgeN]][1+Total[sqrt\[Lambda]ScaleList[[edgeN]][[;;t-1]]]/Total[sqrt\[Lambda]ScaleList[[edgeN]]] (Subscript[length, edgeN]-1)],{t,Subscript[length, edgeN]}];,
{edgeN,Length@edgeList}];
];];
(*No Subscript[sqrtZScale, edgeN]: disallows sqrtZScale, for conformal deformation purposes.*)
If[zScale==False,
wconf\[Lambda]=0;
Do[Subscript[constraintsAllciso, edgeN]=Chop[Subscript[constraintsAllciso, edgeN]/.<|findMinimumList[[edgeN,2,-1]]|>],{edgeN,Complement[Range[Length@edgeList],stillKeepZScaleForEdges]}];
];
];


(* ::Input::Initialization:: *)
(*doFindMinimum[wfair0,maxIterations,wSurface,wCircle,wPathPerFair,wIso]: Minimization of the total cost function.*)
doFindMinimum[wfair0_,maxIterations_:5000,wSurface_:10^-2,wCircle_:10^3,wPathPerFair_:0(*10^0*),wIso_:1]:=Module[{},
wfair=wfair0;
Do[
constraintsAll=
10^3 Total[(cposGlueAll)^2]+10^3 Total[cposTerminalMean^2]+wCircle*Total[cposTerminalCircle^2]+(10^-2)*Total[edgeLengthList]*wfair*Total[cfairb^2]
+Sum[
wIso*Subscript[constraintsAllciso, edgeN]
+wfair*Total[(Subscript[cfair, edgeN])^2](*Penalty for fairness.*)
(*+10^3(Total[(Subscript[cposCoordinateS, edgeN])^2]+Total[(Subscript[cposCoordinateT, edgeN])^2])(*Boundary conditions.*)*)
(*+wconf\[Lambda] Total[(Subscript[cconf\[Lambda], edgeN])^2](*Penalty for conformal deformation.*)*)
(*+10^3(Total[(Subscript[cposLengthS, edgeN])^2]+Total[(Subscript[cposLengthT, edgeN])^2])+10^3Total[(Subscript[ccircle, edgeN])^2]+10^-2Total[(Subscript[cposCoordinatesMeanS, edgeN])^2]+10^1Total[(Subscript[cposCoordinatesMeanRotateT, edgeN])^2](*Extras.*)*)
+wSurface*Subscript[edgeLength, edgeN] Subscript[edgeWidth, edgeN]Mean[1+Subscript[cconf\[Lambda], edgeN]](Subscript[sqrtZScale, edgeN])^2(*Surface area minimizer.*)
+wPathPerFair*wfair*Total[(Subscript[cposPathMean, edgeN])^2](*Path mean.*)
,{edgeN,Length@edgeList}];

constraintsAll=globalFindMinimumFactor*constraintsAll;(*Normalization.*)

Print["wfair=",wfair];
Print[AbsoluteTiming[
findMinimumAll=FindMinimum[constraintsAll,Transpose[{parametersAll,findMinimumAll[[2,;;,2]]}],MaxIterations->maxIterations];
]];
wfair=wfair/10;,
1];

(*Check.*)
wfair=wfair0;
Print["iso: ",globalFindMinimumFactor*Sum[Subscript[constraintsAllciso, edgeN],{edgeN,Length@edgeList}]/.Association[findMinimumAll[[2]]]];
Print["fairness: ",globalFindMinimumFactor*Sum[wfair Total[(Subscript[cfair, edgeN])^2],{edgeN,Length@edgeList}]/.Association[findMinimumAll[[2]]]];
Print["boundary fairness: ",globalFindMinimumFactor*(10^-2*Total[edgeLengthList]*wfair Total[cfairb^2])/.Association[findMinimumAll[[2]]]];
Print["surface: ",globalFindMinimumFactor*Sum[wSurface Subscript[edgeLength, edgeN] Subscript[edgeWidth, edgeN]Mean[1+Subscript[cconf\[Lambda], edgeN]](Subscript[sqrtZScale, edgeN])^2,{edgeN,Length@edgeList}]/.Association[findMinimumAll[[2]]]];
Print["path: ",globalFindMinimumFactor*Sum[wPathPerFair*wfair*Total[(Subscript[cposPathMean, edgeN])^2],{edgeN,Length@edgeList}]/.Association[findMinimumAll[[2]]]];

(*Output.*)
findMinimumList=Transpose[{ConstantArray[findMinimumAll[[1]],Length@findMinimumList],TakeList[findMinimumAll[[2]],Length/@findMinimumList[[;;,2]]]}];
];
