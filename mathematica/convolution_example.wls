#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll
sigma = 0.3;
f[x_] := (HeavisideTheta[x+5]-0.5*HeavisideTheta[x-5])-0.5*HeavisideTheta[x-10]
g[x_] := (1/Sqrt[2 Pi sigma^2]) Exp[-0.5 x^2/sigma^2]
h[y] = Convolve[f[x], g[x], x, y]

Plot[f[x],{x,-12,12}]
Plot[g[x],{x,-2,2}]
Plot[h[y],{y,-12,12}]

F[w] = FourierTransform[f[x],x,w]
G[w] = FourierTransform[g[x],x,w]
H[w] = F[w] G[w]
h[y] = InverseFourierTransform[H[w],w,y]
Plot[h[y],{y,-12,12}]





(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
