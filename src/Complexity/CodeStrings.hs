{- |
Module   : CodeStrings
Description: Strings for generated minimal Haskell code--functions for Algorithmic (Kolmogorov) complexity
Copyright : (c) 2019-2023 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License   :

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project

Maintainer : Ward Wheeler <wheeler@amnh.org>
Stability : unstable
Portability: portable (I hope)

-}
module Complexity.CodeStrings where

import      Complexity.Constants

--In general, function names are letter number and internal variable are letters only
-- exception \aG, \bG, cG, dG for \aGraph etc since 10 \a functions
-- \fM for \fModel when 0-9 used as for graph ones above

--remove prelude stuff for production
--only for debugging
programStartString :: String
programStartString="module Main where\n" --import Debug.Trace\n"

recursivePrintString :: String
recursivePrintString="\
\p0 a b=\n\
\ if nU b then putStr a\n\
\ else\n\
\ let ((c,d):e)=b\n\
\ in p0(a++(show c++\" \"++show d++\"\\n\")) e\n"

programStartStringGraph :: String
programStartStringGraph="module Main where\n"++recursivePrintString

getSingletonEdgesString :: String
getSingletonEdgesString="\
\aG b c=\n\
\ if c==0 then []\n\
\ else (2*b,2*b+1):aG(b+1)(c-1)\n"

minimalTreesString :: String
minimalTreesString="\
\bG a c d=\n\
\ if (d<2)||(c==0) then []\n\
\ else [(a,a+1),(a,a+2)]++bG(a+3)(c-2)(d-1)\n"

fullTreeString :: String
fullTreeString="\
\cG a b d\n\
\ |d==0=[]\n\
\ |a=[(b+2,b),(b+2,b+1)]++cG False (b+2)(d-2)\n\
\ |otherwise=[(b+2,b),(b+2,b+1)]++cG False (b+2)(d-1)\n"

addEdgeString :: String
addEdgeString="\
\dG a ((e,f):(g,h):i) c=\n\
\ if c==0 then (e,f):(g,h):i\n\
\ else\n\
\ let j=a\n\
\     k=a+1\n\
\     l=i++[(j,k),(e,j),(j,f),(g,k),(k,h)]\n\
\ in dG(a+2) l (c-1)\n"

factorialString :: String
factorialString="\
\fM 0=1\n\
\fM 1=1\n\
\fM n=n*fM (n-1)\n"

--power :: Double -> Int -> Int -> Double
--power base exponent counter=
powerString :: String
powerString="\
\p1 b x c\n\
\ |x==0=1\n\
\ |x==c=1\n\
\ |otherwise=b*(p1 b x(c+1))\n"
{-
\p1 b x c=\n\
\ if x==0 then 1\n\
\ else if x==c then 1\n\
\ else b*(p1 b x (c+1))\n"
-}

--expE :: Double -> Int -> Int -> Double
--expE exponent iterations counter=
-- This needs alot of iterations to behave well--100 for sure
-- Is n^2 due to use of factorial could be made linear
-- by saving previous factorial value and just mulgtipolying
-- the use of p1 also makes n^2--same preocedure--need a new value to pass
expEString :: String
expEString="\
\e0 x i c=\n\
\ if c=="++ fixedPrecisionString ++" then 0\n\
\ else ((p1 x c 0)/(fM$fromIntegral c))+e0 x i (c+1)\n"

{-
-- \e x i c=exp x\n"  test with library
expEString="\
\e x i c=p(1+(x/(fromIntegral "++ fixedPrecisionString ++")))" ++ fixedPrecisionString ++ " 0\n"
-}

--logE :: Double -> Int -> Int -> Double -> Double
--logE value iterations counter curValue=
logEString :: String
logEString="\
\l0 v i c w=\n\
\ if c==i then 2*w\n\
\ else l0 v i (c+1) (w+((p1 ((v-1)/(v+1)) (1+2*c) 0)/(fromIntegral$1+2*c)))\n"

-- log2 :: Double -> Int-> Int -> Double -> Double
-- log2 value iterations blah bleh=(logE value iterations 0 0.0)/(logE 2.0 iterations 0 0.0)
-- the 'c' and 'w' params are not used but there for function as parame with logE
log2String :: String
log2String="\
\b0 v i c w=(l0 v i 0 0)/(l0 2 i 0 0)\n"

-- start with 0
-- \r0 :: Int -> a -> Int -> [a]\n\
replicateString :: String
replicateString="\
\r0 n v c=\n\
\ if n==c then []\n\
\ else v:r0 n v (c+1)\n"

--makeSimpleMatrix :: Int -> Double -> Double -> Int -> String -> [[Double]]
--makeSimpleMatrix size diag nondiag rowCounter lastElement=
-- \m :: Int -> Double -> Double -> Int -> String -> [[Double]]\n\
makeMatrixString :: String
makeMatrixString="\
\m0 s d n c e\n\
\ |s==c=[]\n\
\ |(e/=\"-\")||(c<s-1)=((r0 c n 0)++[d]++(r0(s-c-1)n 0)):m0 s d n(c+1)e\n\
\ |otherwise=[(r0 c n 0)++[0]]\n"
{-
\m s d n c e=\n\
\ if s==c then []\n\
\ else if (e/=\"-\")||(c<s-1) then [(r c n 0)++[d]++(r(s-c-1)n 0)]++m s d n(c+1) e\n\
\ else [(r c n 0)++[0]]\n"
-}

-- is the first arg element of second (list)
-- eS :: a -> [a] -> Bool
elemString :: String
elemString="\
\eS a b\n\
\ |nU b=False\n\
\ |a==a5 b=True\n\
\ |otherwise=eS a (a6 b)\n" 

-- is the first arg not an element of second (list)
-- eS :: a -> [a] -> Bool
notElemString :: String
notElemString="\
\nE a b\n\
\ |nU b=True\n\
\ |a==a5 b=False\n\
\ |otherwise=nE a (a6 b)\n" 

-- getRepeatedElements in list
getRepeatedElementsString :: String
getRepeatedElementsString="\
\rE l@(x:y)\n\
\ |nU l=[]\n\
\ |nU y=[]\n\
\ |eS x y=x:(rE y)\n\
\ |otherwise=rE y\n"

-- childrenParentsOfNode call with 2 nU lists
-- childrenParentsOfNodeString :: Int -> [(Int,Int)] -> [(Int,Int)] -> [(Int,Int)] -> ([(Int,Int)], [(Int,Int)])
childrenParentsOfNodeString :: String
childrenParentsOfNodeString="\
\cN n l c p=\n\
\  if nU l then (c,p)\n\
\  else\n\
\  let (e,u)=a5 l\n\
\  in\n\
\  if n==e then cN n (a6 l) (u:c) p\n\
\  else if n==u then cN n (a6 l) c (e:p)\n\
\  else cN n (a6 l) c p\n"

-- get display edges from full input edge list (n)
-- not checking for nU list
-- first arg is 0/1 bit list 1 per net node
-- uses binary strem fo detemine left right
displayEdgesString :: String
displayEdgesString="\
\dE b n=\n\
\ let v=gM iM n\n\
\     r=rE v\n\
\ in\n\
\ if nU r then n\n\
\ else\n\
\ let h=a5 r\n\
\     (c,p)=cN h n [] []\n\
\     (_,pL')=cN (a5 p) n [] []\n\
\     (_,pR')=cN (a3 p) n [] []\n\
\     (pL, pR)=if(a5 b)==0 then (pL',pR') else (pR',pL')\n\
\     gR=a5 pR\n\
\     (cR, _)=cN (a3 p) n [] []\n\
\     sR=a5$g2 (/=h) cR\n\
\     d=[(h,a5 c),(a5 p, h),(a3 p, h),(gR,a3 p),(a3 p,sR)]\n\
\     a=[(a5 p, a5 c),(gR, sR)]\n\
\     w=(g2(`nE` d)n)++a\n\
\ in dE (a6 b) w\n"

fmapString :: String
fmapString="\
\gM f l=\n\
\ if nU l then []\n\
\ else f (a5 l):gM f (a6 l)\n"

nullString :: String
nullString="\
\nU a=(a==[])\n"

matrix2StringString :: String
matrix2StringString="\
\s0 m=\n\
\ if nU m then \"\\n\"\n\
\ else (a4$gM (++ \" \")$gM show$a5 m)++\"\\n\"++s0 (a6 m)\n"

neymanUniformString :: String
neymanUniformString="\
\u0 n a i=((a-((((fromIntegral n)-1))*(((e0 (-1*a) i 0)-1))))/((fromIntegral n)*a),(a-1+(e0 (-1*a) i 0))/((fromIntegral n)*a))\n"

--p argument here is ignored by exponential
neymanExponentialString :: String
neymanExponentialString="\
\y0 n a p=((1+(a*(fromIntegral n)))/((a+1)*(fromIntegral n)),1/((a+1)*(fromIntegral n)))\n"

-- |this passes function of Neyman Uniform or NeymanExponential as argument
-- f is functino, n=alpabet size, a is branch distribution param, p is precision integer,
-- \c9 :: (Int -> Double -> Int -> (Double, Double)) -> Int -> Double -> Int -> String -> [[Double]]\n\
makeTCMBitsString :: String
makeTCMBitsString="\
\c9 f n a p s=\n\
\ let (d,e)=f n a p\n\
\ in m0 n (-1*b0 d p 0 0)(-1*b0 e p 0 0) 0 s\n"

fstString :: String
fstString="\
\hM(a,_)=a\n"

sndString :: String
sndString="\
\iM(_,b)=b\n"

-- |does the fmap/sum operation over multiple pii/pij for weights if [(1,1)]
-- same as simple Neyman functions
--makeNeymanMatrix :: (Double -> Int-> Int -> Double -> Double) -> Distribution -> Int -> DistributionParameter -> Int -> [(Double, Double)] -> [[Double]]
--makeNeymanMatrix logType distribution cbetSize cParam iterations modifiers
-- codes in log2 unlike general cde that passes log type
-- \n :: Int -> Int -> Double -> Int -> [(Double, Double)] -> String -> [[Double]]\n\
makeNeymanGeneralMatrixString :: String
makeNeymanGeneralMatrixString="\
\n0 r l a p o s=\n\
\ let k=gM(jM r l a p)o\n\
\     q=c3(+)0$gM hM k\n\
\     c=c3(+)0$gM iM k\n\
\     e=(-1)*b0 q p 0 0\n\
\     f=(-1)*b0 c p 0 0\n\
\ in m0 l e f 0 s\n"

-- inner part of Neym,ena general--get the pii/pij for weight and fractions
-- d=0 is "uniform" distribution
-- \j :: Int -> Int -> Double -> Int -> (Double, Double) -> (Double, Double)\n\
neymanGeneralWithKString :: String
neymanGeneralWithKString="\
\jM d r a i (w, f)\n\
\ |w<" ++ epsilonString ++ "=(f,0)\n\
\ |d==0=\n\
\  let r1=fromIntegral r\n\
\      b=(e0(-1 *a*w)i 0)\n\
\      l=((b+(a*w)-1)/(a*w*r1))\n\
\      m=1-((r1-1)*l)\n\
\  in (f*m,f*l)\n\
\ |otherwise=\n\
\  let r1=fromIntegral r\n\
\      l=(w/(r1*(w+a)))\n\
\      m=1-((r1-1)*l)\n\
\  in (f*m,f*l)\n"
{-
\j d r a i (w,f)=\n\
\ if w<" ++ epsilonString ++ " then (f,0)\n\
\ else if d==0 then\n\
\   let r1=fromIntegral r\n\
\     b=(e(-1*a*w)i 0)\n\
\     l=((b+(a*w)-1)/(a*w*r1))\n\
\     m=1-((r1-1)*l)\n\
\   in (f*m,f*l)\n\
\ else\n\
\   let r1=fromIntegral r\n\
\     l=(w/(r1*(w+a)))\n\
\     m=1-((r1-1)*l)\n\
\   in (f*m,f*l)\n"
-}

-- \x :: Double -> Int -> Double -> Int -> Int -> [Double]\n\
discreteGammaString :: String
discreteGammaString="\
\x0 a n m i e=\n\
\ if n==1 then [1]\n\
\ else\n\
\  let b=n*e\n\
\      v=m/(fromIntegral b)\n\
\      z=(k2(r0 b v 0)0)\n\
\      c=k2(gM(t0 a i v)z)0\n\
\      j=gM(+(1-(a3 c)))c\n\
\      s=w0 n j z 1 True\n\
\  in gM(*((fromIntegral n)/(c3(+)0 s)))s\n"

-- \k2:: [Double] -> Double -> [Double]\n\
cumulativeSumString :: String
cumulativeSumString="\
\k2(a:c)b=\n\
\ if nU c then [b+a]\n\
\ else (b+a):k2 c(b+a)\n"

-- \t :: Double -> Int -> Double -> Double -> Double\n\
gammaPDFString :: String
gammaPDFString="\
\t0 a i l r=l*(q0 a a i)*(q0 r(a-1)i)*(e0(-1*a*r)i 0)/(v0 a i 0 1)\n"

-- \w0 :: Int -> [Double] -> [Double] -> Int -> Bool-> [Double]\n\
getNTilesString :: String
getNTilesString="\
\w0 c d g t o=\n\
\ if nU d then []\n\
\ else\n\
\   let f=a5 d\n\
\       e=a5 g\n\
\       v=(fromIntegral t)/(2*(fromIntegral c))\n\
\   in\n\
\   if f>=v then\n\
\    if o then e:w0 c(a6 d)(a6 g)(t+1)False\n\
\    else w0 c(a6 d)(a6 g)(t+1)True\n\
\   else w0 c(a6 d)(a6 g)t o\n"

-- \q :: Double -> Double -> Int -> Double\n\
-- maybe replace with expE when used
expX2YString :: String
expX2YString="\
\q0 a x i=e0(x*(l0 a i 0 0))i 0\n"

-- \v0 :: Double -> Int -> Int -> Double -> Double\n\
gammaFunString :: String
gammaFunString="\
\v0 u i c l=\n\
\ if c>i then (fM$fromIntegral i)*(q0(fromIntegral i)u i)/l\n\
\ else v0 u i(c+1)(l*(u+(fromIntegral c)))\n"

-- \a9 :: Int -> Double -> Int -> (Double, Double) -> (Double, Double)\n\
neymanUniformWithKString :: String
neymanUniformWithKString="\
\a9 r c i(w,f)=\n\
\ if w<" ++ epsilonString ++ " then (f,0)\n\
\ else\n\
\  let r1=fromIntegral r\n\
\      b=(e0(-1*c*w)i 0)\n\
\      l=((b+(c*w)-1)/(c*w*r1))\n\
\      m=1-((r1-1)*l)\n\
\   in (f*m,f*l)\n"

-- \z :: Int -> Double -> Int -> (Double, Double) -> (Double, Double)\n\
neymanExponentialWithKString :: String
neymanExponentialWithKString="\
\z0 r a i(w,f)=\n\
\ if w<" ++ epsilonString ++ " then (f,0)\n\
\ else\n\
\ let r1=fromIntegral r\n\
\     l=(w/(r1*(w+a)))\n\
\     m=1-((r1-1)*l)\n\
\ in (f*m,f*l)\n"

-- \a1 :: Int -> Double -> Int -> Int -> Double -> Int -> [(Double, Double)]\n\
getModifierListSmallString :: String
getModifierListSmallString="\
\a1 v t h n a i\n\
\ |v==0 && h==0=[(1,1)]\n\
\ |v==1 && h==0=[(0,t),(1/(1-t),1-t)]\n\
\ |otherwise=\n\
\  let w=x0 a n 10.0 i(n*i)\n\
\      f=r0 n(1.0/(fromIntegral n))0\n\
\  in\n\
\  if t==0 then a2 w f else\n\
\  a2(0:(gM(*(1/(1-t)))w))(t:(gM(*(1-t))f))\n"
{-
\a1 v t h n a i=\n\
\ if v==0 && h==0 then [(1,1)]\n\
\ else if v==1 && h==0 then [(0,t),(1/(1-t),(1-t))]\n\
\ else\n\
\  let w=x a n " ++ maxGammaRateString ++ " i(n*i)\n\
\    f=r n (1.0/(fromIntegral n))0\n\
\  in\n\
\  if t==0 then a2 w f\n\
\  else a2(0:(g(*(1/(1-t)))w))(t:(g(*(1-t))f))\n"
-}

-- \a2 :: [a] -> [b] -> [(a,b)]\n\
zipString :: String
zipString="\
\a2(l:m)(n:o)=\n\
\ if nU m then [(l,n)]\n\
\ else (l,n):a2 m o\n"

-- \a3 :: [a] -> a\n\
lastString :: String
lastString="\
\a3(a:b)=\n\
\ if nU b then a\n\
\ else a3 b\n"

-- \a4 :: [[a]] -> [a]\n\
concatString :: String
concatString="\
\a4(a:b)=\n\
\ if nU b then a\n\
\ else a++a4 b\n"

headString :: String
headString="\
\a5(a:_)=a\n"
-- \a5[]=error"BadHead"\n"


tailString :: String
tailString="\
\a6(_:b)=b\n"

-- \a7 :: Int -> Double -> Int -> [(Double, Double)] -> String -> [[Double]]\n\
makeNeymanUniformMatrixString :: String
makeNeymanUniformMatrixString="\
\a7 l j p o s=\n\
\ let k=gM(a9 l j p)o\n\
\     q=c3(+)0$gM hM k\n\
\     c=c3(+)0$gM iM k\n\
\     e=(-1)*b0 q p 0 0\n\
\     f=(-1)*b0 c p 0 0\n\
\ in m0 l e f 0 s\n"

-- \a8 :: Int -> Double -> Int -> [(Double, Double)] -> String -> [[Double]]\n\
makeNeymanExponentialMatrixString :: String
makeNeymanExponentialMatrixString="\
\a8 l j p o s=\n\
\ let k=gM(z0 l j p)o\n\
\     q=c3(+)0$gM hM k\n\
\     c=c3(+)0$gM iM k\n\
\     e=(-1)*b0 q p 0 0\n\
\     f=(-1)*b0 c p 0 0\n\
\ in m0 l e f 0 s\n"

--- \b7 :: [[Double]] -> [[Double]]\n\
transposeMatrixString :: String
transposeMatrixString="\
\b7 a=\n\
\ if nU$a5 a then []\n\
\ else gM a5 a:b7(gM a6 a)\n"

absString :: String
absString="\
\c8 a=\n\
\ if a<0 then -1*a\n\
\ else a\n"

-- start at 0
-- \d5 :: Int -> [a] -> Int\n\
lengthString :: String
lengthString="\
\d5 b a=\n\
\ if nU a then b\n\
\ else d5(b+1)(a6 a)\n"

-- \c0 :: [[Double]] -> [[Double]]\n\
invertMatrixString :: String
invertMatrixString="\
\c0 a=d8(1/(d7 a))(d6 a 1)\n"

-- \zipWith :: (Double -> Double -> Double) -> [Double] -> [Double] -> [Double]\n\
zipWithString :: String
zipWithString="\
\e5 f a b=\n\
\ if nU a ||nU b then []\n\
\ else f(a5 a)(a5 b):e5 f(a6 a)(a6 b)\n"

-- \b2 :: Int -> [Double] -> [[Double]]\n\
split2MatrixString :: String
split2MatrixString="\
\b2 a b=\n\
\ if nU b then []\n\
\ else (d3 0 a b):b2 a(d4 0 a b)\n"

-- \c4 :: [[Double]] -> [[Double]] -> [[Double]]\n\
addMatricesString :: String
addMatricesString="\
\c4 a b=\n\
\ if nU a && nU b then []\n\
\ else (e5(+)(a5 a)(a5 b)):c4(a6 a)(a6 b)\n"

-- \c3 :: (a -> b -> a) -> a -> [b] -> a\n\
foldlString :: String
foldlString="\
\c3 f a b=\n\
\ if nU b then a\n\
\ else c3 f(f a$a5 b)(a6 b)\n"

-- \c2 :: [Double] -> [[Double]] -> [[Double]]-> Int -> Int -> Double -> Double -> Int -> Int -> (Double -> Int -> Double -> Double) -> (Double, Double) -> [Double]\n\
--2.0 is maxTime--should be String
integrateGTRMatrixWithKString :: String
integrateGTRMatrixWithKString="\
\c2 a b c d e f g h i j(k,l)\n\
\ |d==i=[]\n\
\ |e==i=c2 a b c(d+1)0 f g h i j(k,l)\n\
\ |k<" ++ epsilonString ++ "=0:(c2 a b c d(e+1)f g h i j(k,l))\n\
\ |otherwise=(l*(e6 e7 j a b c d e f 2.0 h 0 k)):(c2 a b c d(e+1)f g h i j(k,l))\n"


-- \c5 :: [[Double]] -> Int -> Int -> [Double]\n\
adjustSymString :: String
adjustSymString="\
\c5 a b c=\n\
\ let e=c5 a b(c+1)\n\
\ in\n\
\ if b==d5 0 a then []\n\
\ else if c==d5 0 a then c5 a(b+1)0\n\
\ else if c==b then ((a!!b)!!c):e\n\
\ else\n\
\ let d=((a!!b)!!c+(a!!c)!!b)/2\n\
\ in\n\
\ if d>0 then d:e\n\
\ else " ++ epsilonString ++ ":e\n"

-- \c6 :: [[Double]] -> [[Double]] -> Int -> [[Double]]\n\
adjustDiagString :: String
adjustDiagString="\
\c6 a b c=\n\
\ if c==d5 0 b then []\n\
\ else\n\
\ let e=a5 a\n\
\ in ((d3 0 c e)++[1-((c3(+)0 e)-((b!!c)!!c))]++(d4 0(c+1)e)):c6(a6 a)b(c+1)\n"

-- \e7 :: [Double] -> [[Double]] -> [[Double]] -> Double -> Int -> Int -> Int -> Int -> Double\n\
getPijString :: String
getPijString="\
\e7 a b c d f g h i=\n\
\ if f==d5 0 b then 0\n\
\ else (e0((a!!f)*d)i 0*((b!!g)!!f)*((c!!f)!!h))+e7 a b c d(f+1)g h i\n"

-- removeColumn :: Int -> [[a]] -> [[a]]\n\
removeColumnString :: String
removeColumnString="\
\e9 a b=\n\
\ if nU b then []\n\
\ else\n\
\ let c=a5 b\n\
\ in ((d3 0(a-1)c)++(d4 0 a c)):(e9 a$a6 b)\n"

-- \f0 :: Int -> Int -> [[a]] -> [[a]]\n\
removeRowAndColumnString :: String
removeRowAndColumnString="\
\f0 a b c=\n\
\ if nU c then []\n\
\ else e9 b ((d3 0 (a -1) c) ++ (d4 0 a c))\n"

-- \d8 :: Double -> [[Double]] -> [[Double]]\n\
matrixMultiplyScalarString :: String
matrixMultiplyScalarString="\
\d8 a b=\n\
\ if nU b then []\n\
\ else gM(a *)(a5 b):d8 a(a6 b)\n"

-- \d7 :: [[Double]] -> Double\n\
determinantNumericalString :: String
determinantNumericalString="\
\d7 e=\n\
\ if (d5 0$a5 e)==2 then\n\
\ let [a,b]=a5 e\n\
\     [c,d]=a3 e\n\
\ in a*d-b*c\n\
\ else f2 1 e\n"

-- \f2 :: Int -> [[Double]] -> Double\n\
getCofactor1String :: String
getCofactor1String="\
\f2 a b=\n\
\ if (a>(d5 0$a5 b)) then 0\n\
\ else (((a5 b)!!(a-1))*(p1(-1)(1+a)0)*(d7(e9 a(a6 b))))+f2(a+1)b\n"

-- \d6 :: [[Double]] -> Int -> [[Double]]\n\
cofactorTMatrixString :: String
cofactorTMatrixString="\
\d6 a b=\n\
\ if b>(d5 0$a5 a) then []\n\
\ else f3 b 1 a:d6 a(b+1)\n"

-- \e8 :: Double -> Int -> Double -> Double\n\
getUniformPdfString :: String
getUniformPdfString="\
\e8 a b c=1/a\n"

-- numerical issue?
-- a*exp (-1*a*c)\n"
-- gM[a*(e(-1 *a*c)b 0)]0\n"
-- \h1 a b c=a*(e(-1*a*c)b 0)\n"
-- \getExponentialPdf :: Double -> Int -> Double -> Double\n\
getExponentialPdfString :: String
getExponentialPdfString="\
\h1 a b c=a*(e0(-1*a*c)b 0)\n"

-- \f7 :: (a -> Bool) -> [a] -> [a]\n\
takeWhileString :: String
takeWhileString="\
\f7 f a\n\
\ |nU a=[]\n\
\ |(not$f(a5 a))=[]\n\
\ |otherwise=(a5 a):(f7 f(a6 a))\n"
{-
\f7 f a=\n\
\ if nU a then []\n\
\ else if (not$f(a5 a)) then []\n\
\ else (a5 a):(f7 f(a6 a))\n"
-}

-- \b1 :: [[Double]] -> [[Double]]\n\
regularizeRString :: String
regularizeRString="\
\b1 a=\n\
\ let b=fromIntegral$d5 0 a\n\
\ in gM(gM(*((((b*b)-b)/2)/(c3(+)0(gM(c3(+)0)(gM(f7(>0))a))))))a\n"

-- need to start with 0 as first arg
-- \d3 :: Int -> Int -> [a] -> [a]\n\
takeString :: String
takeString="\
\d3 m n a\n\
\ |nU a=[]\n\
\ |m==n=[]\n\
\ |otherwise=a5 a:d3(m+1)n(a6 a)\n"
{-
\d3 m n a=\n\
\ if nU a then []\n\
\ else if m==n then []\n\
\ else a5 a:d3(m+1)n(a6 a)\n"
-}

-- need to start with first arg 0
-- \d4 :: Int -> Int -> [a] -> [a]\n\
dropString :: String
dropString="\
\d4 m n a\n\
\ |nU a=[]\n\
\ |m==n=a\n\
\ |otherwise=d4(m+1)n(a6 a)\n"
{-
\d4 m n a=\n\
\ if nU a then []\n\
\ else if m==n then a\n\
\ else d4(m+1)n(a6 a)\n"
-}

{-
-- \reverse :: [a] -> [a]\n\
reverseString :: String
reverseString="\
\f8 a=\n\
\ if nU a then []\n\
\ else\n\
\ a3 a:f8(f6 a)\n"
-}

-- \b4 :: [[Double]] -> Int -> Int -> Int -> [Double]\n\
addDiagValuesString :: String
addDiagValuesString="\
\b4 a b c e\n\
\ |c==b=[]\n\
\ |e==b=b4 a b(c+1)0\n\
\ |c/=e=((a!!c)!!e):b4 a b c(e+1)\n\
\ |otherwise=(-1 *(c3(+)0(a!!c))):b4 a b c(e+1)\n"
{-
\b4 a b c e=\n\
\ if c==b then []\n\
\ else if e==b then b4 a b(c+1)0\n\
\ else if c/=e then ((a!!c)!!e):b4 a b c(e+1)\n\
\ else (-1 *(c3(+)0(a!!c))):b4 a b c(e+1)\n"
-}

-- \b3 :: [[Double]] -> [Double] -> Int -> Int -> Int -> [Double]\n\
makeQString :: String
makeQString="\
\b3 a b c d e\n\
\ |d==c=[]\n\
\ |e==c=b3 a b c(d+1)0\n\
\ |d==e=0:b3 a b c d(e+1)\n\
\ |otherwise=((a!!d)!!e)*(b!!e):b3 a b c d(e+1)\n"
{-
\b3 a b c d e=\n\
\ if d==c then []\n\
\ else if e==c then b3 a b c(d+1)0\n\
\ else if d==e then 0:b3 a b c d(e+1)\n\
\ else ((a!!d)!!e)*(b!!e):b3 a b c d(e+1)\n"
-}

-- \f3 :: Int -> Int -> [[Double]] -> [Double]\n\
getRowString :: String
getRowString="\
\f3 a b c=\n\
\ if b>(d5 0$a5 c) then []\n\
\ else (p1(-1)(b+a)0)*(d7$f0 b a c):f3 a(b+1)c\n"

-- /f9 :: Double -> Int -> Int -> Double
-- value iterations counter (start at 0)
sqrtString :: String
sqrtString="\
\f9 x i c=e0(0.5*l0 x i c 0)i c\n"


-- \g0 ::(Ord a)=> [a] -> a -> a\n\
maximumString :: String
maximumString="\
\g0 a c\n\
\ |nU a=c\n\
\ |c>a5 a=g0(a6 a)c\n\
\ |otherwise=g0(a6 a)(a5 a)\n"

{-
\g0 a c=\n\
\ if nU a then c\n\
\ else if c>a5 a then g0(a6 a)c\n\
\ else g0(a6 a)(a5 a)\n"
-}


-- \g1 ::(Ord a)=> [a] -> a -> a\n\
minimumString :: String
minimumString="\
\g1 a c\n\
\ |nU a=c\n\
\ |c<a5 a=g1(a6 a)c\n\
\ |otherwise=g1(a6 a)(a5 a)\n"


{-
\ if nU a then c\n\
\ else if c<(a5 a) then g1(a6 a)c\n\
\ else g1(a6 a)(a5 a)\n"
-}

-- \b8 :: Int -> [Double]-> [Double]\n\
normalizeVectorWithSignString :: String
normalizeVectorWithSignString="\
\b8 i a=\n\
\ let s=f9(c3(+)0$e5(*)a a)i 0\n\
\ in\n\
\ if g0(a6 a)(a5 a)>=(c8$g1(a6 a)(a5 a)) then gM(/s)a\n\
\ else gM(/(-1 *s))a\n"

-- \g2 :: (a -> Bool) -> [a] -> [a]\n\
filterString :: String
filterString="\
\g2 f b\n\
\ |nU b=[]\n\
\ |f(a5 b)=(a5 b):g2 f(a6 b)\n\
\ |otherwise=g2 f(a6 b)\n"


-- \g3 :: [[Double]] -> Int -> [[Double]]\n\
reduceRowString :: String
reduceRowString="\
\g3 c r=\n\
\ let d=(g4 c r$a5$g2(\\x -> c!!x!!r /=0)[r..(d5 0 c)-1])\n\
\     e=d!!r\n\
\     f=g(\\x -> x/(e!!r))e\n\
\     h nr=let k=nr!!r in e5(\\a b -> k*a -b)f nr\n\
\     i=g h$d4 0(r+1)d\n\
\ in d3 0 r d++[f]++i\n"

-- \fixlastrow :: [[Double]] -> [[Double]]\n\
fixlastrowString :: String
fixlastrowString="\
\g5 b=(f6 b)++[f6(f6(a3 b))++[1,(a3(a3 b))/a3(f6(a3 b))]]\n"


-- \substitute :: [[Double]] -> [Double]\n\
substituteString :: String
substituteString="\
\g6 a=h0 b [a3(a3 a)] (f6 a) where\n\
\ b c e=a3 c -c3(+)0(e5(*)e$f6$d4 0(d5 0 a -d5 0 e)c):e\n"


-- \g9 :: Int -> [a] -> ([a],[a])\n\
splitAtString :: String
splitAtString="\
\g9 a b=(d3 0 a b,d4 0 a b)\n"


-- \h0 :: (a -> b -> b) -> b -> [a] -> b\n\
foldrString :: String
foldrString="\
\h0 f b a=\n\
\ if nU a then b\n\
\ else h0 f (f(a5 a)b) (a6 a)\n"


--error \"\"\n\ if maxiterations exceeded precision issue?
-- \h3 :: [[Double]] -> [[Double]] -> Int -> Int -> ([[Double]], [[Double]], [[Double]])\n\
qrFactorizationString :: String
qrFactorizationString="\
\h3 a u c i=\n\
\ let (q,r)=h4 a i\n\
\     n=h2 r q\n\
\     o=h2 u q\n\
\     l=a4 a\n\
\     m=a4 n\n\
\ in\n\
\ if (((c3(+)0$gM c8$e5(-)l m)/fromIntegral (d5 0 a * d5 0 a))<" ++ epsilonString ++ ") ||(c<" ++ maxIterationsString ++ ") then (q,r,o)\n\
\ else h3 n o (c+1)i\n"

-- \h5 :: [[Double]] -> Int -> [Double]\n\
getDiagValuesString :: String
getDiagValuesString="\
\h5 i c\n\
\ |nU i=error \"\"\n\
\ |c==d5 0 i=[]\n\
\ |otherwise=((i!!c)!!c):h5 i(c+1)\n"
{-
\h5 i c=\n\
\ if nU i then error \"\"\n\
\ else if c==d5 0 i then []\n\
\ else ((i!!c)!!c):h5 i (c+1)\n"
-}

-- \h4 :: [[Double]] -> Int -> ([[Double]], [[Double]])\n\
qrDecompositionString :: String
qrDecompositionString="\
\h4 a i=\n\
\ let c=gM b7(h6 a a(d5 0 a)0 i)\n\
\     d=c3 h2(a5 c)(a6 c)\n\
\ in (d,h2(b7 d)a)\n"

-- \h2 :: [[Double]] -> [[Double]] -> [[Double]]\n\
matrixMultiplyString :: String
matrixMultiplyString="\
\h2 a b=h7 0(d5 0 a)(d5 0$a5 b)a b\n"

-- \h6 :: [[Double]] -> [[Double]] -> Int -> Int -> Int -> [[[Double]]]\n\
getHouseholderListString :: String
getHouseholderListString="\
\h6 a b d c i=\n\
\ let j=h8 a i\n\
\     k=h9 j d\n\
\ in\n\
\ if d5 0 j==2 then [k]\n\
\ else k:h6(i0 c c(h2 k b))b d(c+1)i\n"

-- \h8 :: [[Double]]-> Int -> [[Double]]\n\
getHouseholderString :: String
getHouseholderString="\
\h8 a i=\n\
\ let v=i2 0 a\n\
\     e=i5(i4 v(d8(i3(a4 v)i)(b7[i1(d5 0 a)0 1])))i\n\
\ in i4(i6(d5 0 a)1)(d8 2$h2 e(b7 e))\n"

-- \h9 :: [[Double]] ->Int -> [[Double]]\n\
padOutMinorString :: String
padOutMinorString="\
\h9 a d=\n\
\ if nU a then i6 d 1 --identity matrix\n\
\ else\n\
\  let m=d5 0 a\n\
\  in\n\
\  if m==d then a\n\
\  else\n\
\   let r=d-m\n\
\       l=i6 d 1\n\
\   in d3 0 r l++e5(++)(gM(d3 0 r)$d4 0 r l)a\n"

-- \i0 :: Int -> Int -> [[Double]] -> [[Double]]\n\
makeMatrixMinorString :: String
makeMatrixMinorString="\
\i0 r c a=gM(d4 0(c+1))(d4 0(r+1)a)\n"

-- \i6 :: Int -> Double -> [[Double]]\n\
makeDiagMatrixString :: String
makeDiagMatrixString="\
\i6 d v=i7 d v 0\n"

-- \i2 :: Int -> [[Double]] -> [[Double]]\n\
getColumnVectorString :: String
getColumnVectorString="\
\i2 c a=gM(:[])$gM(!!c)a\n"

-- \i1 :: Int -> Int -> Double -> [Double]\n\
makeEVectorString :: String
makeEVectorString="\
\i1 d p v=r0 p 0 0++[v]++r0(d-p-1)0 0\n"

-- \i4 :: [[Double]] -> [[Double]] -> [[Double]]\n\
subtractMatricesString :: String
subtractMatricesString="\
\i4 a b=\n\
\ if nU a && nU b then []\n\
\ else e5(-)(a5 a)(a5 b):i4(a6 a)(a6 b)\n"

-- \i5 :: [[Double]] -> Int -> [[Double]]\n\
normalizeColumnVectorString :: String
normalizeColumnVectorString="\
\i5 c i=gM(:[])(i8(a4 c)i)\n"

-- \i7 :: Int -> Double -> Int -> [[Double]]\n\
makeDiagRowString :: String
makeDiagRowString="\
\i7 d v c=\n\
\ if c==d then []\n\
\ else (r0 c 0 0++[v]++r0(d-c-1)0 0):i7 d v(c+1)\n"

-- \h7 :: Int -> Int -> Int -> [[Double]] -> [[Double]] -> [[Double]]\n\
getRowsString :: String
getRowsString="\
\h7 c a b d e=\n\
\ if c==a then []\n\
\ else i9 c 0 b d e:h7(c+1)a b d e\n"

-- \i9 :: Int -> Int -> Int -> [[Double]] -> [[Double]] -> [Double]\n\
getElementString :: String
getElementString="\
\i9 a c b d e=\n\
\ if c==b then []\n\
\ else c3(+)0(e5(*)(d!!a)(gM(!!c)e)):i9 a(c+1)b d e\n"

-- \i8 :: [Double] -> Int -> [Double]\n\
normalizeVectorString :: String
normalizeVectorString="\
\i8 a i=gM(/i3 a i)a\n"

-- \i3 :: [Double] -> Int -> Double\n\
euclidNormString :: String
euclidNormString="\
\i3 a i=f9(c3(+)0(e5(*)a a))i 0\n"

-- \j0 :: [[Double]] -> Int -> Int -> Int -> [Double]\n\
isPosRString :: String
isPosRString="\
\j0 i a r c\n\
\ |r==a=[]\n\
\ |c==a=j0 i a(r+1)0\n\
\ |r==c=0:j0 i a r(c+1)\n\
\ |otherwise=((i!!r)!!c):j0 i a r(c+1)\n"
{-
\j0 i a r c=\n\
\ if r==a then []\n\
\ else if c==a then j0 i a(r+1)0\n\
\ else if r==c then 0:j0 i a r(c+1)\n\
\ else ((i!!r)!!c):j0 i a r(c+1)\n"
-}

-- \makeGTRMatrixLocal :: Int -> [[Double]] -> [Double]-> Int -> ([Double], [[Double]], [[Double]])\n\
-- \a0 :: Int -> [[Double]] -> [Double]-> Int -> ([Double], [[Double]], [[Double]])\n\
makeGTRMatrixLocalString :: String
makeGTRMatrixLocalString="\
\a0 a r p i=\n\
\ let q=(b2 a$b4(b2 a$b3(b1(b2 a$j0 r a 0 0))p a 0 0)a 0 0)\n\
\     (j,k,l)=h3 q(i6(d5 0 q)1)0 i\n\
\ in (h5(h2 k j)0,l,c0 l)\n"

-- \c1 :: (Double -> Int-> Int -> Double -> Double) -> String -> [Double] -> [[Double]] -> [[Double]] -> Int -> Double -> Double -> Int -> (Double -> Int -> Double -> Double) -> [(Double, Double)] -> [[Double]]\n\
makeGTRLogMatrixString :: String
makeGTRLogMatrixString="\
\c1 l m e u v a p x i d f=\n\
\ let y=b2 a$c5(c3 c4(r0 a(r0 a 0.0 0)0)(gM(b2 a)$gM(c2 e u v 0 0 p x i a d)f))0 0\n\
\ in b2 a$gM(*(-1))$c7 l(c6 y y 0)a m 0 0 i\n"

-- \makeGTRLogMatrix4State :: (Double -> Int-> Int -> Double -> Double) -> ([Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]) -> Double -> Int -> [Double] -> [Double] -> [(Double, Double)] -> [[Double]]\n
makeGTRLogMatrix4StateString :: String
makeGTRLogMatrix4StateString="\
\j1 l m p i o q f=\n\
\ let y=b2 4$c5 (c3 c4(r0 4(r0 4 0.0 0)0)(gM(m o q p i)f)) 0 0\n\
\ in b2 4$gM(* (-1))$c7 l(c6 y y 0)4 \"T\" 0 0 i\n"

-- \c7 :: (Double -> Int-> Int -> Double -> Double) -> [[Double]] -> Int -> String -> Int -> Int -> Int -> [Double]\n\
getLogMatrixString :: String
getLogMatrixString="\
\c7 l m a e r c i\n\
\ |r==a=[]\n\
\ |c==a=c7 l m a e(r+1)0 i\n\
\ |(r==(a-1))&&(c==(a-1))=if e==\"-\" then [0] else (l((m!!r)!!c)i 0 0):c7 l m a e r(c+1)i\n\
\ |otherwise=(l((m!!r)!!c)i 0 0):c7 l m a e r(c+1)i\n"
{-
\c7 l m a e r c i=\n\
\ if r==a then []\n\
\ else if c==a then c7 l m a e(r+1)0 i\n\
\ else if (r==(a-1))&&(c==(a-1)) then\n\
\ if e==\"-\" then [0]\n\
\ else (l((m!!r)!!c)i 0 0):c7 l m a e r(c+1)i\n\
\ else (l((m!!r)!!c)i 0 0):c7 l m a e r(c+1)i\n"
-}

-- \e6 :: ([Double] -> [[Double]] -> [[Double]] -> Double -> Int -> Int -> Int -> Int -> Double) -> (Double -> Int -> Double -> Double) -> [Double] -> [[Double]] -> [[Double]]-> Int -> Int -> Double -> Double -> Int -> Int -> Double -> Double\n\
trapezoidIntegrationString :: String
trapezoidIntegrationString="\
\e6 f s e u v r j p x i c k=\n\
\ if c==2*i then 0\n\
\ else\n\
\ let l=x/fromIntegral(2*i)\n\
\     t=fromIntegral c*l\n\
\ in (l*((f e u v(t*k)0 r j i*s p i t)+(f e u v((t*k)+l)0 r j i*s p i (t+l)))/2)+e6 f s e u v r j p x i(c+1)k\n"


{-
The variable names with letter number and the 'l' like k2l are to avoid name shadowing with functions.
some simpler models have unused parasm to allow for function as parameter
-}

-- \tn93ExponentialWithK :: [Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]\n\
tn93ExponentialWithKString :: String
tn93ExponentialWithKString="\
\j2 [a,c,b] [pA,pC,pG,pT] x i (w,k)=\n\
\ if w<" ++ epsilonString ++ " then [[k,0,0,0],[0,k,0,0],[0,0,k,0],[0,0,0,k]]\n\
\ else\n\
\ let pR=pA+pG\n\
\     pY=pC+pT\n\
\     f1l=(pY*a+pR*b)\n\
\     f2l=(pR*c+pY*b)\n\
\     c0l=pA*pY/pR\n\
\     c1l=pC*pR/pY\n\
\     c2l=pG*pY/pR\n\
\     c3l=pT*pR/pY\n\
\     d0l=pG/pR\n\
\     d1l=pT/pY\n\
\     d2l=pA/pR\n\
\     d3l=pC/pY\n\
\     m2l=x*x\n\
\     km=x*w\n\
\     k2l=w*w\n\
\     p00=(((d0l+c0l+pA)*m2l)+((((c0l+pA)*f1l)+(b*d0l)+(pA *b))*km)+(pA*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p01=pC*b*w/(x+(b*w))\n\
\     p02=(-1)*(((d0l-c2l-pG)*m2l)+(((((-1*c2l)-pG)*f1l)+(b*d0l)-(pG *b))*km)-(pG*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p03=pT*b*w/(x+(b*w))\n\
\     p10=pA*b*w/(x+(b*w))\n\
\     p11=(((d1l+c1l+pC)*m2l)+((((c1l+pC)*f2l)+(b*d1l)+(pC *b))*km)+(pC*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\     p12=pG*b*w/(x+(b*w))\n\
\     p13=(-1)*(((d1l-c3l-pT)*m2l)+(((((-1*c3l)-pT)*f2l)+(b*d1l)-(pT *b))*km)-(pT*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\     p20=(-1)*(((d2l-c0l-pA)*m2l)+(((((-1*c0l)-pA)*f1l)+(b*d2l)-(pA *b))*km)-(pA*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p21=pC*b*w/(x+(b*w))\n\
\     p22=(((d2l+c2l+pG)*m2l)+((((c2l+pG)*f1l)+(b*d2l)+(pG *b))*km)+(pG*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p23=pT*b*w/(x+(b*w))\n\
\     p30=pA*b*w/(x+(b*w))\n\
\     p31=(-1)*(((d3l-c1l-pC)*m2l)+(((((-1*c1l)-pC)*f2l)+(b*d3l)-(pC *b))*km)-(pC*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\     p32=pG*b*w/(x+(b*w))\n\
\     p33=(((d3l+c3l+pT)*m2l)+((((c3l+pT)*f2l)+(b*d3l)+(pT *b))*km)+(pT*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\ in d8 k [[p00,p01,p02,p03],[p10,p11,p12,p13],[p20,p21,p22,p23],[p30,p31,p32,p33]]\n"

-- \tn93UniformWithK :: [Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]\n\
tn93UniformWithKString :: String
tn93UniformWithKString="\
\j3 [a,c,b] [pA,pC,pG,pT] x i (w,k)=\n\
\ if w<" ++ epsilonString ++ " then [[k,0,0,0],[0,k,0,0],[0,0,k,0],[0,0,0,k]]\n\
\ else\n\
\ let pR=pA+pG\n\
\     pY=pC+pT\n\
\     f1l=(pY*a+pR*b)\n\
\     f2l=(pR*c+pY*b)\n\
\     c0l=pA*pY/pR\n\
\     c1l=pC*pR/pY\n\
\     c2l=pG*pY/pR\n\
\     c3l=pT*pR/pY\n\
\     d0l=pG/pR\n\
\     d1l=pT/pY\n\
\     d2l=pA/pR\n\
\     d3l=pC/pY\n\
\     ebkm=e0(b*w*x)i 0\n\
\     enbkm=e0(-1*b*w*x)i 0\n\
\     enf1lkm=e0(-1*f1l*w*x)i 0\n\
\     enf2lkm=e0(-1*f2l*w*x)i 0\n\
\     bkmfactor=(enbkm+((b*w*x)-1))/(b*w*x)\n\
\     p00=((enbkm*((ebkm*((pA*b*f1l*w*x)+(c0l*f1l)+(b*d0l)))-(c0l*f1l)))-(b*d0l*enf1lkm))/(b*f1l*w*x)\n\
\     p01=pC*bkmfactor\n\
\     p02=((enbkm*((ebkm*((pG*b*f1l*w*x)+(c2l*f1l)-(b*d0l)))-(c2l*f1l)))+(b*d0l*enf1lkm))/(b*f1l*w*x)\n\
\     p03=pT*bkmfactor\n\
\     p10=pA*bkmfactor\n\
\     p11=((enbkm*((ebkm*((pC*b*f2l*w*x)+(c1l*f2l)+(b*d1l)))-(c1l*f2l)))-(b*d1l*enf2lkm))/(b*f2l*w*x)\n\
\     p12=pG*bkmfactor\n\
\     p13=((enbkm*((ebkm*((pT*b*f2l*w*x)+(c3l*f2l)-(b*d1l)))-(c3l*f2l)))+(b*d1l*enf2lkm))/(b*f2l*w*x)\n\
\     p20=((enbkm*((ebkm*((pA*b*f1l*w*x)+(c0l*f1l)-(b*d2l)))-(c0l*f1l)))+(b*d2l*enf1lkm))/(b*f1l*w*x)\n\
\     p21=pC*bkmfactor\n\
\     p22=((enbkm*((ebkm*((pG*b*f1l*w*x)+(c2l*f1l)+(b*d2l)))-(c2l*f1l)))-(b*d2l*enf1lkm))/(b*f1l*w*x)\n\
\     p23=pG*bkmfactor\n\
\     p30=pA*bkmfactor\n\
\     p31=((enbkm*((ebkm*((pC*b*f2l*w*x)+(c1l*f2l)-(b*d3l)))-(c1l*f2l)))+(b*d3l*enf2lkm))/(b*f2l*w*x)\n\
\     p32=pG*bkmfactor\n\
\     p33=((enbkm*((ebkm*((pT*b*f2l*w*x)+(c3l*f2l)+(b*d3l)))-(c3l*f2l)))-(b*d3l*enf2lkm))/(b*f2l*w*x)\n\
\ in d8 k [[p00,p01,p02,p03],[p10,p11,p12,p13],[p20,p21,p22,p23],[p30,p31,p32,p33]]\n"

-- \f84ExponentialWithK :: [Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]\n\
f84ExponentialWithKString :: String
f84ExponentialWithKString="\
\j4 [p, b] [pA,pC,pG,pT] x i (w,k)=\n\
\ if w<" ++ epsilonString ++ " then [[k,0,0,0],[0,k,0,0],[0,0,k,0],[0,0,0,k]]\n\
\ else\n\
\ let pR=pA+pG\n\
\     pY=pC+pT\n\
\     a=(1+(p/pY))*b\n\
\     c=(1+(p/pR))*b\n\
\     f1l=(pY*a+pR*b)\n\
\     f2l=(pR*c+pY*b)\n\
\     c0l=pA*pY/pR\n\
\     c1l=pC*pR/pY\n\
\     c2l=pG*pY/pR\n\
\     c3l=pT*pR/pY\n\
\     d0l=pG/pR\n\
\     d1l=pT/pY\n\
\     d2l=pA/pR\n\
\     d3l=pC/pY\n\
\     m2l=x*x\n\
\     km=x*w\n\
\     k2l=w*w\n\
\     p00=(((d0l+c0l+pA)*m2l)+((((c0l+pA)*f1l)+(b*d0l)+(pA *b))*km)+(pA*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p01=pC*b*w/(x+(b*w))\n\
\     p02=(-1)*(((d0l-c2l-pG)*m2l)+(((((-1*c2l)-pG)*f1l)+(b*d0l)-(pG *b))*km)-(pG*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p03=pT*b*w/(x+(b*w))\n\
\     p10=pA*b*w/(x+(b*w))\n\
\     p11=(((d1l+c1l+pC)*m2l)+((((c1l+pC)*f2l)+(b*d1l)+(pC *b))*km)+(pC*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\     p12=pG*b*w/(x+(b*w))\n\
\     p13=(-1)*(((d1l-c3l-pT)*m2l)+(((((-1*c3l)-pT)*f2l)+(b*d1l)-(pT *b))*km)-(pT*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\     p20=(-1)*(((d2l-c0l-pA)*m2l)+(((((-1*c0l)-pA)*f1l)+(b*d2l)-(pA *b))*km)-(pA*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p21=pC*b*w/(x+(b*w))\n\
\     p22=(((d2l+c2l+pG)*m2l)+((((c2l+pG)*f1l)+(b*d2l)+(pG *b))*km)+(pG*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p23=pT*b*w/(x+(b*w))\n\
\     p30=pA*b*w/(x+(b*w))\n\
\     p31=(-1)*(((d3l-c1l-pC)*m2l)+(((((-1*c1l)-pC)*f2l)+(b*d3l)-(pC *b))*km)-(pC*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\     p32=pG*b*w/(x+(b*w))\n\
\     p33=(((d3l+c3l+pT)*m2l)+((((c3l+pT)*f2l)+(b*d3l)+(pT *b))*km)+(pT*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\ in d8 k [[p00,p01,p02,p03],[p10,p11,p12,p13],[p20,p21,p22,p23],[p30,p31,p32,p33]]\n"

-- \f84UniformWithK :: [Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]\n\
f84UniformWithKString :: String
f84UniformWithKString="\
\j5 [p, b] [pA,pC,pG,pT] x i (w,k)=\n\
\ if w<" ++ epsilonString ++ " then [[k,0,0,0],[0,k,0,0],[0,0,k,0],[0,0,0,k]]\n\
\ else\n\
\ let pR=pA+pG\n\
\     pY=pC+pT\n\
\     a=(1+(p/pY))*b\n\
\     c=(1+(p/pR))*b\n\
\     f1l=(pY*a+pR*b)\n\
\     f2l=(pR*c+pY*b)\n\
\     c0l=pA*pY/pR\n\
\     c1l=pC*pR/pY\n\
\     c2l=pG*pY/pR\n\
\     c3l=pT*pR/pY\n\
\     d0l=pG/pR\n\
\     d1l=pT/pY\n\
\     d2l=pA/pR\n\
\     d3l=pC/pY\n\
\     ebkm=e0(b*w*x)i 0\n\
\     enbkm=e0(-1*b*w*x)i 0\n\
\     enf1lkm=e0(-1*f1l*w*x)i 0\n\
\     enf2lkm=e0(-1*f2l*w*x)i 0\n\
\     bkmfactor=(enbkm+((b*w*x)-1))/(b*w*x)\n\
\     p00=((enbkm*((ebkm*((pA*b*f1l*w*x)+(c0l*f1l)+(b*d0l)))-(c0l*f1l)))-(b*d0l*enf1lkm))/(b*f1l*w*x)\n\
\     p01=pC*bkmfactor\n\
\     p02=((enbkm*((ebkm*((pG*b*f1l*w*x)+(c2l*f1l)-(b*d0l)))-(c2l*f1l)))+(b*d0l*enf1lkm))/(b*f1l*w*x)\n\
\     p03=pT*bkmfactor\n\
\     p10=pA*bkmfactor\n\
\     p11=((enbkm*((ebkm*((pC*b*f2l*w*x)+(c1l*f2l)+(b*d1l)))-(c1l*f2l)))-(b*d1l*enf2lkm))/(b*f2l*w*x)\n\
\     p12=pG*bkmfactor\n\
\     p13=((enbkm*((ebkm*((pT*b*f2l*w*x)+(c3l*f2l)-(b*d1l)))-(c3l*f2l)))+(b*d1l*enf2lkm))/(b*f2l*w*x)\n\
\     p20=((enbkm*((ebkm*((pA*b*f1l*w*x)+(c0l*f1l)-(b*d2l)))-(c0l*f1l)))+(b*d2l*enf1lkm))/(b*f1l*w*x)\n\
\     p21=pC*bkmfactor\n\
\     p22=((enbkm*((ebkm*((pG*b*f1l*w*x)+(c2l*f1l)+(b*d2l)))-(c2l*f1l)))-(b*d2l*enf1lkm))/(b*f1l*w*x)\n\
\     p23=pG*bkmfactor\n\
\     p30=pA*bkmfactor\n\
\     p31=((enbkm*((ebkm*((pC*b*f2l*w*x)+(c1l*f2l)-(b*d3l)))-(c1l*f2l)))+(b*d3l*enf2lkm))/(b*f2l*w*x)\n\
\     p32=pG*bkmfactor\n\
\     p33=((enbkm*((ebkm*((pT*b*f2l*w*x)+(c3l*f2l)+(b*d3l)))-(c3l*f2l)))-(b*d3l*enf2lkm))/(b*f2l*w*x)\n\
\ in d8 k [[p00,p01,p02,p03],[p10,p11,p12,p13],[p20,p21,p22,p23],[p30,p31,p32,p33]]\n"

-- \hky85ExponentialWithK :: [Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]\n\
hky85ExponentialWithKString :: String
hky85ExponentialWithKString="\
\j6 [h, b] [pA,pC,pG,pT] x i (w,k)=\n\
\ if w<" ++ epsilonString ++ " then [[k,0,0,0],[0,k,0,0],[0,0,k,0],[0,0,0,k]]\n\
\ else\n\
\ let pR=pA+pG\n\
\     pY=pC+pT\n\
\     f1l=(pY*h+pR*b)\n\
\     f2l=(pR*h+pY*b)\n\
\     c0l=pA*pY/pR\n\
\     c1l=pC*pR/pY\n\
\     c2l=pG*pY/pR\n\
\     c3l=pT*pR/pY\n\
\     d0l=pG/pR\n\
\     d1l=pT/pY\n\
\     d2l=pA/pR\n\
\     d3l=pC/pY\n\
\     m2l=x*x\n\
\     km=x*w\n\
\     k2l=w*w\n\
\     p00=(((d0l+c0l+pA)*m2l)+((((c0l+pA)*f1l)+(b*d0l)+(pA *b))*km)+(pA*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p01=pC*b*w/(x+(b*w))\n\
\     p02=(-1)*(((d0l-c2l-pG)*m2l)+(((((-1*c2l)-pG)*f1l)+(b*d0l)-(pG *b))*km)-(pG*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p03=pT*b*w/(x+(b*w))\n\
\     p10=pA*b*w/(x+(b*w))\n\
\     p11=(((d1l+c1l+pC)*m2l)+((((c1l+pC)*f2l)+(b*d1l)+(pC *b))*km)+(pC*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\     p12=pG*b*w/(x+(b*w))\n\
\     p13=(-1)*(((d1l-c3l-pT)*m2l)+(((((-1*c3l)-pT)*f2l)+(b*d1l)-(pT *b))*km)-(pT*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\     p20=(-1)*(((d2l-c0l-pA)*m2l)+(((((-1*c0l)-pA)*f1l)+(b*d2l)-(pA *b))*km)-(pA*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p21=pC*b*w/(x+(b*w))\n\
\     p22=(((d2l+c2l+pG)*m2l)+((((c2l+pG)*f1l)+(b*d2l)+(pG *b))*km)+(pG*b*f1l*k2l))/((x+(b*w))*(x+(f1l*w)))\n\
\     p23=pT*b*w/(x+(b*w))\n\
\     p30=pA*b*w/(x+(b*w))\n\
\     p31=(-1)*(((d3l-c1l-pC)*m2l)+(((((-1*c1l)-pC)*f2l)+(b*d3l)-(pC *b))*km)-(pC*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\     p32=pG*b*w/(x+(b*w))\n\
\     p33=(((d3l+c3l+pT)*m2l)+((((c3l+pT)*f2l)+(b*d3l)+(pT *b))*km)+(pT*b*f2l*k2l))/((x+(b*w))*(x+(f2l*w)))\n\
\ in d8 k [[p00,p01,p02,p03],[p10,p11,p12,p13],[p20,p21,p22,p23],[p30,p31,p32,p33]]\n"

-- \hky85UniformWithK :: [Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]\n\
hky85UniformWithKString :: String
hky85UniformWithKString="\
\j7 [h, b] [pA,pC,pG,pT] x i (w,k)=\n\
\ if w<" ++ epsilonString ++ " then [[k,0,0,0],[0,k,0,0],[0,0,k,0],[0,0,0,k]]\n\
\ else\n\
\ let pR=pA+pG\n\
\     pY=pC+pT\n\
\     f1l=(pY*h+pR*b)\n\
\     f2l=(pR*h+pY*b)\n\
\     c0l=pA*pY/pR\n\
\     c1l=pC*pR/pY\n\
\     c2l=pG*pY/pR\n\
\     c3l=pT*pR/pY\n\
\     d0l=pG/pR\n\
\     d1l=pT/pY\n\
\     d2l=pA/pR\n\
\     d3l=pC/pY\n\
\     ebkm=e0(b*w*x)i 0\n\
\     enbkm=e0(-1*b*w*x)i 0\n\
\     enf1lkm=e0(-1*f1l*w*x)i 0\n\
\     enf2lkm=e0(-1*f2l*w*x)i 0\n\
\     bkmfactor=(enbkm+((b*w*x)-1))/(b*w*x)\n\
\     p00=((enbkm*((ebkm*((pA*b*f1l*w*x)+(c0l*f1l)+(b*d0l)))-(c0l*f1l)))-(b*d0l*enf1lkm))/(b*f1l*w*x)\n\
\     p01=pC*bkmfactor\n\
\     p02=((enbkm*((ebkm*((pG*b*f1l*w*x)+(c2l*f1l)-(b*d0l)))-(c2l*f1l)))+(b*d0l*enf1lkm))/(b*f1l*w*x)\n\
\     p03=pT*bkmfactor\n\
\     p10=pA*bkmfactor\n\
\     p11=((enbkm*((ebkm*((pC*b*f2l*w*x)+(c1l*f2l)+(b*d1l)))-(c1l*f2l)))-(b*d1l*enf2lkm))/(b*f2l*w*x)\n\
\     p12=pG*bkmfactor\n\
\     p13=((enbkm*((ebkm*((pT*b*f2l*w*x)+(c3l*f2l)-(b*d1l)))-(c3l*f2l)))+(b*d1l*enf2lkm))/(b*f2l*w*x)\n\
\     p20=((enbkm*((ebkm*((pA*b*f1l*w*x)+(c0l*f1l)-(b*d2l)))-(c0l*f1l)))+(b*d2l*enf1lkm))/(b*f1l*w*x)\n\
\     p21=pC*bkmfactor\n\
\     p22=((enbkm*((ebkm*((pG*b*f1l*w*x)+(c2l*f1l)+(b*d2l)))-(c2l*f1l)))-(b*d2l*enf1lkm))/(b*f1l*w*x)\n\
\     p23=pG*bkmfactor\n\
\     p30=pA*bkmfactor\n\
\     p31=((enbkm*((ebkm*((pC*b*f2l*w*x)+(c1l*f2l)-(b*d3l)))-(c1l*f2l)))+(b*d3l*enf2lkm))/(b*f2l*w*x)\n\
\     p32=pG*bkmfactor\n\
\     p33=((enbkm*((ebkm*((pT*b*f2l*w*x)+(c3l*f2l)+(b*d3l)))-(c3l*f2l)))-(b*d3l*enf2lkm))/(b*f2l*w*x)\n\
\ in d8 k [[p00,p01,p02,p03],[p10,p11,p12,p13],[p20,p21,p22,p23],[p30,p31,p32,p33]]\n"

-- \f81ExponentialWithK :: [Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]\n\
f81ExponentialWithKString :: String
f81ExponentialWithKString="\
\j8 z [pA,pC,pG,pT] x i (w,k)=\n\
\ if w<" ++ epsilonString ++ " then [[k,0,0,0],[0,k,0,0],[0,0,k,0],[0,0,0,k]]\n\
\ else\n\
\ let pR=pA+pG\n\
\     pY=pC+pT\n\
\     c0l=pA*pY/pR\n\
\     c1l=pC*pR/pY\n\
\     c2l=pG*pY/pR\n\
\     c3l=pT*pR/pY\n\
\     d0l=pG/pR\n\
\     d1l=pT/pY\n\
\     d2l=pA/pR\n\
\     d3l=pC/pY\n\
\     m2l=x*x\n\
\     km=x*w\n\
\     k2l=w*w\n\
\     p00=(((d0l+c0l+pA)*m2l)+((((c0l+pA))+d0l+pA)*km)+(pA*k2l))/((x+w)*(x+w))\n\
\     p01=pC*w/(x+w)\n\
\     p02=(-1)*(((d0l-c2l-pG)*m2l)+(((((-1*c2l)-pG))+d0l-pG)*km)-(pG*k2l))/((x+w)*(x+w))\n\
\     p03=pT*w/(x+w)\n\
\     p10=pA*w/(x+w)\n\
\     p11=(((d1l+c1l+pC)*m2l)+((((c1l+pC))+d1l+pC)*km)+(pC*k2l))/((x+w)*(x+w))\n\
\     p12=pG*w/(x+w)\n\
\     p13=(-1)*(((d1l-c3l-pT)*m2l)+(((((-1*c3l)-pT))+d1l-pT)*km)-(pT*k2l))/((x+w)*(x+w))\n\
\     p20=(-1)*(((d2l-c0l-pA)*m2l)+(((((-1*c0l)-pA))+(d2l)-pA)*km)-(pA*k2l))/((x+w)*(x+w))\n\
\     p21=pC*w/(x+w)\n\
\     p22=(((d2l+c2l+pG)*m2l)+((((c2l+pG))+(d2l)+pG)*km)+(pG*k2l))/((x+w)*(x+w))\n\
\     p23=pT*w/(x+w)\n\
\     p30=pA*w/(x+w)\n\
\     p31=(-1)*(((d3l-c1l-pC)*m2l)+(((((-1*c1l)-pC))+(d3l)-pC)*km)-(pC*k2l))/((x+w)*(x+w))\n\
\     p32=pG*w/(x+w)\n\
\     p33=(((d3l+c3l+pT)*m2l)+((((c3l+pT))+(d3l)+pT)*km)+(pT*k2l))/((x+w)*(x+w))\n\
\ in d8 k [[p00,p01,p02,p03],[p10,p11,p12,p13],[p20,p21,p22,p23],[p30,p31,p32,p33]]\n"

-- \f81UniformWithK :: [Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]\n\
f81UniformWithKString :: String
f81UniformWithKString="\
\j9 z [pA,pC,pG,pT] x i (w,k)=\n\
\ if w<" ++ epsilonString ++ " then [[k,0,0,0],[0,k,0,0],[0,0,k,0],[0,0,0,k]]\n\
\ else\n\
\ let pR=pA+pG\n\
\     pY=pC+pT\n\
\     c0l=pA*pY/pR\n\
\     c1l=pC*pR/pY\n\
\     c2l=pG*pY/pR\n\
\     c3l=pT*pR/pY\n\
\     d0l=pG/pR\n\
\     d1l=pT/pY\n\
\     d2l=pA/pR\n\
\     d3l=pC/pY\n\
\     ebkm=e0(w*x)i 0\n\
\     enbkm=e0(-1*w*x)i 0\n\
\     enf1lkm=e0(-1*w*x)i 0\n\
\     enf2lkm=e0(-1*w*x)i 0\n\
\     bkmfactor=(enbkm+((w*x)-1))/(w*x)\n\
\     p00=((enbkm*((ebkm*((pA*w*x)+(c0l)+d0l))-(c0l)))-(d0l*enf1lkm))/(w*x)\n\
\     p01=pC*bkmfactor\n\
\     p02=((enbkm*((ebkm*((pG*w*x)+(c2l)-d0l))-(c2l)))+(d0l*enf1lkm))/(w*x)\n\
\     p03=pT*bkmfactor\n\
\     p10=pA*bkmfactor\n\
\     p11=((enbkm*((ebkm*((pC*w*x)+(c1l)+d1l))-(c1l)))-(d1l*enf2lkm))/(w*x)\n\
\     p12=pG*bkmfactor\n\
\     p13=((enbkm*((ebkm*((pT*w*x)+(c3l)-d1l))-(c3l)))+(d1l*enf2lkm))/(w*x)\n\
\     p20=((enbkm*((ebkm*((pA*w*x)+(c0l)-(d2l)))-(c0l)))+(d2l*enf1lkm))/(w*x)\n\
\     p21=pC*bkmfactor\n\
\     p22=((enbkm*((ebkm*((pG*w*x)+(c2l)+(d2l)))-(c2l)))-(d2l*enf1lkm))/(w*x)\n\
\     p23=pG*bkmfactor\n\
\     p30=pA*bkmfactor\n\
\     p31=((enbkm*((ebkm*((pC*w*x)+(c1l)-(d3l)))-(c1l)))+(d3l*enf2lkm))/(w*x)\n\
\     p32=pG*bkmfactor\n\
\     p33=((enbkm*((ebkm*((pT*w*x)+(c3l)+(d3l)))-(c3l)))-(d3l*enf2lkm))/(w*x)\n\
\ in d8 k [[p00,p01,p02,p03],[p10,p11,p12,p13],[p20,p21,p22,p23],[p30,p31,p32,p33]]\n"

-- \k80ExponentialWithK :: [Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]\n\
k80ExponentialWithKString :: String
k80ExponentialWithKString="\
\k0 [h, b] z x i (w,k)=\n\
\ if w<" ++ epsilonString ++ " then [[k,0,0,0],[0,k,0,0],[0,0,k,0],[0,0,0,k]]\n\
\ else\n\
\ let f=(h +b)/2\n\
\     c=0.25\n\
\     d=0.5\n\
\     m2l=x*x\n\
\     km=x*w\n\
\     k2l=w*w\n\
\     o=0.25*b*w/(x+(b*w))\n\
\     p=(((d+c+0.25)*m2l)+((((c+0.25)*f)+(b*d)+(0.25 *b))*km)+(0.25*b*f*k2l))/((x+(b*w))*(x+(f*w)))\n\
\     m=(-1)*(((d-c-0.25)*m2l)+(((((-1*c)-0.25)*f)+(b*d)-(0.25*b))*km)-(0.25*b*f*k2l))/((x+(b*w))*(x+(f*w)))\n\
\ in d8 k [[p,o,m,o],[o,p,o,m],[m,o,p,o],[o,m,o,p]]\n"

-- \k80UniformWithK :: [Double] -> [Double] -> Double -> Int -> (Double, Double) -> [[Double]]\n\
k80UniformWithKString :: String
k80UniformWithKString="\
\k1 [h, b] z x i (w,k)=\n\
\ if w<" ++ epsilonString ++ " then [[k,0,0,0],[0,k,0,0],[0,0,k,0],[0,0,0,k]]\n\
\ else\n\
\ let f=(h +b)/2\n\
\     c=0.25\n\
\     d=0.5\n\
\     ebkm=e0(b*w*x)i 0\n\
\     enbkm=e0(-1*b*w*x)i 0\n\
\     enfkm=e0(-1*f*w*x)i 0\n\
\     o=0.25*(enbkm+((b*w*x)-1))/(b*w*x)\n\
\     p=((enbkm*((ebkm*((0.25*b*f*w*x)+(c*f)+(b*d)))-(c*f)))-(b*d*enfkm))/(b*f*w*x)\n\
\     m=((enbkm*((ebkm*((0.25*b*f*w*x)+(c*f)-(b*d)))-(c*f)))+(b*d*enfkm))/(b*f*w*x)\n\
\ in d8 k [[p,o,m,o],[o,p,o,m],[m,o,p,o],[o,m,o,p]]\n"

