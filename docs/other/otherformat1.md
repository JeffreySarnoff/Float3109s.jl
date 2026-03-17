Let

*   KKK = total bitwidth
    
*   PPP = precision in bits, counting the implicit leading 1 for normals
    
*   W\=K−PW = K-PW\=K−P = exponent-field width
    
*   explicit fraction width F\=P−1F = P-1F\=P−1
    

and assume

*   one zero code point at 000
    
*   one NaN code point at 2K−12^K-12K−1
    
*   all other codes are finite nonzero numbers
    
*   sign symmetry holds
    
*   the represented values are dyadic
    
*   there are equally many magnitudes below and above 111
    

This determines a surprisingly large amount.

* * *

## 1\. Basic parameter constraints

Because there is a sign bit and K−PK-PK−P exponent bits, the remaining explicit significand bits are

F\=P−1.F = P-1.F\=P−1.

So the total payload structure is

1 sign bit+W exponent bits+F fraction bits\=K,1 \\text{ sign bit} + W \\text{ exponent bits} + F \\text{ fraction bits} = K,1 sign bit+W exponent bits+F fraction bits\=K,

with

W\=K−P,F\=P−1.W = K-P, \\qquad F=P-1.W\=K−P,F\=P−1.

Hence necessarily

1≤P≤K−1.1 \\le P \\le K-1.1≤P≤K−1.

The case P\=1P=1P\=1 is the degenerate “powers of two only” case: there are no explicit fraction bits and no subnormals.

* * *

## 2\. Code-point counts

There are 2K2^K2K total code points.

Two are special:

*   zero at code point 000
    
*   NaN at code point 2K−12^K-12K−1
    

So the number of finite nonzero codes is

2K−2.2^K-2.2K−2.

By sign symmetry, these split evenly into positive and negative:

#{positive nonzero}\=#{negative nonzero}\=2K−1−1.\\#\\{\\text{positive nonzero}\\}=\\#\\{\\text{negative nonzero}\\}=2^{K-1}-1.#{positive nonzero}\=#{negative nonzero}\=2K−1−1.

Thus the positive nonzero magnitudes are in one-to-one correspondence with the integers

q\=1,2,…,2K−1−1.q = 1,2,\\dots,2^{K-1}-1.q\=1,2,…,2K−1−1.

This is the natural positive-code index.

A convenient canonical code-point assignment is:

*   c(0)\=0c(0)=0c(0)\=0
    
*   for positive nonzero values, c(x)\=qc(x)=qc(x)\=q, where q∈\[1,2K−1−1\]q\\in\[1,2^{K-1}-1\]q∈\[1,2K−1−1\]
    
*   for negative values, c(−x)\=(2K−1−1)+qc(-x)= (2^{K-1}-1)+qc(−x)\=(2K−1−1)+q
    
*   c(NaN)\=2K−1c(\\text{NaN})=2^K-1c(NaN)\=2K−1
    

So negative finite codes occupy

2K−1,…,2K−2.2^{K-1},\\dots,2^K-2.2K−1,…,2K−2.

* * *

## 3\. Positive-code decomposition

For each positive nonzero code index q∈{1,…,2K−1−1}q\\in\\{1,\\dots,2^{K-1}-1\\}q∈{1,…,2K−1−1}, write

q\=e 2P−1+tq = e\\,2^{P-1}+tq\=e2P−1+t

with

0≤e≤2W−1,0≤t≤2P−1−1,0 \\le e \\le 2^W-1, \\qquad 0 \\le t \\le 2^{P-1}-1,0≤e≤2W−1,0≤t≤2P−1−1,

and not both e\=t\=0e=t=0e\=t\=0.

Then:

*   eee is the exponent-field value
    
*   ttt is the trailing-significand field
    

This works because

(2W−1)2P−1+(2P−1−1)\=2W+P−1−1\=2K−1−1.(2^W-1)2^{P-1} + (2^{P-1}-1) = 2^{W+P-1}-1 = 2^{K-1}-1.(2W−1)2P−1+(2P−1−1)\=2W+P−1−1\=2K−1−1.

So the positive nonzero codes are exactly exhausted by

*   subnormals: e\=0, t\=1,…,2P−1−1e=0,\\ t=1,\\dots,2^{P-1}-1e\=0, t\=1,…,2P−1−1
    
*   normals: e\=1,…,2W−1, t\=0,…,2P−1−1e=1,\\dots,2^W-1,\\ t=0,\\dots,2^{P-1}-1e\=1,…,2W−1, t\=0,…,2P−1−1
    

Counts:

*   subnormals: 2P−1−12^{P-1}-12P−1−1
    
*   normals: (2W−1)2P−1(2^W-1)2^{P-1}(2W−1)2P−1
    

Total:

(2P−1−1)+(2W−1)2P−1\=2K−1−1.(2^{P-1}-1)+(2^W-1)2^{P-1}=2^{K-1}-1.(2P−1−1)+(2W−1)2P−1\=2K−1−1.

* * *

## 4\. The exponent bias is forced

Let the bias be BBB. Then normal values decode as

x\=(1+t2P−1)2e−B,e≥1,x = \\left(1+\\frac{t}{2^{P-1}}\\right)2^{e-B}, \\qquad e\\ge 1,x\=(1+2P−1t​)2e−B,e≥1,

and subnormals decode as

x\=(t2P−1)21−B,e\=0, t≥1.x = \\left(\\frac{t}{2^{P-1}}\\right)2^{1-B}, \\qquad e=0,\\ t\\ge 1.x\=(2P−1t​)21−B,e\=0, t≥1.

Now impose:

> the number of magnitudes below 111 equals the number above 111.

The value 111 occurs when

e\=B,t\=0.e=B,\\quad t=0.e\=B,t\=0.

Count positive magnitudes <1<1<1:

*   subnormals: 2P−1−12^{P-1}-12P−1−1
    
*   normal binades with e\=1,…,B−1e=1,\\dots,B-1e\=1,…,B−1: (B−1)2P−1(B-1)2^{P-1}(B−1)2P−1
    

So

N<1\=(B−1)2P−1+(2P−1−1)\=B 2P−1−1.N\_{<1}= (B-1)2^{P-1}+(2^{P-1}-1)=B\\,2^{P-1}-1.N<1​\=(B−1)2P−1+(2P−1−1)\=B2P−1−1.

Count positive magnitudes \>1\>1\>1:

*   in the e\=Be=Be\=B binade, all except t\=0t=0t\=0: 2P−1−12^{P-1}-12P−1−1
    
*   binades e\=B+1,…,2W−1e=B+1,\\dots,2^W-1e\=B+1,…,2W−1: (2W−1−B)2P−1(2^W-1-B)2^{P-1}(2W−1−B)2P−1
    

So

N\>1\=(2W−1−B)2P−1+(2P−1−1)\=(2W−B)2P−1−1.N\_{>1}=(2^W-1-B)2^{P-1}+(2^{P-1}-1)=(2^W-B)2^{P-1}-1.N\>1​\=(2W−1−B)2P−1+(2P−1−1)\=(2W−B)2P−1−1.

Equating them gives

B\=2W−1.B = 2^{W-1}.B\=2W−1.

So the bias is uniquely determined:

B\=2W−1\=2K−P−1.\\boxed{B=2^{W-1}=2^{K-P-1}}.B\=2W−1\=2K−P−1​.

This is not the IEEE-754 usual bias 2W−1−12^{W-1}-12W−1−1; your symmetry condition forces the shifted value.

* * *

## 5\. Exponent range

Since normal exponent field values are e\=1,…,2W−1e=1,\\dots,2^W-1e\=1,…,2W−1, the unbiased exponent range is

E\=e−B.E = e-B.E\=e−B.

Therefore

Emin⁡\=1−B\=1−2W−1,E\_{\\min} = 1-B = 1-2^{W-1},Emin​\=1−B\=1−2W−1, Emax⁡\=(2W−1)−B\=2W−1−1.E\_{\\max} = (2^W-1)-B = 2^{W-1}-1.Emax​\=(2W−1)−B\=2W−1−1.

So

Emin⁡\=1−2W−1,Emax⁡\=2W−1−1.\\boxed{E\_{\\min}=1-2^{W-1}},\\qquad \\boxed{E\_{\\max}=2^{W-1}-1}.Emin​\=1−2W−1​,Emax​\=2W−1−1​.

These are symmetric around 000 in the sense that

Emin⁡+Emax⁡\=0.E\_{\\min}+E\_{\\max}=0.Emin​+Emax​\=0.

There are

2W−12^W-12W−1

normal exponent levels.

* * *

## 6\. Exact decode formulas

## Positive values

For q∈{1,…,2K−1−1}q\\in\\{1,\\dots,2^{K-1}-1\\}q∈{1,…,2K−1−1}, write

q\=e 2P−1+t,0≤e≤2W−1,0≤t≤2P−1−1.q=e\\,2^{P-1}+t,\\qquad 0\\le e\\le 2^W-1,\\quad 0\\le t\\le 2^{P-1}-1.q\=e2P−1+t,0≤e≤2W−1,0≤t≤2P−1−1.

Then

### Subnormal case

If e\=0e=0e\=0, necessarily 1≤t≤2P−1−11\\le t\\le 2^{P-1}-11≤t≤2P−1−1, and

x(q)\=t 2Emin⁡−(P−1)\\boxed{x(q)= t\\,2^{E\_{\\min}-(P-1)}}x(q)\=t2Emin​−(P−1)​

since

(t2P−1)2Emin⁡\=t 2Emin⁡−(P−1).\\left(\\frac{t}{2^{P-1}}\\right)2^{E\_{\\min}} = t\\,2^{E\_{\\min}-(P-1)}.(2P−1t​)2Emin​\=t2Emin​−(P−1).

### Normal case

If e≥1e\\ge 1e≥1, then

x(q)\=(2P−1+t)2e−B−(P−1).\\boxed{x(q)=\\left(2^{P-1}+t\\right)2^{e-B-(P-1)}}.x(q)\=(2P−1+t)2e−B−(P−1)​.

Equivalently, if

M\=2P−1+t∈{2P−1,…,2P−1},M = 2^{P-1}+t \\in \\{2^{P-1},\\dots,2^P-1\\},M\=2P−1+t∈{2P−1,…,2P−1},

then

x\=M 2E−(P−1),E\=e−B.\\boxed{x = M\\,2^{E-(P-1)}},\\qquad E=e-B.x\=M2E−(P−1)​,E\=e−B.

So every normal value has an integer PPP\-bit significand MMM.

## Negative values

If qqq is the positive index of ∣x∣|x|∣x∣, then the negative code point is

c(−x)\=(2K−1−1)+qc(-x)= (2^{K-1}-1)+qc(−x)\=(2K−1−1)+q

and the value is

x(c)\=−x(q).\\boxed{x(c)=-x(q)}.x(c)\=−x(q)​.

* * *

## 7\. Exact encode formulas

Let x\>0x>0x\>0 be representable.

### Subnormal encode

If

0<x<2Emin⁡,0 < x < 2^{E\_{\\min}},0<x<2Emin​,

then xxx must be of the form

x\=t 2Emin⁡−(P−1)x = t\\,2^{E\_{\\min}-(P-1)}x\=t2Emin​−(P−1)

for an integer

1≤t≤2P−1−1.1 \\le t \\le 2^{P-1}-1.1≤t≤2P−1−1.

Then

q\=t.\\boxed{q=t}.q\=t​.

### Normal encode

If x≥2Emin⁡x\\ge 2^{E\_{\\min}}x≥2Emin​, then it must be of the form

x\=M 2E−(P−1)x=M\\,2^{E-(P-1)}x\=M2E−(P−1)

with

2P−1≤M≤2P−1,Emin⁡≤E≤Emax⁡.2^{P-1}\\le M\\le 2^P-1,\\qquad E\_{\\min}\\le E\\le E\_{\\max}.2P−1≤M≤2P−1,Emin​≤E≤Emax​.

Then

e\=B+E,t\=M−2P−1,e = B+E, \\qquad t=M-2^{P-1},e\=B+E,t\=M−2P−1,

and

q\=(B+E)2P−1+(M−2P−1).\\boxed{q=(B+E)2^{P-1} + (M-2^{P-1})}.q\=(B+E)2P−1+(M−2P−1)​.

For a negative representable value x<0x<0x<0,

c(x)\=(2K−1−1)+q(∣x∣).\\boxed{c(x)=(2^{K-1}-1)+q(|x|)}.c(x)\=(2K−1−1)+q(∣x∣)​.

* * *

## 8\. Special values and landmark code points

Because B\=2W−1B=2^{W-1}B\=2W−1,

q(1)\=B 2P−1\=2W−12P−1\=2K−2.q(1)=B\\,2^{P-1}=2^{W-1}2^{P-1}=2^{K-2}.q(1)\=B2P−1\=2W−12P−1\=2K−2.

So

c(+1)\=2K−2.\\boxed{c(+1)=2^{K-2}}.c(+1)\=2K−2​.

And therefore

c(−1)\=2K−1−1+2K−2\=3⋅2K−2−1.\\boxed{c(-1)=2^{K-1}-1+2^{K-2}=3\\cdot 2^{K-2}-1}.c(−1)\=2K−1−1+2K−2\=3⋅2K−2−1​.

More generally, for exact powers of two 2n2^n2n in range Emin⁡≤n≤Emax⁡E\_{\\min}\\le n\\le E\_{\\max}Emin​≤n≤Emax​,

q(2n)\=2K−2+n 2P−1.\\boxed{q(2^n)=2^{K-2}+n\\,2^{P-1}}.q(2n)\=2K−2+n2P−1​.

Hence

q(1/2)\=2K−2−2P−1,q(2)\=2K−2+2P−1.\\boxed{q(1/2)=2^{K-2}-2^{P-1}}, \\qquad \\boxed{q(2)=2^{K-2}+2^{P-1}}.q(1/2)\=2K−2−2P−1​,q(2)\=2K−2+2P−1​.

And for negatives:

c(−2n)\=2K−1−1+2K−2+n 2P−1.\\boxed{c(-2^n)=2^{K-1}-1+2^{K-2}+n\\,2^{P-1}}.c(−2n)\=2K−1−1+2K−2+n2P−1​.

In particular,

c(−1/2)\=2K−1−1+2K−2−2P−1,\\boxed{c(-1/2)=2^{K-1}-1+2^{K-2}-2^{P-1}},c(−1/2)\=2K−1−1+2K−2−2P−1​, c(−2)\=2K−1−1+2K−2+2P−1.\\boxed{c(-2)=2^{K-1}-1+2^{K-2}+2^{P-1}}.c(−2)\=2K−1−1+2K−2+2P−1​.

Also:

*   smallest positive subnormal:
    
    xmin⁡+\=2Emin⁡−(P−1)\=2 2−P−2W−1\\boxed{x\_{\\min+}=2^{E\_{\\min}-(P-1)} = 2^{\\,2-P-2^{W-1}}}xmin+​\=2Emin​−(P−1)\=22−P−2W−1​
    
    at code point 111
    
*   largest subnormal:
    
    xsub,max\=(2P−1−1)2Emin⁡−(P−1)\\boxed{x\_{\\text{sub,max}}=(2^{P-1}-1)2^{E\_{\\min}-(P-1)}}xsub,max​\=(2P−1−1)2Emin​−(P−1)​
    
    at code point 2P−1−12^{P-1}-12P−1−1
    
*   smallest positive normal:
    
    xnorm,min\=2Emin⁡\\boxed{x\_{\\text{norm,min}}=2^{E\_{\\min}}}xnorm,min​\=2Emin​​
    
    at code point 2P−12^{P-1}2P−1
    
*   largest finite positive:
    
    xmax⁡\=(2−2−(P−1))2Emax⁡\\boxed{x\_{\\max}=(2-2^{-(P-1)})2^{E\_{\\max}}}xmax​\=(2−2−(P−1))2Emax​​
    
    at code point 2K−1−12^{K-1}-12K−1−1
    

* * *

## 9\. Value set structure

The positive values are exactly

{ t 2Emin⁡−(P−1):1≤t≤2P−1−1 }\\{\\,t\\,2^{E\_{\\min}-(P-1)} : 1\\le t\\le 2^{P-1}-1\\,\\}{t2Emin​−(P−1):1≤t≤2P−1−1}

for subnormals, together with

{ M 2E−(P−1):Emin⁡≤E≤Emax⁡, 2P−1≤M≤2P−1 }\\{\\,M\\,2^{E-(P-1)} : E\_{\\min}\\le E\\le E\_{\\max},\\ 2^{P-1}\\le M\\le 2^P-1\\,\\}{M2E−(P−1):Emin​≤E≤Emax​, 2P−1≤M≤2P−1}

for normals.

So the full finite value set is

{0} ∪ {±x:x∈V+}\\boxed{ \\{0\\}\\ \\cup\\ \\{\\pm x : x\\in\\mathcal V\_+\\} }{0} ∪ {±x:x∈V+​}​

where V+\\mathcal V\_+V+​ is the positive set above.

Every finite nonzero value is dyadic, i.e. of the form

±m 2r\\pm m\\,2^r±m2r

with integer mmm.

* * *

## 10\. Binades and lattice structure

Each normal exponent EEE gives one binade

\[2E, (2−2−(P−1))2E\]\[2^E,\\ (2-2^{-(P-1)})2^E\]\[2E, (2−2−(P−1))2E\]

containing exactly

2P−12^{P-1}2P−1

values:

2E, (1+2−(P−1))2E, (1+2⋅2−(P−1))2E, …, (2−2−(P−1))2E.2^E,\\ \\left(1+2^{-(P-1)}\\right)2^E,\\ \\left(1+2\\cdot 2^{-(P-1)}\\right)2^E,\\ \\dots,\\ (2-2^{-(P-1)})2^E.2E, (1+2−(P−1))2E, (1+2⋅2−(P−1))2E, …, (2−2−(P−1))2E.

The subnormal region is

{ 2Emin⁡−(P−1), 2⋅2Emin⁡−(P−1), …, (2P−1−1)2Emin⁡−(P−1) }.\\{\\,2^{E\_{\\min}-(P-1)},\\ 2\\cdot 2^{E\_{\\min}-(P-1)},\\ \\dots,\\ (2^{P-1}-1)2^{E\_{\\min}-(P-1)}\\,\\}.{2Emin​−(P−1), 2⋅2Emin​−(P−1), …, (2P−1−1)2Emin​−(P−1)}.

Thus:

*   number of normal binades:
    
    2W−1\\boxed{2^W-1}2W−1​
*   number of values per normal binade:
    
    2P−1\\boxed{2^{P-1}}2P−1​
*   number of positive subnormals:
    
    2P−1−1\\boxed{2^{P-1}-1}2P−1−1​

* * *

## 11\. Spacing, ulps, and adjacency

Within a normal binade of exponent EEE, adjacent values differ by

ulp⁡(E)\=2E−(P−1).\\boxed{\\operatorname{ulp}(E)=2^{E-(P-1)}}.ulp(E)\=2E−(P−1)​.

In the subnormal region, adjacent values also differ by

2Emin⁡−(P−1).\\boxed{2^{E\_{\\min}-(P-1)}}.2Emin​−(P−1)​.

The transition from largest subnormal to smallest normal also has exactly that gap:

2Emin⁡−(2P−1−1)2Emin⁡−(P−1)\=2Emin⁡−(P−1).2^{E\_{\\min}}-(2^{P-1}-1)2^{E\_{\\min}-(P-1)} =2^{E\_{\\min}-(P-1)}.2Emin​−(2P−1−1)2Emin​−(P−1)\=2Emin​−(P−1).

Likewise, the jump from the largest value in binade EEE to the first in binade E+1E+1E+1 is

2E+1−(2−2−(P−1))2E\=2E−(P−1).2^{E+1}-(2-2^{-(P-1)})2^E =2^{E-(P-1)}.2E+1−(2−2−(P−1))2E\=2E−(P−1).

So the code-index successor function on positive values is very regular:

*   increasing qqq by 111 always moves to the next positive representable value
    
*   the step size doubles exactly when the exponent binade increments
    

* * *

## 12\. Symmetry around 1

Since

c(+1)\=2K−2,c(+1)=2^{K-2},c(+1)\=2K−2,

there are exactly

2K−2−12^{K-2}-12K−2−1

positive magnitudes below 111, and exactly the same number above 111.

So among positive nonzero values,

#{x:0<x<1}\=#{x:x\>1}\=2K−2−1.\\#\\{x:0<x<1\\}=\\#\\{x:x>1\\}=2^{K-2}-1.#{x:0<x<1}\=#{x:x\>1}\=2K−2−1.

This is equivalent to your balance condition and is exactly why the bias is 2W−12^{W-1}2W−1.

* * *

## 13\. Dynamic range

The smallest positive value is

xmin⁡+\=2Emin⁡−(P−1).x\_{\\min+}=2^{E\_{\\min}-(P-1)}.xmin+​\=2Emin​−(P−1).

The largest positive finite value is

xmax⁡\=(2−2−(P−1))2Emax⁡.x\_{\\max}=(2-2^{-(P-1)})2^{E\_{\\max}}.xmax​\=(2−2−(P−1))2Emax​.

Their ratio is

xmax⁡xmin⁡+\=(2−2−(P−1))2Emax⁡−Emin⁡+P−1.\\frac{x\_{\\max}}{x\_{\\min+}} = (2-2^{-(P-1)})2^{E\_{\\max}-E\_{\\min}+P-1}.xmin+​xmax​​\=(2−2−(P−1))2Emax​−Emin​+P−1.

Since

Emax⁡−Emin⁡\=(2W−1−1)−(1−2W−1)\=2W−2,E\_{\\max}-E\_{\\min} = (2^{W-1}-1)-(1-2^{W-1}) = 2^W-2,Emax​−Emin​\=(2W−1−1)−(1−2W−1)\=2W−2,

this becomes

xmax⁡xmin⁡+\=(2P−1) 2 2W−2.\\boxed{ \\frac{x\_{\\max}}{x\_{\\min+}} = (2^P-1)\\,2^{\\,2^W-2} }.xmin+​xmax​​\=(2P−1)22W−2​.

That is the exact positive dynamic range factor.

* * *

## 14\. Representability criterion

A positive dyadic number xxx is representable iff exactly one of the following holds.

### Subnormal representability

0<x<2Emin⁡0 < x < 2^{E\_{\\min}}0<x<2Emin​

and

x 2P−1−Emin⁡x\\,2^{P-1-E\_{\\min}}x2P−1−Emin​

is an integer in

{1,…,2P−1−1}.\\{1,\\dots,2^{P-1}-1\\}.{1,…,2P−1−1}.

### Normal representability

There exists an integer EEE with

Emin⁡≤E≤Emax⁡E\_{\\min}\\le E\\le E\_{\\max}Emin​≤E≤Emax​

such that

x 2P−1−Ex\\,2^{P-1-E}x2P−1−E

is an integer in

{2P−1,…,2P−1}.\\{2^{P-1},\\dots,2^P-1\\}.{2P−1,…,2P−1}.

Equivalently: the normalized binary significand of xxx has at most PPP bits, and the unbiased exponent lies in range.

* * *

## 15\. The P\=1P=1P\=1 special case

When P\=1P=1P\=1,

*   there are no explicit fraction bits
    
*   there are no subnormals
    
*   every positive finite value is an exact power of two
    

Then

W\=K−1,B\=2K−2,W=K-1,\\qquad B=2^{K-2},W\=K−1,B\=2K−2, Emin⁡\=1−2K−2,Emax⁡\=2K−2−1.E\_{\\min}=1-2^{K-2},\\qquad E\_{\\max}=2^{K-2}-1.Emin​\=1−2K−2,Emax​\=2K−2−1.

Positive values are exactly

2Emin⁡,2Emin⁡+1,…,2Emax⁡.2^{E\_{\\min}}, 2^{E\_{\\min}+1}, \\dots, 2^{E\_{\\max}}.2Emin​,2Emin​+1,…,2Emax​.

There are

2K−1−12^{K-1}-12K−1−1

of them, as required.

So P\=1P=1P\=1 is the pure logarithmic signed-power-of-two member of the family.

* * *

## 16\. Clean summary of the core formulas

Let

W\=K−P,B\=2W−1.W=K-P,\\qquad B=2^{W-1}.W\=K−P,B\=2W−1.

Then:

Emin⁡\=1−B,Emax⁡\=B−1.\\boxed{E\_{\\min}=1-B},\\qquad \\boxed{E\_{\\max}=B-1}.Emin​\=1−B​,Emax​\=B−1​.

Positive code indices are

q\=1,…,2K−1−1,q\=e 2P−1+t.q=1,\\dots,2^{K-1}-1, \\qquad q=e\\,2^{P-1}+t.q\=1,…,2K−1−1,q\=e2P−1+t.

Decode:

x(q)\={t 2Emin⁡−(P−1),e\=0, t≥1,(2P−1+t) 2e−B−(P−1),e≥1.\\boxed{ x(q)= \\begin{cases} t\\,2^{E\_{\\min}-(P-1)}, & e=0,\\ t\\ge 1,\\\\\[1ex\] (2^{P-1}+t)\\,2^{e-B-(P-1)}, & e\\ge 1. \\end{cases}}x(q)\={t2Emin​−(P−1),(2P−1+t)2e−B−(P−1),​e\=0, t≥1,e≥1.​​

Negative code points:

c(−x)\=2K−1−1+c(∣x∣).\\boxed{c(-x)=2^{K-1}-1+c(|x|)}.c(−x)\=2K−1−1+c(∣x∣)​.

Special code points:

c(0)\=0,c(NaN)\=2K−1,\\boxed{c(0)=0},\\qquad \\boxed{c(\\text{NaN})=2^K-1},c(0)\=0​,c(NaN)\=2K−1​, c(+1)\=2K−2,c(−1)\=3⋅2K−2−1.\\boxed{c(+1)=2^{K-2}},\\qquad \\boxed{c(-1)=3\\cdot 2^{K-2}-1}.c(+1)\=2K−2​,c(−1)\=3⋅2K−2−1​.

Exact powers of two:

c(2n)\=2K−2+n 2P−1(Emin⁡≤n≤Emax⁡),\\boxed{c(2^n)=2^{K-2}+n\\,2^{P-1}} \\quad (E\_{\\min}\\le n\\le E\_{\\max}),c(2n)\=2K−2+n2P−1​(Emin​≤n≤Emax​),

and

c(−2n)\=2K−1−1+2K−2+n 2P−1.\\boxed{c(-2^n)=2^{K-1}-1+2^{K-2}+n\\,2^{P-1}}.c(−2n)\=2K−1−1+2K−2+n2P−1​.

Counts:

#positive nonzero\=2K−1−1,\\boxed{\\#\\text{positive nonzero}=2^{K-1}-1},#positive nonzero\=2K−1−1​, #subnormals\=2P−1−1,\\boxed{\\#\\text{subnormals}=2^{P-1}-1},#subnormals\=2P−1−1​, #normal binades\=2K−P−1,\\boxed{\\#\\text{normal binades}=2^{K-P}-1},#normal binades\=2K−P−1​, #values per normal binade\=2P−1.\\boxed{\\#\\text{values per normal binade}=2^{P-1}}.#values per normal binade\=2P−1​.

* * *

## 17\. What is determined uniquely by your axioms

Your assumptions force all of the following:

*   the exponent bias is not arbitrary; it must be
    
    2K−P−12^{K-P-1}2K−P−1
*   the positive value 111 sits exactly at code point
    
    2K−22^{K-2}2K−2
*   the normal exponent range is exactly
    
    1−2K−P−1 to 2K−P−1−11-2^{K-P-1}\\ \\text{to}\\ 2^{K-P-1}-11−2K−P−1 to 2K−P−1−1
*   the positive codes decompose exactly into one subnormal strip plus 2K−P−12^{K-P}-12K−P−1 equal-width normal binades
    
*   the code-to-value map on positive codes is strictly increasing
    
*   the spacing is uniform within each binade and doubles from one binade to the next
    
*   the transition gaps at subnormal→normal and binade→binade boundaries match the local ulp exactly
    

So, up to your chosen sign/codepoint convention, this family is essentially fully determined.