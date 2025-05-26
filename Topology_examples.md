# Some examples in topology

#### Continuous map that sends an open subset to a non-open subset

Let $f(x)=x^2$. $f[(-1,1)]=[0,1)$.

#### An image of the closure of a set is strictly contained within the closure of the image of that set under a continuous map ($f[{\rm Cl}(A)]\subsetneq {\rm Cl}[f(A)]$)


Let $f(x)=e^{-x}$ and $A:=(0,+\infty)$. 
$f[{\rm Cl}(A)]=(0,1]\subset [0,1]={\rm Cl}[f(A)]$.

#### Factorial map that is not open

Let $X=[0,1]\subset \mathbb{R}$ and let $x\sim y$, whenever $x=y$ or $\{x,y\}=\{0,1\}$ (0 and 1 glued together). Then $X/\sim$ is homeomorphic to a circle $S^1$. The corresponding factorial map $f: X\rightarrow S^1$ is not open: The image of an open subset $[0,0.1)\subset X$ is half-open interval in $S^1$, which is not open in the standard topology of $S^1$.