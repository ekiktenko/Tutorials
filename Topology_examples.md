# Some examples in topology

Here we condider some examples from topology that may provide additional insight (at least, they were helpful for the author :-)).

#### Continuous map that sends an open set to a non-open set

Let $f(x)=x^2$. $f[(-1,1)]=[0,1)$.

#### An image of the closure of a set is strictly contained in the closure of the image of that set under a continuous map: $f[{\rm Cl}(A)]\subsetneq {\rm Cl}[f(A)]$


Let $f(x)=e^{-x}$ and $A=(0,+\infty)$. 
$f[{\rm Cl}(A)]=(0,1]\subset [0,1]={\rm Cl}[f(A)]$.

#### Factorial map that is not open

Let $X=[0,1]\subset \mathbb{R}$ and let $x\sim y$, whenever $x=y$ or $\{x,y\}=\{0,1\}$ (0 and 1 glued together). Then $X/\sim$ is homeomorphic to a circle $S^1$. The corresponding factorial map $f: X\rightarrow S^1$ is not open: The image of an open subset $[0,0.1)\subset X$ is half-open interval in $S^1$, which is not open in the standard topology of $S^1$.

#### Discontinuous open map

Let $X=\{a,b,c\}$, $T_X=\{X,\{\emptyset\},\{a\}\}$, $Y=\{0,1\}$ and $f:X\rightarrow Y$ acts as
$$
    f(x) = \begin{cases}
        0, x=a,b\\
        1, x=c
    \end{cases}.
$$
Consider the topology $T_Y=\{f(U)\subset Y: U\in T_X\}=\{Y, \{\emptyset\}, \{0\}\}$.
By construction, $f$ is open on $(X, T_X)$ and $(Y, T_Y)$. 
However, it is not continuous, since  $f^{-1}(\{0\})=\{a,b\}\not\in T_X$.