# Notes

<!--toc:start-->
- [Notes](#notes)
  - [Univariate representation](#univariate-representation)
  - [Polynomials as evaluations](#polynomials-as-evaluations)
  - [Evaluation and interpolation](#evaluation-and-interpolation)
    - [Evaluation](#evaluation)
    - [Interpolation](#interpolation)
    - [Barycentric evaluation](#barycentric-evaluation)
  - [Optimization approach](#optimization-approach)
    - [Optimizing divisor arithmetic](#optimizing-divisor-arithmetic)
    - [Batching line computation](#batching-line-computation)
  - [Potential improvements](#potential-improvements)
    - [Within rules](#within-rules)
    - [Beyond rules](#beyond-rules)
<!--toc:end-->

In this notes I describe the main changes done to the original implementation,
together with some background in the polynomial arithmetic techniques used.  
Some additional optimizations which I didn't implement are also mentioned.

## Univariate representation

Divisors can be represented as 2 univariate polynomials $a$ and $b$ such that:
$$f(x,y) = a(x) - yb(x)$$  
That's also already the form used for ECIP, then multiplying divisors is just
univariate polynomial arithmetic.
$$f_1(x,y) = a_1(x) - yb_1(x)$$
$$f_2(x,y) = a_2(x) - yb_2(x)$$
$$ f_1 f_2 = (a_1 - yb_1)(a_2 - yb_2) $$
$$ f_1 f_2 = a_1(a_2 - yb_2) - yb_1(a_2 - yb_2)$$
$$ f_1 f_2 = a_1a_2 - ya_1b_2 - yb_1a_2 + y^2b_1b_2$$
$$ f_1 f_2 = a_1a_2 - y(a_1b_2 + b_1a_2) + y^2b_1b_2$$

Reduction is done by replacing $y^2$ by $x^3+Ax+B$ (from the curve equation).
Returning to the original $a - yb$ form.
$$ f_1 f_2 = a_1a_2 - y(a_1b_2 + b_1a_2) + (x^3+Ax+B)b_1b_2$$
$$ f_1 f_2 = (a_1a_2 + (x^3+Ax+B)b_1b_2) - y(a_1b_2 + b_1a_2) $$
$$ f_1 f_2 = a' - yb' $$

## Polynomials as evaluations

The most common form of representing a polynomial $p(x)$ of degree $d$ is as
$d+1$ coefficients, each a field element. With this representation addition is
$O(d)$, but multiplication $O(d^2)$.  
Another way is as $d+1$ evaluation over a given set of points (the domain).
Addition and multiplication become just $O(d)$ element-wise operations.  
The main issue is that switching between the 2 forms is $O(d^2)$ in general,
but particular cases can overcome this issue.  
It's worth mentioning that $d + 1$ is a requirement based on the degree of the
result, as long as there are enough evaluations, an arbitrary number of
operations can be performed.  
Most of the time, it is possible to predict a higher bound for the final
result and evaluate operands only once at the start.

## Evaluation and interpolation

These are the operations of switching between evaluations and coefficients.

### Evaluation

- The general algorithm for evaluating a degree $d$ polynomial over $n$ points
is $O(nd)$ multiplications and additions.
- Evaluating a degree 1 polynomial is the same, but due to $d=1$ it can be done
in $O(n).
- If your domain is $D = {0..n}$, then evaluation of a degree $1$ polynomial
can be done with only additions, as $p(i+1) = p(i) + c_1$ for the only degree 
1 coefficient $c_1$.
- With another special domain, most of the time of size $2^k$, interpolation is just
$O(n \log n)$ using FFT. But fields in general, and ed25519's in particular,
don't have such domain.
- ECFFT works similarly, but a suitable domain can be found for any field.
Performance is worse, around $O(n \log^2 n)$.

### Interpolation

- FFT and ECFFT work generally the same as for evaluation.
- For a general algorithm, it isn't as trivial as just evaluating the polynomial
in $n$ points like evaluation.
- But there is barycentric lagrange interpolation, which is $O(d^2)$ and works
on any domain. It requires a per domain, $O(d^2)$ precomputation.
- In some domains, precomputation can be just $O(n)$, but those domains are
generally the same you need for FFT.
- In the domain $D = 0..n$, part of the precomputation can be made $O(n)$, but
there is still a $O(d^2)$ component.

### Barycentric evaluation

The same values precomputed for interpolation can be used to evaluate a
polynomial, defined as the evaluations over a domain of size $n$ in
just $O(n)$.  
As such, if the only reason to interpolate a polynomial is to be able to
evaluate it on an arbitrary point, it can be done in just $O(n)$, and not
$O(n^2)$ as interpolation would require.

## Optimization approach

There are 2 main optimizations, described below. I didn't measure them
separately, but I remember going from $~1s$ to $90ms$ with the first one,
and then to the final $~15ms$.

### Optimizing divisor arithmetic

- All points are grouped in pairs. Each pair is interpolated into a divisor,
which being 2, is a line. Lines being degree 1 can be manipulated arbitrarily
with little cost. So far, this is the same as the original implementation.
- Now, the lines are converted into the representation described above, the
divisor is $a(x) - yb(x)$, with $a$ and $b$ degree 1.
- Both $a$ and $b$ are converted into evaluations and remain that way until
the end. Being degree 1, it can be done in just $O(n). They are all evaluated
in a domain big enough for the final result, which I set to 130, enough for
the typical 256 bits field.
- Then the algorithm remains mostly the same, as I wanted to get the exact
same result, the loop goes over the same partial divisors in the same order,
just a different representation. There is also a `SmallDivisor` with slightly
cheaper arithmetic when one of the operands is a line.
- At the end, the final divisor is interpolated and converted into the
bivariate representation expected by the tests.  
This is the only interpolation in the entire computation. It is $O(n^2)$,
but $n= 130$, which isn't that bad compared to an equivalent ECFFT. An
FFT would be better, but we can't have it for any field.

### Batching line computation

This is simpler in principle, but the changes to the original code end up being
more drastic than the point above. The problem is, that when divisor arithmetic
gets fast enough, interpolating the lines starts taking most of the time.
It can be attributed to a couple inverses, and to computing the $X$ and $Y$ of
each point, which unless you use affine coordinates, requires at least an
inverse to compute.
The solution to this is to batch all the coordinate translations and other
inverses so that they can be done with a single inversion.  
The new algorithm essentially runs the loop twice, with the first one
collecting all expensive operations, and computing the results in a single
batch. Then the original loop just takes those results instead of computing them.
An additional trait is required to be able to handle ed25519, as the library
makes it impossible to safely convert coordinates without expensive operations,
even with batching. The trait allows me to add a simple projective arithmetic
implementation to use instead. Curves that already expose a way to do batch
conversion can just use their own representation.

## Potential improvements

The are some possible improvements that could be done within the rules. And some
that wouldn't be compatible with them, but may still be valuable in practice.

### Within rules

- Intermediate interpolations: Currently, I evaluate lines on the final domain,
and interpolate only the final result. This results in $O(n^2)$ divisor arithmetic
in the loop, and a $O(n^2)$ interpolation of the final result.  
As an example, another option would be to evaluate lines over $n/2$ points,
this would make divisor arithmetic around $O({n/2}^2)$, but add 2 
intermediate $O({n/2}^2)$ interpolations/evaluations (you can go from
evaluations over a smaller domain to a bigger one without interpolation,
but still $O(n^2)$).  
My approach was just the simplest, but a different combination may end up with
a better runtime overall.
- Changing the loop: I kept the loop the same to be sure I was ending up with
the same divisor, but with the new representation, there are new options with
potentially better runtime to explore.  
One option would be to, instead of the current binary tree, to keep an
accumulator, and add points to it incrementally, increasing its degree by 1
each step. Complexity shouldn't be that different, but the constant factor
could be smaller.
- Barycentric precomputation: This is mostly already implemented, just not used
in order to comply with the rules without. It ended up being not that expensive,
but a few $ms$ can be saved.
- Batching a few more inverses: Each time 2 divisors are merged, the result has
to be divided by a line to remove the 2 partial sums that end up in the result.  
This requires a division, which is done by inverting the divisor and multiplying
by it.  
All evaluations of the divisor are batch inverted, but it could be possible to
go further and invert all the evaluations of all denominators at once. Which 
could cut the time of that part by half.

### Beyond rules

The main thing I can think here, is that the final interpolation may not be
necessary. ECIP already expects divisors to be of the form
$f(x,y) = a(x) - yb(x)$.  
Depending of what you want to do with a polynomial, you may be able to do it
with its evaluations.

- Evaluation in the domain: This trivial case is $O(1)$, as you already
have the evaluation.
- Single point evaluation outside domain: Can be $O(n)$ with barycentric
evaluation, not need for $O(n^2)$ interpolation.
- Addition of polynomial and scaling can be done $O(n)$, as it is just
element wise, multiplication requires having enough evaluations for the
degree of the result.
- If you want to commit to $a$ or $b$, Pedersen-based commitments allow
to commit directly to evaluations without interpolation, only a
precomputation equivalent in Complexity to interpolation.
