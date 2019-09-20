# Pollard-Kangaroo, C99

Pollard, kangaroo method, solving discrete logarithm problem (DLP) using pseudorandom walks, C99.
Its the kangaroo method (of [11, p. 13]) using distinguished points.
Runtime expected of 2w<sup>1/2</sup> group operations.

[11] P. C. van Oorschot and M. J. Wiener, Parallel collision search with cryptanalytic applications, J. Cryptology, #12 (1999)

# Feature:

 - singlecore
 - bitcoin-core/secp256k1 optimized library
 - Affine points addition: A+A->A

Expected in the future
 - multicore

# Compilation

Its C99: 
```
gcc -O2 -I secp256k1/src/ -I secp256k1/ pollard-kangaroo.c -lgmp
```
The secp256k1 bitcoin library needs to be precompiled : 
```
./autogen.sh
./configure --enable-endomorphism --enable-ecmult-static-precomputation  --with-bignum=gmp --with-scalar=64bit --with-field=64bit --with-asm=x86_64 
make
```

How compile under Windows - can use Cygwin.

For Cygwin, libraries must be installed:
```
cygport
bash, make, perl, gcc-core, autoconf, automake, binutils, libtool, pkg-config
libgmp10, devgmp
openssl-devel

example versions:
  autoconf-13-1
  automake-10-1
  binutils-2.29-1
  cygport-0.31.1-1
  gcc-core-7.3.0-2
  libtool-2.4.6-6
  make-4.2.1-2
  libgmp10-6.1.2-1
  pkg-config-0.29.1-1
 ```

for Windows+Cygwin recommended "-no-undefined" of secp256k1 :
```
./autogen.sh
./configure --enable-endomorphism --enable-ecmult-static-precomputation  --with-bignum=gmp --with-scalar=64bit --with-field=64bit --with-asm=x86_64 
make LDFLAGS=" -no-undefined "
```

# Benchmark libs
Algo: 1 Tame + 1 Wild with distinguished points,  expected of 2w<sup>1/2</sup> group operations

1core i5-2540, win7x64, algo A+A->A
```
0.219Mk/s Mk/s
```

1core i7-6820, win7x64, algo A+A->A
```
0.291Mk/s
```

# Discussion
more info

https://bitcointalk.org/index.php?topic=5173445.msg52473992#msg52473992

https://bitcointalk.org/index.php?topic=5166284.msg52318676#msg52318676

# How pollard-kangaroo works, the Tame and Wild kangaroos, is a simple explanation.

Suppose there is pubkeyX, unknow privkeyX, but privkeyX is in range w=[L..U]. 
The keys have a property - if we increase pubkey by S, then its privkey will also increase by S. 
We start step-by-step to increase pubkeyX by S(i), keeping sumX(S(i)). This is a Wild kangaroo. 
We select a random privkeyT from range [L..U], compute pubkeyT. 
We start step-by-step to increment pubkeyT by S(i) while maintaining sumT(S(i)). This is a Tame kangaroo. 
The size of the jump S(i) is determined by the x coordinate of the current point, so if a Wild or Tame kangaroo lands on one point, their paths will merge. 
(we are concerned with pseudo random walks whose next step is determined by the current position) 
Thanks to the Birthday Paradox (Kruskal's card trick), their paths will one day meet. 
Knowing each traveled path (sumX and sumT), privkeyX is calculated. 
The number of jumps is approximately 2w<sup>1/2</sup> group operations, which is much less than a full search w. 

# Articles

base
https://andrea.corbellini.name/2015/05/17/elliptic-curve-cryptography-a-gentle-introduction/

hand translate to ru/ua (recommend instead translate.google)
https://habr.com/ru/post/335906/


- [1] Best of the best, all in 1, epic,  2012

Chapter 14. Factoring and Discrete Logarithms using Pseudorandom Walks 
https://www.math.auckland.ac.nz/~sgal018/crypto-book/ch14.pdf

- [2] with reference to old

J. M. Pollard, “Kangaroos, monopoly and discrete logarithms,” Journal of Cryptology, #13 (2000) 
https://web.northeastern.edu/seigen/11Magic/KruskalsCount/PollardKangarooMonopoly.pdf

(good dir web.northeastern.edu/seigen/11Magic/KruskalsCount/)

- [3] About parallelism problems

P. C. van Oorschot and M. J. Wiener, Parallel collision search with cryptanalytic applications, J. Cryptology, #12 (1999) 
https://people.scs.carleton.ca/~paulv/papers/JoC97.pdf

Pollard, kangaroo method, solving discrete logarithm problem (DLP) using pseudorandom walks, C99
Its the kangaroo method (of [11, p. 13]) using distinguished points.
Runtime expected of 2w^(1/2) group operations.

[11] P. C. van Oorschot and M. J. Wiener, Parallel collision search with cryptanalytic applications, J. Cryptology, #12 (1999)
