-- Configuring RELIC 0.7.0...

   DEBUG=[off|on] Build with debugging support.
   PROFL=[off|on] Build with profiling support.
   CHECK=[off|on] Build with error-checking support.
   VERBS=[off|on] Build with detailed error messages.
   OVERH=[off|on] Build with overhead estimation.
   DOCUM=[off|on] Build documentation.
   STRIP=[off|on] Build only selected algorithms.
   QUIET=[off|on] Build with printing disabled.
   COLOR=[off|on] Build with colored output.
   BIGED=[off|on] Build with big-endian support.
   SHLIB=[off|on] Build shared library.
   STLIB=[off|on] Build static library.
   STBIN=[off|on] Build static binaries.
   AMALG=[off|on] Build amalgamation for better performance.
   AUSAN=[off|on] Build with ASan and UBSan (gcc/clang only).

   TESTS=n        If n > 0, build automated tests and run them n times.
   BENCH=n        If n > 0, build automated benchmarks and run them n * n times.

   CORES=n        If n > 1, enable multithreading support.

   WITH=BN       Multiple precision arithmetic.
   WITH=DV       Temporary double-precision digit vectors.
   WITH=FP       Prime field arithmetic.
   WITH=FPX      Prime extension field arithmetic.
   WITH=FB       Binary field arithmetic.
   WITH=EP       Elliptic curves over prime fields.
   WITH=EPX      Elliptic curves over extensions of prime fields.
   WITH=EB       Elliptic curves over binary fields.
   WITH=ED       Elliptic Edwards curves over prime fields.
   WTTH=EC       Elliptic curve cryptography.
   WITH=PB       Pairings over binary elliptic curves.
   WITH=PP       Pairings over prime elliptic curves.
   WTTH=PC       Pairing-based cryptography.
   WITH=BC       Block ciphers (symmetric encryption).
   WITH=MD       Message digests (hash functions).
   WITH=CP       Cryptographic protocols.
   WITH=MPC      Multi-party computation primitives.
   WITH=ALL      All of the above.
   Note: the programmer is responsible for not using unselected modules.

   ARITH=easy     Easy-to-understand and portable, but slow backend.
   ARITH=fiat     Backend based on code generated from Fiat-Crypto.
   ARITH=gmp      Backend based on GNU Multiple Precision library.

   ARITH=gmp-sec  Same as above, but using constant-time code.

   ALLOC=AUTO     All memory is automatically allocated.
   ALLOC=DYNAMIC  All memory is allocated dynamically on demand.

   OPSYS=         Undefined/No specific operating system.
   OPSYS=LINUX    GNU/Linux operating system.
   OPSYS=FREEBSD  FreeBSD operating system.
   OPSYS=NETBSD   NetBSD operating system.
   OPSYS=MACOSX   Mac OS X operating system.
   OPSYS=WINDOWS  Windows operating system.
   OPSYS=DROID    Android operating system.
   OPSYS=DUINO    Arduino platform.

   MULTI=         No multithreading support.
   MULTI=OPENMP   Open Multi-Processing.
   MULTI=PTHREAD  POSIX threads.

   TIMER=         No timer.
   TIMER=HREAL    GNU/Linux realtime high-resolution timer.
   TIMER=HPROC    GNU/Linux per-process high-resolution timer.
   TIMER=HTHRD    GNU/Linux per-thread high-resolution timer.
   TIMER=ANSI     ANSI-compatible timer.
   TIMER=POSIX    POSIX-compatible timer.
   TIMER=CYCLE    Cycle-counting timer. (architecture-dependant)
   TIMER=PERF     GNU/Linux performance monitoring framework.

   LABEL=relic

   RLC_DEPTH=w    Width w in [2,8] of table for fixed exponentiation.
   RLC_WIDTH=w    Width w in [2,6] of table for exponentiation.

   ARCH=          No specific architecture (disable some features).
   ARCH=AVR       Atmel AVR ATMega128 8-bit architecture.
   ARCH=MSP       TI MSP430 16-bit architecture.
   ARCH=ARM       ARM 32-bit architecture.
   ARCH=A64       ARM 64-bit architecture.
   ARCH=X86       Intel x86-compatible 32-bit architecture.
   ARCH=X64       AMD x86_64-compatible 64-bit architecture.

   WSIZE=8        Build a 8-bit library.
   WSIZE=16       Build a 16-bit library.
   WSIZE=32       Build a 32-bit library.
   WSIZE=64       Build a 64-bit library.

   ALIGN=1        Do not align digit vectors.
   ALIGN=2        Align digit vectors into 16-bit boundaries.
   ALIGN=8        Align digit vectors into 64-bit boundaries.
   ALIGN=16       Align digit vectors into 128-bit boundaries.

   ** Options for the multiple precision module (default = 1024,DOUBLE,0):

      BN_PRECI=n        The base precision in bits. Let w be n in words.
      BN_MAGNI=DOUBLE   A multiple precision integer can store 2w words.
      BN_MAGNI=CARRY    A multiple precision integer can store w+1 words.
      BN_MAGNI=SINGLE   A multiple precision integer can store w words.
      BN_KARAT=n        The number of Karatsuba steps.

   ** Available multiple precision arithmetic methods (default = COMBA;COMBA;MONTY;SLIDE;BASIC;BASIC):

      Integer multiplication:
      BN_METHD=BASIC    Schoolbook multiplication.
      BN_METHD=COMBA    Comba multiplication.

      Integer squaring:
      BN_METHD=BASIC    Schoolbook squaring.
      BN_METHD=COMBA    Comba squaring.
      BN_METHD=MULTP    Reuse multiplication for squaring.

      Modular reduction:
      BN_METHD=BASIC    Division-based modular reduction.
      BN_METHD=BARRT    Barrett modular reduction.
      BN_METHD=MONTY    Montgomery modular reduction.
      BN_METHD=RADIX    Diminished radix modular reduction.

      Modular exponentiation:
      BN_METHD=BASIC    Binary modular exponentiation.
      BN_METHD=MONTY    Montgomery powering ladder.
      BN_METHD=SLIDE    Sliding window modular exponentiation.

      Greatest Common Divisor:
      BN_METHD=BASIC    Euclid's standard GCD algorithm.
      BN_METHD=BINAR    Binary GCD algorithm.
      BN_METHD=LEHME    Lehmer's fast GCD algorithm.

      Prime generation:
      BN_METHD=BASIC    Basic prime generation.
      BN_METHD=SAFEP    Safe prime generation.
      BN_METHD=STRON    Strong prime generation.

   ** Arithmetic precision of the prime field module (default = 256,0,off,off):

      FP_PRIME=n        The prime modulus size in bits.
      FP_KARAT=n        The number of Karatsuba levels.
      FP_PMERS=[off|on] Prefer Pseudo-Mersenne primes over random primes.
      FP_QNRES=[off|on] Use -1 as quadratic non-residue (make sure that p = 3 mod 8).
   ** Available prime field arithmetic methods (default = BASIC;COMBA;COMBA;MONTY;MONTY;JMPDS;SLIDE):
      Field addition:
      FP_METHD=BASIC    Schoolbook addition.
      FP_METHD=INTEG    Integrated modular addition.

      Field multiplication:
      FP_METHD=BASIC    Schoolbook multiplication.
      FP_METHD=INTEG    Integrated modular multiplication.
      FP_METHD=COMBA    Comba multiplication.

      Field squaring:
      FP_METHD=BASIC    Schoolbook squaring.
      FP_METHD=INTEG    Integrated modular squaring.
      FP_METHD=COMBA    Comba squaring.
      FP_METHD=MULTP    Reuse multiplication for squaring.

      Modular reduction:
      FP_METHD=BASIC    Division-based reduction.
      FP_METHD=QUICK    Fast reduction modulo special form prime (2^t - c, c > 0).
      FP_METHD=MONTY    Montgomery modular reduction.

      Field inversion:
      FP_METHD=BASIC    Inversion by Fermat's Little Theorem.
      FP_METHD=BINAR    Binary Inversion algorithm.
      FP_METHD=MONTY    Montgomery inversion.
      FP_METHD=EXGCD    Inversion by the Extended Euclidean algorithm.
      FP_METHD=DIVST    Constant-time inversion by division steps.
      FP_METHD=LOWER    Pass inversion to the lower level.

      Legendre symbol:
      FP_METHD=BASIC    Computation by Fermat's Little Theorem.
      FP_METHD=DIVST    Constant-time method by division steps.
      FP_METHD=JMPDS    Constant-time method by jump division steps.
      FP_METHD=LOWER    Pass call to the lower level.

      Field exponentiation:
      FP_METHD=BASIC    Binary exponentiation.
      FP_METHD=SLIDE    Sliding window exponentiation.
      FP_METHD=MONTY    Constant-time Montgomery powering ladder.

   ** Available bilinear pairing methods (default = BASIC;BASIC;BASIC):
      Quadratic extension arithmetic:
      FPX_METHD=BASIC    Basic quadratic extension field arithmetic.
      FPX_METHD=INTEG    Quadratic extension field arithmetic with embedded modular reduction.

      Cubic extension arithmetic:
      FPX_METHD=BASIC    Basic cubic extension field arithmetic.
      FPX_METHD=INTEG    Cubic extension field arithmetic with embedded modular reduction.

      Extension field arithmetic:
      FPX_METHD=BASIC    Basic extension field arithmetic.
      FPX_METHD=LAZYR    Lazy-reduced extension field arithmetic.

   ** Options for the binary elliptic curve module (default = 283,0,on,on,on):

      FB_POLYN=n        The irreducible polynomial size in bits.
      FB_KARAT=n        The number of Karatsuba levels.
      FB_TRINO=[off|on] Prefer trinomials.
      FB_SQRTF=[off|on] Prefer square-root friendly polynomials.
      FB_PRECO=[off|on] Precompute multiplication table for sqrt(z).
      RLC_WIDTH=w        Width w in [2,6] of window processing for exponentiation methods.

   ** Available binary field arithmetic methods (default = LODAH;QUICK;QUICK;BASIC;QUICK;QUICK;EXGCD;SLIDE;QUICK):

      Field multiplication:
      FB_METHD=BASIC    Right-to-left shift-and-add multiplication.
      FB_METHD=INTEG    Integrated modular multiplication.
      FB_METHD=LODAH    L�pez-Dahab comb multiplication with window of width 4.

      Field squaring:
      FB_METHD=BASIC    Bit manipulation squaring.
      FB_METHD=INTEG    Integrated modular squaring.
      FB_METHD=QUICK    Table-based squaring.

      Modular reduction:
      FB_METHD=BASIC    Shift-and-add modular reduction.
      FB_METHD=QUICK    Fast reduction modulo a trinomial or pentanomial.

      Field square root:
      FB_METHD=BASIC    Square root by repeated squaring.
      FB_METHD=QUICK    Fast square root extraction.

      Trace computation:
      FB_METHD=BASIC    Trace computation by repeated squaring.
      FB_METHD=QUICK    Fast trace computation.

      Quadratic equation solver:
      FB_METHD=BASIC    Solve a quadratic equation by half-trace computation.
      FB_METHD=QUICK    Fast solving with precomputed half-traces.

      Field inversion:
      FB_METHD=BASIC    Inversion by Fermat's Little Theorem.
      FB_METHD=BINAR    Binary Inversion algorithm.
      FB_METHD=ALMOS    Inversion by the Amost inverse algorithm.
      FB_METHD=EXGCD    Inversion by the Extended Euclidean algorithm.
      FB_METHD=ITOHT    Inversion by Itoh-Tsuji.
      FB_METHD=CTAIA    Constant-time almost inversion algorithm.
      FB_METHD=BRUCH    Hardware-friendly inversion by Brunner et al.
      FB_METHD=LOWER    Pass inversion to the lower level.

      Field exponentiation:
      FB_METHD=BASIC    Binary exponentiation.
      FB_METHD=SLIDE    Sliding window exponentiation.
      FB_METHD=MONTY    Constant-time Montgomery powering ladder.

      Iterated squaring/square-root:
      FB_METHD=BASIC    Iterated squaring/square-root by consecutive squaring/square-root.
      FB_METHD=QUICK    Iterated squaring/square-root by table-based method.

   ** Options for the prime elliptic curve module (default = all on):

      EP_PLAIN=[off|on] Support for ordinary curves without endomorphisms.
      EP_SUPER=[off|on] Support for supersingular curves.
      EP_ENDOM=[off|on] Support for ordinary curves with endomorphisms.
      EP_MIXED=[off|on] Use mixed coordinates.
      EP_CTMAP=[off|on] Use contant-time SSWU and isogeny map for hashing.

      EP_PRECO=[off|on] Build precomputation table for generator.
      RLC_DEPTH=w        Width w in [2,8] of precomputation table for fixed point methods.
   ** Available prime elliptic curve methods (default = PROJC;LWNAF;COMBS;INTER;SSWUM):

      Point representation:
      EP_METHD=BASIC    Affine coordinates.
      EP_METHD=PROJC    Homogeneous projective coordinates (complete formula).
      EP_METHD=JACOB    Jacobian projective coordinates.

      Variable-base scalar multiplication:
      EB_METHD=BASIC    Binary double-and-add method.
      EP_METHD=SLIDE    Sliding window method.
      EP_METHD=MONTY    Montgomery ladder method.
      EP_METHD=LWNAF    Left-to-right window NAF method.
      EP_METHD=LWREG    Left-to-right regular recoding method (GLV for curves with endomorphisms).

      Fixed-base scalar multiplication:
      EP_METHD=BASIC    Binary method for fixed point multiplication.
      EP_METHD=COMBS    Single-table Comb method for fixed point multiplication.
      EP_METHD=COMBD    Double-table Comb method for fixed point multiplication.
      EP_METHD=LWNAF    Left-to-right window NAF method (GLV for curves with endomorphisms).

      Variable-base simultaneous scalar multiplication:
      EP_METHD=BASIC    Multiplication-and-addition simultaneous multiplication.
      EP_METHD=TRICK    Shamir's trick for simultaneous multiplication.
      EP_METHD=INTER    Interleaving of window NAFs (GLV for Koblitz curves).
      EP_METHD=JOINT    Joint sparse form.

      Hash to point on the elliptic curve:
      EP_METHD=BASIC    Hash to x-coordinate and increment.
      EP_METHD=SSWUM    Simplified Shallue-van de Woestijne-Ulas method.
      EP_METHD=SWIFT    SwiftEC hashing method.

   ** Options for the binary elliptic curve module (default = on, w = 4):

      EB_PLAIN=[off|on] Support for ordinary curves without endomorphisms.
      EB_KBLTZ=[off|on] Support for Koblitz anomalous binary curves.
      EB_MIXED=[off|on] Use mixed coordinates.
      EB_PRECO=[off|on] Build precomputation table for generator.
   ** Available binary elliptic curve methods (default = PROJC;LWNAF;COMBS;INTER):

      Point representation:
      EB_METHD=BASIC    Affine coordinates.
      EB_METHD=PROJC    Projective coordinates (L�pez-Dahab for ordinary curves).

      Variable-base scalar multiplication:
      EB_METHD=BASIC    Binary double-and-add method.
      EB_METHD=LODAH    Lopez-Dahab constant-time point multiplication.
      EB_METHD=LWNAF    Left-to-right window (T)NAF method.
      EB_METHD=RWNAF    Right-to-left window (T)NAF method.
      EB_METHD=HALVE    Halving method.

      Fixed-base scalar multiplication:
      EB_METHD=BASIC    Binary method for fixed point multiplication.
      EB_METHD=COMBS    Single-table Comb method for fixed point multiplication.
      EB_METHD=COMBD    Double-table Comb method for fixed point multiplication.
      EB_METHD=LWNAF    Left-to-right window (T)NAF method.

      Variable-base simultaneous scalar multiplication:
      EB_METHD=BASIC    Multiplication-and-addition simultaneous multiplication.
      EB_METHD=TRICK    Shamir's trick for simultaneous multiplication.
      EB_METHD=INTER    Interleaving of window (T)NAFs.
      EB_METHD=JOINT    Joint sparse form.

   ** Options for the prime elliptic Edwards curve module (default = all on):
      ED_PRECO=[off|on] Build precomputation table for generator.
      RLC_DEPTH=w        Width w in [2,6] of precomputation table for fixed point methods.
   ** Available prime elliptic Edwards curve methods (default = PROJC;LWNAF;COMBS;INTER):
      ED_METHD=BASIC    Affine coordinates.
      EP_METHD=PROJC  	 Simple projective twisted Edwards coordinates.
      EP_METHD=EXTND 	 Extended projective twisted Edwards coordinates.

      *** variable-base multiplication method ***
      ED_METHD=BASIC    Binary method.
      ED_METHD=SLIDE    Sliding window method.
      ED_METHD=MONTY    Montgomery ladder method.
      ED_METHD=LWNAF    Left-to-right window NAF method.
      EP_METHD=LWREG    Left-to-right regular recoding method (GLV for curves with endomorphisms).

      *** fixed-base multiplication method ***
      ED_METHD=BASIC    Binary method for fixed point multiplication.
      ED_METHD=COMBS    Single-table Comb method for fixed point multiplication.
      ED_METHD=COMBD    Double-table Comb method for fixed point multiplication.
      ED_METHD=LWNAF    Left-to-right window NAF method.

      *** variable-base simultaneous multiplication method ***
      ED_METHD=BASIC    Multiplication-and-addition simultaneous multiplication.
      ED_METHD=TRICK    Shamir's trick for simultaneous multiplication.
      ED_METHD=INTER    Interleaving of window NAFs (GLV for Koblitz curves).
      ED_METHD=JOINT    Joint sparse form.

      Note: these methods must be given in order. Ex: ED_METHD="EXTND;LWNAF;COMBD;TRICK"

   ** Options for the binary elliptic curve module (default = on):

      EC_ENDOM=[off|on] Prefer (prime or binary) curves with endomorphisms.

   ** Available elliptic curve methods (default = PRIME):

      EC_METHD=PRIME    Use prime curves.
      EC_METHD=CHAR2    Use binary curves.
      EC_METHD=EDDIE    Use prime Edwards curves.

   ** Available bilinear pairing methods (default = BASIC;OATEP):

      Extension field arithmetic:
      PP_METHD=BASIC    Basic extension field arithmetic.
      PP_METHD=LAZYR    Lazy reduced extension field arithmetic.

      Pairing computation:
      PP_METHD=TATEP    Tate pairing.
      PP_METHD=WEILP    Weil pairing.
      PP_METHD=OATEP    Optimal ate pairing.

   ** Available hash functions (default = SH256):

      MD_METHD=SH224        SHA-224 hash function.
      MD_METHD=SH256        SHA-256 hash function.
      MD_METHD=SH384        SHA-384 hash function.
      MD_METHD=SH512        SHA-512 hash function.
      MD_METHD=B2S160       BLAKE2s-160 hash function.
      MD_METHD=B2S256       BLAKE2s-256 hash function.

   ** Options for the cryptographic protocols module (default = on, PKCS2):

      CP_CRT=[off|on] Support for faster CRT-based exponentiation in factoring-based cryptosystems.

      CP_RSAPD=BASIC    RSA with basic padding.
      CP_RSAPD=PKCS1    RSA with PKCS#1 v1.5 padding.
      CP_RSAPD=PKCS2    RSA with PKCS#1 v2.1 padding.

   RAND=HASHD     Use the HASH-DRBG generator. (recommended)
   RAND=RDRND     Use Intel RdRand instruction directly.
   RAND=UDEV      Use the operating system underlying generator.
   RAND=CALL      Override the generator with a callback.

   SEED=          Use a zero seed. (horribly insecure!)
   SEED=LIBC      Use rand()/random() functions. (insecure!)
   SEED=RDRND     Use Intel RdRand instruction directly.
   SEED=UDEV      Use non-blocking /dev/urandom. (recommended)
   SEED=WCGR      Use Windows' CryptGenRandom. (recommended)

