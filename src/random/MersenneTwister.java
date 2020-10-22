/*
  A C-program for MT19937, with initialization improved 2002/1/26.
  Coded by Takuji Nishimura and Makoto Matsumoto.

  Before using, initialize the state by using initGenRand(seed)
  or initByArray(init_key, key_length).

  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
  All rights reserved.
  Copyright (C) 2005, Mutsuo Saito,
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

   1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   3. The names of its contributors may not be used to endorse or promote
   products derived from this software without specific prior written
   permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


  Any feedback is very welcome.
  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
 */

package random;

public class MersenneTwister {
    private final static int N = 624;
    private final static int M = 397;
    private final static long MATRIX_A = 0x9908b0dfL;
    private final static long UPPER_MASK = 0x80000000L;
    private final static long LOWER_MASK = 0x7fffffffL;

    private final static long[] mt = new long[N];
    private static int mti = N + 1;

    // Initializes mt[N] with a seed
    public static void initGenRand(long s)
    {
        mt[0] = s & 0xffffffffL;
        for (mti = 1; mti < N; mti++) {
            mt[mti] = 1812433253L * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti;
            mt[mti] &= 0xffffffffL;
        }
    }

    public static void initByArray(long[] initKey, int keyLength)
    {
        initGenRand(19650218L);
        int i = 1;
        int j = 0;
        int k = Math.max(N, keyLength);

        for (; k != 0; k--) {
            mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525L)) + initKey[j] + j;
            mt[i] &= 0xffffffffL;
            i++; j++;
            if (i >= N) {
                mt[0] = mt[N-1];
                i = 1;
            }
            if (j >= keyLength) {
                j = 0;
            }
        }

        for (k = N - 1; k != 0; k--) {
            mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941L)) - i;
            mt[i] &= 0xffffffffL;
            i++;
            if (i >= N) {
                mt[0] = mt[N-1];
                i = 1;
            }
        }
        mt[0] = 0x80000000L;
    }

    // generates a random number on [0,0xffffffff]-interval
    public static long genRandInt32()
    {
        int kk;
        long y;
        long[] mag01 = {0x0L, MATRIX_A};

        if (mti >= N) {
            if (mti == N + 1) {
                initGenRand(5489L);
            }

            for (kk = 0; kk < N - M; kk++) {
                y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
                mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[(int)(y & 0x1L)];
            }

            for (; kk < N - 1; kk++) {
                y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
                mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[(int)(y & 0x1L)];
            }

            y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[(int)(y & 0x1L)];

            mti = 0;
        }

        y = mt[mti++];

        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680L;
        y ^= (y << 15) & 0xefc60000L;
        y ^= (y >> 18);

        return y;
    }

    // generates a random number on [0,0x7fffffff]-interval
    public static long genRandInt31()
    {
        return genRandInt32() >> 1;
    }

    // generates a random number on [0,1]-real-interval
    public static double genRandReal1()
    {
        return genRandInt32() * (1.0 / 4294967295.0);
    }

    // generates a random number on [0,1)-real-interval
    public static double genRandReal2()
    {
        return genRandInt32() * (1.0 / 4294967296.0);
    }

    // generates a random number on (0,1)-real-interval
    public static double genRandReal3()
    {
        return (((double)genRandInt32()) + 0.5) * (1.0 / 4294967296.0);
    }

    // generates a random number on [0,1) with 53-bit resolution
    public static double genRandRes53()
    {
        long a = genRandInt31() >> 5;
        long b = genRandInt31() >> 6;
        return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
    }

    public static void main(String[] args)
    {
        long[] init = {0x123L, 0x234L, 0x345L, 0x456L};
        int length = 4;
        initByArray(init, length);

        System.out.println("1000 outputs of genRandInt32()");
        for (int i = 0; i < 1000; i++) {
            System.out.println(genRandInt31());
            if (i % 5 == 4) {
                System.out.println();
            }
        }

        System.out.println("1000 outputs of genRandReal2()");
        double sum = 0.0;
        for (int i = 0; i < 1000; i++) {
            sum += genRandReal2();
        }
        System.out.println(sum / 1000);
    }
}
