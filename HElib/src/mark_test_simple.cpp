
//stripping down the test_general to uderstand the lib better


#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include <NTL/lzz_pXFactoring.h>

#include <cassert>
#include <cstdio>

#define debugCompare(ea,sk,p,c)

/**************

1. c1.multiplyBy(c0)
2. c0 += random constant
3. c2 *= random constant
4. tmp = c1
5. ea.shift(tmp, random amount in [-nSlots/2, nSlots/2])
6. c2 += tmp
7. ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
8. c1.negate()
9. c3.multiplyBy(c2) 
10. c0 -= c3

**************/


void  TestIt(long R, long p, long r, long d, long c, long k, long w, long L, long m, const Vec<long>& gens, const Vec<long>& ords) {
  char buffer[32];

  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);

  FHEcontext context(m, p, r, gens1, ords1);
  buildModChain(context, L, c);

  context.zMStar.printout();
  
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key


  ZZX G;

  if(d == 0) {
    G = context.alMod.getFactorsOverZZ()[0];
  } else {
    G = makeIrredPoly(p, d);
  }

  addSome1DMatrices(secretKey); // compute key-switching matrices that we need

  EncryptedArray ea(context, G);

  cerr << "ea size: " << ea.size() << endl;

  //create empty plaintext
  NewPlaintextArray p0(ea);

  //create random vector (binary?)
  random(ea, p0);

  cerr << "plaintext (initial): " << p0 << endl;

  //create encrypted cypher text
  Ctxt c0(publicKey);
  ea.encrypt(c0, publicKey, p0);

  //cerr << "cyphertext: " << c0.parts() << endl;

  NewPlaintextArray pp0(ea);
  ea.decrypt(c0, secretKey, pp0);

  cerr << "plaintext (decrypt): " << pp0 << endl;

/*
  resetAllTimers();

  FHE_NTIMER_START(Circuit);

  for (long i = 0; i < R; i++) {

    cerr << "*** round " << i << "..." << endl;

     long shamt = RandomBnd(2*(nslots/2) + 1) - (nslots/2);
                  // random number in [-nslots/2..nslots/2]
     long rotamt = RandomBnd(2*nslots - 1) - (nslots - 1);
                  // random number in [-(nslots-1)..nslots-1]

     // two random constants
     NewPlaintextArray const1(ea);
     NewPlaintextArray const2(ea);
     random(ea, const1);
     random(ea, const2);

     ZZX const1_poly, const2_poly;
     ea.encode(const1_poly, const1);
     ea.encode(const2_poly, const2);

     mul(ea, p1, p0);     // c1.multiplyBy(c0)
     c1.multiplyBy(c0);              CheckCtxt(c1, "c1*=c0");
     debugCompare(ea,secretKey,p1,c1);

     add(ea, p0, const1); // c0 += random constant
     c0.addConstant(const1_poly);    CheckCtxt(c0, "c0+=k1");
     debugCompare(ea,secretKey,p0,c0);

     mul(ea, p2, const2); // c2 *= random constant
     c2.multByConstant(const2_poly); CheckCtxt(c2, "c2*=k2");
     debugCompare(ea,secretKey,p2,c2);

     NewPlaintextArray tmp_p(p1); // tmp = c1
     Ctxt tmp(c1);
     sprintf(buffer, "c2>>=%d", (int)shamt);
     shift(ea, tmp_p, shamt); // ea.shift(tmp, random amount in [-nSlots/2,nSlots/2])
     ea.shift(tmp, shamt);           CheckCtxt(tmp, buffer);
     debugCompare(ea,secretKey,tmp_p,tmp);

     add(ea, p2, tmp_p);  // c2 += tmp
     c2 += tmp;                      CheckCtxt(c2, "c2+=tmp");
     debugCompare(ea,secretKey,p2,c2);

     sprintf(buffer, "c2>>>=%d", (int)rotamt);
     rotate(ea, p2, rotamt); // ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
     ea.rotate(c2, rotamt);          CheckCtxt(c2, buffer);
     debugCompare(ea,secretKey,p2,c2);

     ::negate(ea, p1); // c1.negate()
     c1.negate();                    CheckCtxt(c1, "c1=-c1");
     debugCompare(ea,secretKey,p1,c1);

     mul(ea, p3, p2); // c3.multiplyBy(c2) 
     c3.multiplyBy(c2);              CheckCtxt(c3, "c3*=c2");
     debugCompare(ea,secretKey,p3,c3);

     sub(ea, p0, p3); // c0 -= c3
     c0 -= c3;                       CheckCtxt(c0, "c0=-c3");
     debugCompare(ea,secretKey,p0,c0);

  }

  c0.cleanUp();
  c1.cleanUp();
  c2.cleanUp();
  c3.cleanUp();

  FHE_NTIMER_STOP(Circuit);
   
  cerr << endl;
  printAllTimers();
  cerr << endl;
   
  resetAllTimers();
  FHE_NTIMER_START(Check);
   
  NewPlaintextArray pp0(ea);
  NewPlaintextArray pp1(ea);
  NewPlaintextArray pp2(ea);
  NewPlaintextArray pp3(ea);
   
  ea.decrypt(c0, secretKey, pp0);
  ea.decrypt(c1, secretKey, pp1);
  ea.decrypt(c2, secretKey, pp2);
  ea.decrypt(c3, secretKey, pp3);
   
  if (!equals(ea, pp0, p0))  cerr << "oops 0\n";
  if (!equals(ea, pp1, p1))  cerr << "oops 1\n";
  if (!equals(ea, pp2, p2))  cerr << "oops 2\n";
  if (!equals(ea, pp3, p3))  cerr << "oops 3\n";
   
  FHE_NTIMER_STOP(Check);
   
  cerr << endl;
  printAllTimers();
  cerr << endl;
*/
}




/* A general test program that uses a mix of operations over four ciphertexts.
 * Usage: Test_General_x [ name=value ]...
 *   R       number of rounds  [ default=1 ]
 *   p       plaintext base  [ default=2 ]
 *   r       lifting  [ default=1 ]
 *   d       degree of the field extension  [ default=1 ]
 *              d == 0 => factors[0] defines extension
 *   c       number of columns in the key-switching matrices  [ default=2 ]
 *   k       security parameter  [ default=80 ]
 *   L       # of levels in the modulus chain  [ default=heuristic ]
 *   s       minimum number of slots  [ default=0 ]
 *   repeat  number of times to repeat the test  [ default=1 ]
 *   m       use specified value as modulus
 *   mvec    use product of the integers as  modulus
 *              e.g., mvec='[5 3 187]' (this overwrite the m argument)
 *   gens    use specified vector of generators
 *              e.g., gens='[562 1871 751]'
 *   ords    use specified vector of orders
 *              e.g., ords='[4 2 -4]', negative means 'bad'
 */
int main(int argc, char **argv) {
  //setup
  setTimersOn();

  long R = 1;
  long p = 2;
  long r = 1;
  long d = 1;
  long c = 2;
  long k = 80;
  long L = 0;
  long s = 0;
  long repeat = 1;
  long chosen_m = 0;
  Vec<long> mvec, gens, ords;
  long seed = 0;
  long nt = 1;

  SetSeed(ZZ(seed));
  SetNumThreads(nt);
  
  if(L == 0) { // determine L based on R,r
    L = 3*R + 3;
    if(p>2 || r>1) { // add some more primes for each round
      long addPerRound = 2*ceil(log((double)p)*r*3)/(log(2.0)*FHE_p2Size) +1;
      L += R * addPerRound;
    }
  }

  long w = 64; // Hamming weight of secret key, number of ones
  //  long L = z*R; // number of levels

  if(mvec.length() > 0) {
    chosen_m = computeProd(mvec);
  }
  long m = FindM(k, L, c, p, d, s, chosen_m, true);

  TestIt(R, p, r, d, c, k, w, L, m, gens, ords);
}

