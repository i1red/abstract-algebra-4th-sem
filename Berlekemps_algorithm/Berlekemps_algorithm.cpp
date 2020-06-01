#include <iostream>
#include <string>
#include <stdio.h>

using namespace std;

int main()
{

  string s = { 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0 };
  int L, N, m, d;
  int n = s.length();
  int f[11];
  int b[11];
  int t[11];

  //Initialization
  b[0] = f[0] = 1;
  N = L = 0;
  m = -1;

  //Algorithm core
  while (N < n)
  {
    d = s[N];
    for (int i = 1; i <= L; i++)
      d ^= f[i] & s[N - i];            //(d+=c[i]*s[N-i] mod 2)
    if (d == 1)
    {
      memcpy(t, f, n);    //T(D)<-C(D)
      for (int i = 0; (i + N - m) < n; i++)
        f[i + N - m] ^= b[i];
      if (L <= (N >> 1))
      {
        L = N + 1 - L;
        m = N;
        memcpy(b, t, n);    //B(D)<-T(D)
      }
    }
    N++;
  }
  cout << L << endl;
  system("pause");
  return 0;
}
