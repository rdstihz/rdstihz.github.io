---
title: Lucas定理教程
mathjax: true
date: 2021-03-06 19:35:45
tags:
---


# Lucas定理
若$p$是质数，则对于任意的整数 $1 \le m \le n $, 有 
$$ C_n^m \equiv C_{n \ mod \ p}^{m \ mod \ p} * C_{n / p}^{m/p} （mod \ p)$$  
若如果把$n$写成$p$进制数$n_1n_2...n_k$,$m$写成$p$进制数$m_1m_2...m_k$,则  
$$ C_n^m \equiv C_{n_1}^{m_1} * C_{n_2}^{m_2}*...*C_{n_k}^{m_k} (mod \ p) $$  
当需要计算组合数并取模时,可以考虑使用Lucas定理.

## 证明
证明暂时留坑，以后再填。~~我还不会.~~

# 例题
## （1）Luogu P3807, 模板-卢卡斯定理
[题目链接](https://www.luogu.com.cn/problem/P3807)  

### 分析
这是卢卡斯定理模板题.
$$C_{n}^m = \frac{n!}{m!(n - m)!} =  \frac{ \prod_{i = m + 1}^{n+m} i }{(n - m)!}  $$  
考虑分别计算出分子a和分母b（一边累乘一边对p取模）,然后，将a乘以b模p的逆元，即可得到结果。  
但是，当$n-m$大于或等于$p$时，分母$b = (n-m)!$一定p的倍数，这时b不存在模p的逆元。因此不能通过上述方法算出结果。  
考虑Lucas定理，  
$$ C_n^m \equiv C_{n \ mod \ p}^{m \ mod \ p} * C_{n / p}^{m/p} （mod \ p)$$  
右边第一项中的$$n \ mod \ p$$和$$m \ mod \ p$$一定小于p，可以直接按上述方法计算，右边第二项$C_{n/p}^{m/p}$可以继续使用Lucas定理展开。 

### 代码
``` cpp
#include <iostream>

using namespace std;

typedef long long LL;

//快速幂
LL pow_mod(LL a, LL b, LL p) {
    LL res = 1;
    while (b) {
        if (b & 1)
            res = res * a % p;
        a = a * a % p;
        b >>= 1;
    }
    return res % p;
}

//组合数
LL C(LL n, LL m, LL p) {
    if (n < m) return 0;

    LL a = 1; //分子
    LL b = 1; //分母
    for (LL i = m + 1; i <= n; i++)
        a = a * i % p;
    for (LL i = 1; i <= n - m; i++)
        b = b * i % p;
    return a * pow_mod(b, p - 2, p) % p;
}

LL Lucas(LL n, LL m, LL p) { //用lucs定理计算C(n,m) % p 的值
    if (n - m < p)           // n - m < p 时可直接计算
        return C(n, m, p);
    else // n >= p 递归使用lucas定理
        return C(n % p, m % p, p) * Lucas(n / p, m / p, p) % p;
}

int main() {
    int T;
    cin >> T;
    while (T--) {
        LL n, m, p;
        cin >> n >> m >> p;
        cout << Lucas(n + m, m, p) << endl;
    }
    return 0;
}
```

## (2) Luogu P2480 古代猪文 
[题目链接](https://www.luogu.com.cn/problem/P2480) 
### 题意
给定整数$n,q$，求  
$$ q^{\sum_{d|n}C_n^d } \ mod \ 999911659 $$
的值.（999911659是质数）.
### 分析
根据费马小定理：$p$为质数且$a,p$互质时,有
$$ a^{p - 1} \equiv 1 (mod \ p)  $$  
推论
$$ a^x \equiv a^{x \ mod \ (p - 1)} (mod \ p) $$
所以，要计算$ q^{\sum_{d|n}C_n^d } \ mod \ 999911659 $，只需先算出$x = \sum_{d|n}C_n^d \ mod \ (999911658)$,再计算$q^x \ mod \ 999911659$.    
所以现在只需考虑如何去计算$x$.  
对于$n$的每个约数$d$,要计算$C_n^d \ mod \ 999911658$,考虑使用Lucas定理，但999911658不是质数，不满足Lucas定理的条件。  
将999911658分解质因数，得$999911658 = 2 \times 3 \times 4679 \times 35617$.将这4个质因数分别记为$a_1,a_2,a_3,a_4$.  
可用Lucas定理分别求出$\sum_{d|n}C_n^d$对质数$a_1,a_2,a_3,a_4$取模的结果，分别记为$b_1,b_2,b_3,b_4$  
最后，用中国剩余定理求解线性同余方程组：  
$$\begin{cases} x \ mod \ a_1 = b_1 \\ x \ mod \ a_2 = b_2 \\ x \ mod \ a_3 = b_3 \\ x \ mod \ a_4 = b_4
\end{cases}$$
求出其最小正整数解，就是$x = \sum_{d|n}C_n^d \ mod \ 999911658 $的值。 用快速幂算出 $$q^x \ mod \ 999911659$$ 的值即为最终结果.
另外，当q是质数999911659的倍数时，不满足费马小定理的条件，不能使用上述过程求解，需要特殊处理，直接输出0.

实现细节见代码。

``` cpp
#include <iostream>

using namespace std;

typedef long long LL;

//快速幂
LL pow_mod(LL a, LL b, LL p) {
    LL res = 1;
    while (b) {
        if (b & 1)
            res = res * a % p;
        a = a * a % p;
        b >>= 1;
    }
    return res % p;
}

//组合数
LL C(LL n, LL m, LL p) {
    if (n < m) return 0;

    LL a = 1; //分子
    LL b = 1; //分母
    for (LL i = m + 1; i <= n; i++)
        a = a * i % p;
    for (LL i = 1; i <= n - m; i++)
        b = b * i % p;
    return a * pow_mod(b, p - 2, p) % p;
}

LL Lucas(LL n, LL m, LL p) { //用lucs定理计算C(n,m) % p 的值
    if (n - m < p)           // n - m < p 时可直接计算
        return C(n, m, p);
    else // n >= p 递归使用lucas定理
        return C(n % p, m % p, p) * Lucas(n / p, m / p, p) % p;
}

const LL P = 999911659;

//a[1] - a[4]为 999911658的四个质因数。
LL a[] = {0, 2, 3, 4679, 35617};
LL b[5];
int main() {
    LL n, q;
    cin >> n >> q;

    if (q % P == 0) { // q是p的倍数，直接输出0
        cout << 0 << endl;
        return 0;
    }

    for (int i = 1; i <= 4; i++) {
        for (int d = 1; d * d <= n; d++) {
            if (n % d == 0) {
                b[i] = (b[i] + Lucas(n, d, a[i])) % a[i];
                if (d * d != n) {
                    b[i] = (b[i] + Lucas(n, n / d, a[i])) % a[i];
                }
            }
        }
    }
    //用中国剩余定理求解同余方程组
    LL x = 0;
    for (int i = 1; i <= 4; i++) {
        LL M = (P - 1) / a[i];
        x = (x + (b[i] * M % (P - 1)) * pow_mod(M, a[i] - 2, a[i])) % (P - 1);
    }
    cout << pow_mod(q, x, P) << endl;

    return 0;
}
```

