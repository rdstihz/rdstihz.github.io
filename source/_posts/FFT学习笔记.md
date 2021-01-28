---
title: FFT学习笔记
mathjax: true
date: 2021-01-28 13:04:48
tags:
---

# 多项式

### 系数表达
$$A(x) = \sum_{i=0}^{n-1}a_ix^i$$  
### 点值表达
$$ \{(x_0,y_0), (x_1,y_1),(x_2,y_2),...,(x_n-1,y_n-1) \} $$

# 算法概括
FFT可以在$O(nlogn)$的时间复杂度内计算两个一元多项式的乘积。  
主要过程为：
1. 计算原有两个多项式在单位复数根下的点值。（DFT）
2. 将两个多项式的点值相乘。得到乘积的点值。
3. 将点值还原为系数表达。（IDFT）

# 单位复数根

### $n$次单位复数根
$$w_n^k = e^{i\frac{2\pi k}{n}}= cos(\frac{2\pi k}{n})+isin(\frac{2\pi k}{n})\quad k=0,1,2...,n-1$$

### 单位根的性质

#### 消去定理
$$w_{2n}^{2k} = w_n^k $$
#### 折半定理
$n$为偶数时,$n$个$n$次单位根的集合就是$n/2$个$n/2$次单位根的集合。

#### 其它
$$w_n^{k+n/2} = -w_n^k  $$
$$w_n^{k+n} = w_n^{k}  $$


# 离散傅利叶变换（DFT）
设多项式$$A(x) = a_0+a_1x+a_2x^2+a_3x^3+...+a_{n-1}x^{n-1}$$
他的所有的偶数次项系数组成新多项式  
$$A_0(x) = a_0 + a_2x + a_4x^2 + ... $$
$$A_1(x) = a_1 + a_3x + a_4x^2 + ... $$
则显然有$$A(x) = A_0(x^2)+xA_1(x^2) $$
$$A(w_n^k) = A_0(w_n^{2k})+w_n^kA_1(w_n^{2k})= A_0(w_{n/2}^{k})+w_n^kA_1(w_{n/2}^{k}) $$ 

$$A(w_n^{k+n/2}) = A_0(w_n^{2k+n})+w_n^kA_1(w_n^{2k+n})= A_0(w_{n/2}^{k})-w_n^kA_1(w_{n/2}^{k}) $$

利用分治思想，将$A(x)$的奇数项和偶数项分别拆出一个新多多项式，
递归计算两个新多项式的点值后，可在$O(n)$时间内求出$A(x)$的点值。
总时间复杂度为$O(nlogn)$

## 伪代码

### 递归实现
``` 
DFT(A,n) //A为系数数组，n为项数，且n为2的整数次幂（可在）
if n = 1
    return

m = n/2
for i = 0 to m
    A0[i] = A[2*i]
    A1[i] = A[2*i+1]

DFT(A0,n)
DFT(A1,n)

w = 1
wn = cos(2pi/n) + i*sin(2pi/n)
for i = 0 to m-1
    A[i] = A0[i] + w*A1[i]
    A[i+m] = A0[i] - w*A1[i]
    w = w*wn
return
```

### 迭代实现  
将系数数组排成$log_2n$次奇偶分离操作后的顺序  
这将顺序每个数的二进制秋翻转后即为0~n-1从小到大的排列
每个数翻转后的数可以递推求出。
C++代码:
``` cpp

struct Complex{}//实现复数类，此外省略

//1. rev 的计算 （放在main函数中）
for(int i = 0;i<n;i++)
    rev[i] = (rev[i>>1]>>1) | ((i&1) >> l-1)

void DFT(Complex *A,int n ){
    for(int i = 0;i<n;i++)
        if(i<rev[i]) swap(A[i],A[r[i]]);
    
    for(int len = 2;len <= n;len <<= 1){ 
        //从短到长枚举要外理的长度，相当于从递归的底层到顶层
        
        Complex wn(cos(2*pi/len), sin(2*pi/len) );
        for(int j = 0;j<n;j+=len){//枚举要外理的部分的长度
            int m = len>>1;
            Complex w(1,0);
            for(int i = 0;i<m;i++){
                Complex u = a[j+i];
                Complex v = a[j+i+m];
                a[j+i] = u + w*v;
                a[j+i+m] = u-w*v;//蝴蝶操作
                w = w*wn; 
            }
        }
    }
}
```

# 离散傅利叶逆变换（IDFT）

将单位根的枚举顺序倒过来，进行一次类似DFT操作，将最后得到的每个项点值除以$n$即得到每一项系数。
代码可以与DFT合并 

``` cpp
void FFT(Complex *A,int n ,int type){
    //type = 1表示DFT
    //type = -1 表示IDFT

    for(int i = 0;i<n;i++)
        if(i<rev[i]) swap(A[i],A[r[i]]);
    
    for(int len = 2;len <= n;len <<= 1){ 
        //从短到长枚举要外理的长度，相当于从递归的底层到顶层
        
        Complex wn(cos(2*pi/len), type*sin(2*pi/len) );
        for(int j = 0;j<n;j+=len){//枚举要外理的部分的长度
            int m = len>>1;
            Complex w(1,0);
            for(int i = 0;i<m;i++){
                Complex u = A[j+i];
                Complex v = A[j+i+m];
                A[j+i] = u + w*v;
                A[j+i+m] = u-w*v;//蝴蝶操作
                w = w*wn; 
            }
        }
    }
}

```


# 完整代码
``` cpp
#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

const int maxn = 10000000+100;
const double pi = acos(-1);

//实现复数类
struct Complex{
    double r,v;//实部和虚部
    Complex(){
        r = v = 0;
    }
    Complex(double a,double b){
        r = a; v = b;
    }
    
    Complex operator+(const Complex &p) const{
        return Complex(r+p.r,v+p.v);
    }
    Complex operator-(const Complex &p) const{
        return Complex(r-p.r,v-p.v);
    }
    Complex operator*(const Complex &p) const{
        return Complex(r*p.r-v*p.v ,r*p.v + v*p.r ); 
    }
    
};

Complex a[maxn],b[maxn];

int rev[maxn];

//DFT和IDFT
void fft(Complex * A,int n,int type){
    
    for(int i = 0;i<n;i++)
        if(rev[i] < i ) swap(A[i],A[rev[i]]);
    
    for(int len = 2;len <= n;len <<= 1){//the length of doing 
        Complex wn(cos(2*pi/len), type*sin(2*pi/len) );
        for(int j = 0;j<n;j+=len){
            int m = len>>1;
            
            Complex w(1,0);
            
            for(int i = 0;i<m;i++){
                Complex u = A[j+i];
                Complex v = A[j+i+m];
                A[j+i] = u + w*v;
                A[j+i+m] = u-w*v;
                w = w*wn; 
            }
            
        }
    }
}

//读入优化
inline int read(){
    int x = 0;
    char c =getchar();
    while(c<'0'||c>'9') c = getchar();
    while(c>='0' && c<='9' ){
        x = x*10 + c-'0';
        c = getchar();
    }
    return x;
    
}

int main(){
    int n,m;
    n =  read();
    m = read();

    for(int i = 0;i<=n;i++){
        a[i].r = read();
    }
    
    for(int i = 0;i<=m;i++){
        b[i].r = read();
    }
    
    
    //将项数补齐到2的整数次幂
    
    int N = 1;
    int l = 0;
    while(N<=n+m) N <<= 1,l++;
    
    //求rev
    for(int i =0;i<N;i++)
        rev[i] = (rev[i>>1]>>1) | ((i&1) << (l-1) );
    
    fft(a,N,1);
    fft(b,N,1);
    
    for(int i = 0;i<N;i++)
        a[i] = a[i]*b[i];
    fft(a,N,-1);
    
    for(int i =0;i<=n+m;i++)
        printf("%d ",(int)(a[i].r/N+0.5));//+0.5实现四舍五入
    printf("\n");
    return 0;
}

```



