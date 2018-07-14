Func void calc()
{
    Real m,l,lambda,tau,a,g,K,T;
    Real J,c,M,F;

    m=0.03;
    l=0.155;
    lambda=0.02204813;
    tau=0.9175;
    a=0.6237;
    g=9.8;
    K=0.08238;
    T=0.0942;

    J=(m*g*l*tau*tau)/(4*PI*PI+lambda*lambda)-m*l*l;
    c=2*lambda*(J+m*l*l)/tau;

    F=a/K;
    M=T*F;

    print J;
    print c;
    print F;
    print M;
}
Func void calc2()
{
    Real m,l,lam,tau,a,g,K,T,wn,te,kc;
    Real M,F,A,B,C,D;

    m=0.03;
    l=0.155;

    lam=1.1076;
    tau=0.1675;

    a=0.6237;
    g=9.8;
    K=0.08238;
    T=0.0942;
    kc=1500;

    A=(4*PI*PI+lam*lam)^0.5;
    te=lam/A;
 
    B=(1-te*te)^0.5;
    C=tau*B;
    wn=2*PI/C;

    D=wn*wn;
    M=kc*a/D;
    F=2*te*wn*M;

    print te;
    print wn;
    print M;
    print F;
}