Func void main()
{
    Real t0,t1;
    Real h;
    Matrix x0,TC,XC,UC;
    void diff_eqs(),link_eqs();

    t0=0.0;
    t1=10.0;
    x0=[0 0]';//'
    
    print "simulating\n";

    {TC,XC,UC}=Ode45Auto(t0,t1,x0,diff_eqs,link_eqs,1e-8);

    mgplot(1,TC,XC(1,:));

    print [[TC][XC][UC]] >> "daisyasimu1.mat";

    read data1 << "zawa12.mat";
    T1=data1(1,:);
    X1=data1(3,:);

    mgreplot(1,T1,X1);
}
Func void diff_eqs(DX,t,X,UY)
Real t;
Matrix X,DX,UY;
{
    Matrix xp,up,dxp;
    Real r,dr;
    Matrix A,B;
    Real m,l,M,F,J,c,a,g;

    m=0.030;
    l=0.155;
    M=0.713189;
    F=7.57101;
    J=2.443e-4;
    c=4.818e-5;
    a=0.6237;
    g=9.8;

    r=X(1,1);
    dr=X(2,1);

    up=UY;
    xp=X;

    A=[[0 1]
       [0 ,-F/M]];
    
    B=[[0][a/M]];

    dxp=A*xp+B*up;
    DX=[dxp];
}
Func void link_eqs(UY,t,X)
Real t;
Matrix UY,X;
{
    Matrix xp,up;

    xp=X;

    UY=[12.0];
}