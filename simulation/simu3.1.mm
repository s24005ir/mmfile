Func void main()
{
    Real t0,t1;
    Real h;
    Matrix x0,TC,XC,UC;
    void diff_eqs(),link_eqs();

    t0=0.0;
    t1=10.0;
    x0=[0 10.0/180*PI 0 0]';
    h=0.01;
    
    print "simulating\n";

    {TC,XC,UC}=Ode(t0,t1,x0,diff_eqs,link_eqs,h);

    mgplot(1,TC,XC(1:2,:));

    print [[TC][XC][UC]] >> "TXUenchokuSenkei.mat";
}
Func void diff_eqs(DX,t,X,UY)
Real t;
Matrix X,DX,UY;
{
    Matrix xp,up,dxp;
    Matrix K,A,A21,A22,B,B2;
    Real m,l,M,f,J,c,a,c1,c2,g;

    m=0.038;
    l=0.13;
    M=1.49;
    f=15.10;
    J=4.5e-4;
    c=2.1e-4;
    a=0.73;
    c1=1.0;
    c2=1.0;
    g=9.8;

    K=[[M+m, -m*l][-m*l, J+m*l^2]];
    A21=K~*[[0, 0][0, -m*g*l]];
    A22=K~*[[-f, 0][0, -c]];
    A=[[Z(2), I(2)][A21, A22]];
    B2=K~*[a, 0]';
    B=[[0] [0] [B2]];


    xp=X;
    up=UY;
    dxp=A*xp+B*up;
    DX=[dxp];
}
Func void link_eqs(UY,t,X)
Real t;
Matrix UY,X;
{
    Matrix xp,up;

    xp=X;

    UY=[0];
}