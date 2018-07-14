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

    {TC,XC,UC}=Ode45(t0,t1,x0,diff_eqs,link_eqs,h);

    mgplot(1,TC,XC(1:2,:));

    print [[TC][XC][UC]] >> "TXU45.mat";
}
Func void diff_eqs(DX,t,X,UY)
Real t;
Matrix X,DX,UY;
{
    Matrix xp,up,dxp;
    Real r,th,dr,dth;
    Matrix K,A;
    Real m,l,M,f,J,c,a,g;

    m=0.038;
    l=0.13;
    M=1.49;
    f=15.10;
    J=4.5e-4;
    c=2.1e-4;
    a=0.73;
    g=9.8;

    r=X(1,1);
    th=X(2,1);
    dr=X(3,1);
    dth=X(4,1);

    up=UY;

    K=[[M+m,m*l*cos(th)][m*l*cos(th),J+m*l^2]];
    A=K~*[[-f*dr+m*l*sin(th)*dth*dth+a*up(1,1)][m*g*l*sin(th)-c*dth]];
    dxp=[[dr]
         [dth]
         [A]];

    xp=X;
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