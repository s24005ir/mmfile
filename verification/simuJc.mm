Func void main()
{
    Real t0,t1;
    Matrix x0,TC,XC,UC;
    void diff_eqs(),link_eqs();


    t0=0.0;
    t1=10.0;
    x0=[PI-2.63108385 0]';//'
    
    print "simulating\n";

    {TC,XC,UC}=Ode45Auto(t0,t1,x0,diff_eqs,link_eqs,1e-10);

    mgplot(1,TC,XC(1,:));

    print [[TC][XC][UC]] >> "pendulum1.mat";

    Matrix t,x,d;

    read d << "dataJc4.mat";
    t=d(1,1:4540);
    x=d(4,461:5000).+PI;
    mgreplot(1,t,x);
}
Func void diff_eqs(DX,t,X,UY)
Real t;
Matrix X,DX,UY;
{
    Matrix xp,up,dxp;
    Real th,dth;
    Matrix A,B;
    Real m,l,M,f,J,c,a,g;

    m=0.03;
    l=0.155;
    M=0.713;
    f=7.57;
    J=2.4e-4;
    c=4.6e-5;
    a=0.624;
    g=9.8;

    xp=X;
    th=xp(1,1);
    dth=xp(2,1);

    A=[[0 1]
        [0,-c/(J+m*l*l)]];
    
    B=[[0]
        [-m*g*l*sin(th)/(J+m*l*l)]];

    dxp=[[dth]
            [-m*g*l*sin(th)/(J+m*l*l)-c*dth/(J+m*l*l)]];
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