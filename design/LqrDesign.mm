Matrix A,B,F;
Matrix Q,R,P;
Func void main()
{
    Real t0,t1,r0,th0,tol;
    Matrix x0,T,X,U;
    Matrix diff_eqs(),link_eqs();


    t0=0.0;
    t1=10.0;
    r0=-0.1;
    th0=20;
    x0=[r0 th0/180*PI 0 0]'; //'
    tol = 1.0E-8;
    
    print "simulating\n";

    read A << "A.mat";
    read B << "B.mat";

    //LQ最適制御
    Q=diag(1000,10,1,1);//重み行列
    R=[1];
    {F,P}=lqr(A,B,Q,R);

    {T,X,U}=Ode45Auto(t0,t1,x0,diff_eqs,link_eqs,tol);

    mgplot(1,T,X(1,:),{"r1"});
    mgplot(2,T,X(2,:),{"th1"});
    mgplot(3,T,X(3,:),{"dr1"});
    mgplot(4,T,X(4,:),{"dth1"});
    mgplot(5,T,U(1,:),{"u1"});

    print [[T][X][U]] >> "ctrlPendulum1.mat";

    //LQ最適制御
    Q=diag(100,100,1,1);//重み行列
    {F,P}=lqr(A,B,Q,R);

    {T,X,U}=Ode45Auto(t0,t1,x0,diff_eqs,link_eqs,tol);

    mgreplot(1,T,X(1,:),{"r2"});
    mgreplot(2,T,X(2,:),{"th2"});
    mgreplot(3,T,X(3,:),{"dr2"});
    mgreplot(4,T,X(4,:),{"dth2"});
    mgreplot(5,T,U(1,:),{"u2"});

    print [[T][X][U]] >> "ctrlPendulum2.mat";

    //LQ最適制御
    Q=diag(100,100,30,1);//重み行列
    {F,P}=lqr(A,B,Q,R);

    {T,X,U}=Ode45Auto(t0,t1,x0,diff_eqs,link_eqs,tol);

    mgreplot(1,T,X(1,:),{"r3"});
    mgreplot(2,T,X(2,:),{"th3"});
    mgreplot(3,T,X(3,:),{"dr3"});
    mgreplot(4,T,X(4,:),{"dth3"});
    mgreplot(5,T,U(1,:),{"u3"});

    print [[T][X][U]] >> "ctrlPendulum3.mat";

    {T,X,U}=Ode45Auto(t0,t1,x0,diff_eqs,link_eqs,tol);

    //LQ最適制御
    Q=diag(1,1,1,100);//重み行列
    {F,P}=lqr(A,B,Q,R);

    mgreplot(1,T,X(1,:),{"r4"});
    mgreplot(2,T,X(2,:),{"th4"});
    mgreplot(3,T,X(3,:),{"dr4"});
    mgreplot(4,T,X(4,:),{"dth4"});
    mgreplot(5,T,U(1,:),{"u4"});

    print [[T][X][U]] >> "ctrlPendulum4.mat";
}

Func Matrix diff_eqs(t,x,u)
Real t;
Matrix x,u;
{
    Matrix dxp;
    Real m,l,M,f,J,c,a,c1,c2,g;
    Real r,th,dr,dth;

    m=0.030;
    l=0.16;
    M=0.73;
    f=7.6;
    J=2.5e-4;
    c=4.7e-5;
    a=0.62;
    c1=1.0;
    c2=1.0;
    g=9.8;

    r=x(1,1);
    th=x(2,1);
    dr=x(3,1);
    dth=x(4,1);

    //倒立振子の状態の微分（非線形モデル）
    dxp = [[dr]
            [dth]
         [[[m+M,m*l*cos(th)]
           [m*l*cos(th),J+m*l*l]]~*[[a*u(1,1)-f*dr+m*g*dth*dth*sin(th)]
                                            [m*g*l*sin(th)-c*dth]]]];

    return dxp;
}
Func Matrix link_eqs(t,x)
Real t;
Matrix x;
{
    Matrix u,xref;

    //台車の可動範囲に関する制限
    if(x(1,1) <= -0.16 || 0.16 <= x(1,1)){
        OdeStop();
    }

    xref = [0 0 0 0]';//'
    u = F*(xref - x);

    //入力の大きさに関する制限
    if(u(1,1) <= -15){
        u(1,1) = -15;
    }else if (u(1,1) >= 15){
        u(1,1) = 15;
    }
    return u;
}