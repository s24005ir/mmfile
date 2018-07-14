Matrix A,B,C;
Matrix Ah,Bh,Ch,Dh,Jh;
    CoMatrix obs_p;
    Matrix Q,R,P,F;
Func void main()
{
    Real t0,t1,r0,th0,tol;
    Matrix xp0,z0,x0,T,X,U;
    Matrix diff_eqs(),link_eqs();


    t0=0.0;
    t1=5.0;
    r0=0.0;
    th0=15.0;
    xp0=trans([r0 th0/180*PI 0 0]);
    z0=trans([0 0]);
    x0=[[xp0][z0]];
    tol = 1.0E-8;
    
    print "simulating\n";

    //LQ最適制御による設計
    Q=diag(1000,10,1,1);
    R=[1];
    {F,P}=lqr(A,B,Q,R);

    obs_p = trans([(-2,0),(-2,0)]);
    {Ah, Bh, Ch, Dh,Jh}=obsg(A,B,C,obs_p);

    read A << "A.mat";
    read B << "B.mat";
    read C << "C.mat";

    {T,X,U}=Ode45Auto(t0,t1,x0,diff_eqs,link_eqs,tol);

    mgplot(1,T,X(1,:),{"r"});
    mgplot(2,T,X(2,:),{"th"});
    mgplot(3,T,X(3,:),{"dr"});
    mgplot(4,T,X(4,:),{"dth"});
    mgplot(5,T,U(1,:),{"u"});

    print [[T][X][U]] >> "observerDesign.mat";
}

Func Matrix diff_eqs(t,x,u)
Real t;
Matrix x,u;
{
    Matrix xp,y,z,dx,dxp,dz;
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

    xp = x(1:4,1);
    z = x(5:6,1);
    y = C*xp;

    r=xp(1,1);
    th=xp(2,1);
    dr=xp(3,1);
    dth=xp(4,1);

    //倒立振子の状態の微分（非線形モデル）
    dxp = [[dr]
            [dth]
         [[[m+M,m*l*cos(th)]
           [m*l*cos(th),J+m*l*l]]~*[[a*u(1,1)-f*dr+m*g*dth*dth*sin(th)]
                                            [m*g*l*sin(th)-c*dth]]]];
    
    dz = Ah*z+Bh*y+Jh*u;

    dx = [[dxp][dz]];

    return dx;
}
Func Matrix link_eqs(t,x)
Real t;
Matrix x;
{
    Matrix xp,z,u,y,xref,xh;

    xp=x(1:4,1); //倒立振子の状態
    z = x(5:6,1);//オブザーバの状態



    //台車の可動範囲に関する制限
    if(xp(1,1) <= -0.16 || 0.16 <= xp(1,1)){
        OdeStop();
    }

    y = C*xp; //出力の計算

    xh = Ch*z + Dh*y; //状態の推定値
    xref = trans([0 0 0 0]); //状態の目標値
    u = F*(xref - xh);

    //入力の大きさに関する制限
    if(u(1,1) <= -15){
        u(1,1) = -15;
    }else if (u(1,1) >= 15){
        u(1,1) = 15;
    }
    return u;
}