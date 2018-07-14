Matrix  z,A,B,C,xref;
Real dt,cref;
Matrix Ah,Bh,Jh,Ch,Dh,F;
Matrix Ahd,Bhd,Jhd,Hhd;
CoMatrix obs_p;
Matrix Q,P,R;
Func void main()
{
    Real t0,t1,r0,th0,tol,dtsav;
    Matrix x0,T,X,U,data,n;
    Matrix diff_eqs(),link_eqs();

    read A << "A.mat";
    read B << "B.mat";
    read C << "C.mat";

    read data << "data857";
    //Fの設計
    Q=diag(1E4,1E3,1E0,1E0);
    R=[1];
    {F,P}=lqr(A,B,Q,R);
    //オブザーバの設計
    obs_p = trans([(-100,0),(-100,0)]);
    {Ah, Bh, Ch, Dh,Jh}=obsg(A,B,C,obs_p);
    dt = 0.01;//サンプリング周期

    t0 = 0.0;
    t1 = 16.0;
    r0 = 0.0;
    th0 =data(1,4);
    x0 = trans([r0 th0/180*PI 0 0]); //倒立振子の初期状態
    z = trans([0 0]); //オブザーバの初期状態
    tol = 1E-8;
    dtsav = 0.01;
    xref = trans([0.1 0 0 0]);
    cref = -0.1;

    {Ahd,Hhd} = c2d(Ah, [Bh Jh], dt);
    Bhd = Hhd(:,1:2);
    Jhd = Hhd(:,3);

    {T,X,U} = Ode45HybridAuto(t0,t1,dt,x0,diff_eqs,link_eqs,tol,dtsav);


    mgplot(1,data(1,:),data(3,:),{"actual(r)"});
    mgreplot(1,T,X(1,:),{"simulation(r)"});
    mgplot_xlabel(1,"t[s]");
    mgplot_ylabel(1,"r[m]");
    mgplot(2,data(1,:),data(4,:),{"actual(th)"});
    mgreplot(2,T,X(2,:),{"simulation(th)"});
    mgplot_xlabel(2,"th[rad]");
    mgplot_ylabel(2,"t[s]");

    mgplot_eps(1,"r857");
    mgplot_eps(2,"th857");
}
Func Matrix diff_eqs(t,x,u)
    Matrix x,u;
    Real t;
{
    Matrix dx;
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
    dx = [[dr]
            [dth]
         [[[m+M,m*l*cos(th)]
           [m*l*cos(th),J+m*l*l]]~*[[a*u(1,1)-f*dr+m*g*dth*dth*sin(th)]
                                            [m*g*l*sin(th)-c*dth]]]];

    return dx;
}
Func Matrix link_eqs(t,x)
    Matrix x;
    Real t;
{
    Matrix u,y,xh;

    //台車の可動範囲に関する制限
    if(x(1,1) <= -0.16 || 0.16 <= x(1,1)){
        OdeStop();
    }

    y = C*x;
    xh = Ch*z + Dh*y;

    if(rem(t,5)==0){
        xref(1,1) = xref(1,1)+cref;
        cref = -cref;
    }

    u = F*(xref - xh);

    //入力の大きさに関する制限
    if(u(1,1) <= -15){
        u(1,1) = -15;
    }else if (u(1,1) >= 15){
        u(1,1) = 15;
    }

    z = Ahd*z + Bhd*y + Jhd*u;//オブザーバの状態の更新

    return u;
}