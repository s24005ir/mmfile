Matrix  z,A,B,C;
Real dt;
Matrix Ah,Bh,Jh,Ch,Dh,F;
Matrix Ahd,Bhd,Jhd,Hhd;
Func void main()
{
    Real t0,t1,r0,th0,tol,dtsav;
    Matrix x0,T,X,U;
    CoMatrix pc,obs_p;
    Matrix diff_eqs(),link_eqs();

    t0 = 0.0;
    t1 = 5.0;
    r0 = 0.0;
    th0 = 15.0;
    x0 = trans([r0 th0/180*PI 0 0]); //倒立振子の初期状態
    z = trans([0 0]); //オブザーバの初期状態
    dt = 0.005;
    tol = 1E-8;
    dtsav = 0.05;

    read A << "A.mat";
    read B << "B.mat";
    read C << "C.mat";

    //Fの設計
    pc = [(-100,0),(-100,0),(-2,0),(-2,0)]'; //'
    F = pplace(A,B,pc);
    
    //オブザーバの設計
    obs_p = trans([(-10,0),(-10,0)]);
    {Ah, Bh, Ch, Dh,Jh}=obsg(A,B,C,obs_p);

    {T,X,U} = Ode45HybridAuto(t0,t1,dt,x0,diff_eqs,link_eqs,tol,dtsav);

    mgplot(1,T,X(1:2,:));

    print [[T][X][U]] >> "discreteObsPP.mat";
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
    Matrix u,y,xh,xref;

    //台車の可動範囲に関する制限
    if(x(1,1) <= -0.16 || 0.16 <= x(1,1)){
        OdeStop();
    }

    {Ahd,Hhd} = c2d(Ah, [Bh Jh], dt);
    Bhd = Hhd(:,1:2);
    Jhd = Hhd(:,3);

    y = C*x;
    xh = Ch*z + Dh*y;
    xref = trans([0 0 0 0]);
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