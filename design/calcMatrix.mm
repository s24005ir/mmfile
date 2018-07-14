Func void main()
{
    Matrix K,A,A21,A22,B,B2,C,D,E;
    Matrix Nc,No;
    Real m,l,M,f,J,c,a,g,c1,c2;

    m=0.030;
    l=0.155;
    a=0.624;
    g=9.8;
    M=0.713;
    f=7.57;
    J=2.51e-4;
    c=4.67e-5;
    c1=1.0;
    c2=1.0;

    K=[[M+m,m*l]
        [m*l,J+m*l*l]];

    A21=K~*[[0 0]
            [0 m*g*l]];
    
    A22=K~*[[-f 0]
            [0 ,-c]];

    A=[[Z(2,2) I(2)]
        [A21 A22]];
    
    B2=[a 0]';//'

    B=[[0][0][K~*B2]];

    C=[[c1 0 0 0]
        [0 c2 0 0]];
    
    Nc=[B A*B A*A*B A*A*A*B];
    No=[[C] [C*A] [C*A*A] [C*A*A*A]];

    print A >> "A.mat";
    print A;
    print eigval(A);
    print B >> "B.mat";
    print B;
    print C >> "C.mat";
    print C;
    print Nc;
    print rank(Nc);
    print No;
    print rank(No);
}