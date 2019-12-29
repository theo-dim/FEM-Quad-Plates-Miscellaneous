      function [P]=princ(S,nsd)
%-----------------------------------------------------------------------
%.... PROGRAM TO COMPUTE PRINCIPAL VALUES OF SYMMETRIC 2ND-RANK TENSOR
%
%        S = SYMMETRIC SECOND-RANK TENSOR STORED AS A VECTOR
%      nsd = NUMBER OF DIMENSIONS (2 OR 3)
%        P = PRINCIPAL VALUES
%
%.... THE COMPONENTS OF S MUST BE STORED IN THE FOLLOWING ORDERS
%
%        2-D PROBLEMS, S11,S22,S12
%        3-D PROBLEMS, S11,S22,S33,S12,S23,S31
%-----------------------------------------------------------------------
      RT2=1.41421356237309;
      PI23=2.09439510239321;
      P=zeros(4,1);
 
%.... 2-D PROBLEM

      if( nsd == 2 )
        a=22.5/atan(1.0);
        X=0.5*(S(1)+S(2));
        Y=0.5*(S(1)-S(2));
        R=sqrt(Y*Y+S(3)*S(3));
        P(1)=X+R;
        P(2)=X-R;
        P(3)=R;
        if( Y == 0 | S(3) == 0 )
          P(4)=45.0;
        else
          P(4)=a*atan(S(3)/Y);
        end

%.... 3-D PROBLEM

      elseif( nsd == 3)
        R=0;
        X=(S(1)+S(2)+S(3))/3.;
        Y=S(1)*(S(2)+S(3))+S(2)*S(3)-S(4)*S(4)-S(6)*S(6)-S(5)*S(5);
        Z=S(1)*S(2)*S(3)+2.d0*S(4)*S(6)*S(5)-S(1)*S(5)*S(5)-S(2)*S(6)*S(6)-S(3)*S(4)*S(4);
        T=3*X*X-Y;
        U=0;
        if ( T ~= 0 )
           U=sqrt(2*T/3);
           UCUBED=U*U*U;
           if( UCUBED == 0 )
              U=0;
           else
              a=(Z+(T-X*X)*X)*RT2/UCUBED;
              R=sqrt(abs(1-a*a));
              if(R == 0 | a == 0 )
                R=0;
              else
                R=atan(R/a)/3;
              end
            end
        end
        P(1)=X+U*RT2*cos(R);
        P(2)=X+U*RT2*cos(R-PI23);
        P(3)=X+U*RT2*cos(R+PI23);

      end
