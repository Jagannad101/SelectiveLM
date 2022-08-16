% global NI %Number of nodes in X direction
% global NJ %Number of nodes in Y direction
% global NK %Number of noded in Z direction
% 
% global L1 %I value of the last node in the X direction
% global L2
% global L3
% 
% global M1 %J value of the last node in the Y direction
% global M2
% global M3
% 
% global N1 %K value of the last node in the Z direction
% global N2
% global N3

% global DT

% global LAST
% 
% LAST = 1000;



%global XCV


% NIJ = NI;
% NFMAX = 10;
% NFX4 = NFMAX + 4;

% global EPS



% global RHO CON AP CP

% PC U V W P T 

% global mu mus mul
% 
%  mus = 1.0E4;
%  mul = 5.0E-3;

% global XU XP XDIF XCV XCVS XCVI XCVIP
% global YV YP YDIF YCV YCVS YCVJ YCVJP
% global ZW ZP ZDIF ZCV ZCVS ZCVK ZCVKP
% global AX AY AZ
% global FX FY FZ FXM FYM FZM
% global RELAX
% global TIME ITER




 






% global Ub %Laser velocity
% 
% for i = 1:L1
%     for j = 1:M1
%         for k = 1:N1
%             Ub(i,j,k) = 0.2; %Set the value of laser speed
%         end
%     end
% end


%----------------------------------------------





% for i = 1:L1
%     for j = 1:M1
%         for k = 1:N1
% 
%             RHO(i,j,k) = 1990;
%         end
%     end
% end
  

%------------------Z-------------------------------- 

% Initialize_grid();
%  
% Total_Variation_Diminishing();

%  dbstop if naninf;

[U,V,W,P,T, Res_max, Res_U, Res_V, Res_W, Res_T, ITER, TIME] = Total_Variation_Diminishing();

if(Res_max < 1.0e-6)
    isConverged = 1;
else
    isConverged = 0;
end

if isConverged == 0
    fprintf('Solution failed to converge in %d iterations!!!', 1000);
end

if isConverged == 1
    fprintf('Solution converged in %d iterations!!!', ITER);
end

P_new = P;

[P] = Extrapolate(P_new);

% uc = 0.5*( u( 2 : end-1, 2: end-1) + u( 3 : end, 2: end-1));
% vc = 0.5*( v( 2 : end-1, 2: end-1) + v( 2 : end-1, 3: end));

%disp(P_new);

%  fp1 = fopen('velocity.dat' , 'w' ) ;
% fprintf (fp1, 'TITLE = "Field Data"\n' ) ;
% fprintf (fp1, ' variables="x (m)", "y (m)", "z(m)", "u(m/s)", "v(m/s)", "w(m/s)"\n');
% % fprintf (fp1, ' zone T="n=%d"\n' , Nc) ;
% fprintf (fp1, ' I= %d J= %d K=%d\n' ,202, 92, 92) ;
% for i =1: 22
%     for j = 1:12
%         for k = 1:12
% 
%             if (i == 1)
%                 uc(i,j,k) = U(2,j,k);
%             elseif (i == 22)
%                 uc(i,j,k) = U(22,j,k);
%             else
%                 uc(i,j,k) = 0.5*(U(i,j,k) + U(i+1,j,k));
%             end
% 
%              if (j == 1)
%                 vc(i,j,k) = V(i,2,k);
%             elseif (j == 12)
%                 vc(i,j,k) = V(i,12,k);
%             else
%                 vc(i,j,k) = 0.5*(V(i,j,k) + V(i,j+1,k));
%              end
% 
%              if (k == 1)
%                 wc(i,j,k) = W(i,j,2);
%             elseif (k == 12)
%                 wc(i,j,k) = W(i,j,12);
%             else
%                 wc(i,j,k) = 0.5*(W(i,j,k) + W(i,j,k+1));
%              end
% 
% 
% 
% 
% 
% 
%    
%         fprintf (fp1, '%e, %e, %e, %e, %e, %e\n' , i , j, k,  uc(i,j,k) , vc(i,j,k), wc(i,j,k) ) ;
% 
%         end
%     end
% end
% fclose(fp1);
% 
% fp2 = fopen('temperature.dat' , 'w' ) ;
% fprintf (fp2, 'TITLE = "Field temperature"\n' ) ;
% fprintf (fp2, ' variables="x (m)", "y (m)", "z(m)", "T(k)"\n');
% % fprintf (fp1, ' zone T="n=%d"\n' , Nc) ;
% fprintf (fp2, ' I= %d J= %d K=%d\n' ,202, 92, 92) ;
% for i =1: 22
%     for j = 1:12
%         for k = 1:12
% 
%             fprintf (fp2, '%e, %e, %e, %e\n' , i , j, k, T(i,j,k)) ;
% 
%         end
%     end
% end
% fclose(fp2);
% 
% fp4 = fopen('pressure.dat' , 'w' ) ;
% fprintf (fp4, 'TITLE = "Field pressure"\n' ) ;
% fprintf (fp4, ' variables="x (m)", "y (m)", "z(m)", "P(N/m^2)"\n');
% % fprintf (fp1, ' zone T="n=%d"\n' , Nc) ;
% fprintf (fp4, ' I= %d J= %d K=%d\n' ,22, 12, 12) ;
% for i =1: 22
%     for j = 1:12
%         for k = 1:12
% 
%             fprintf (fp4, '%e, %e, %e, %e\n' , i , j, k, P(i,j,k)) ;
% 
%         end
%     end
% end
% fclose(fp4);
%  
%  

 
    
%%%%%%%%%%-------------Initialize the grid and related parameters here---------------%%%%%%%%%%
function [XDIF, XCV, XCVS, XCVI, XCVIP, YDIF, YCV, YCVS, YCVJ, YCVJP, ZDIF, ZCV, ZCVS, ZCVK, ZCVKP, AX, AY, AZ, FX, FXM, FY, FYM, FZ, FZM, XP, YP, ZP] = Initialize_grid()

NI = 20;
NJ = 10;
NK = 10;

L1 = 22;
L2 = L1 - 1;
L3 = L2 - 1;

M1 = 12;
M2 = M1 - 1;
M3 = M2 - 1;

N1 = 12;
N2 = N1 - 1;
N3 = N2 - 1;




XP = zeros(L1,1);
XU = zeros(L1,1);
XDIF = zeros(L1,1);
XCV = zeros(L1,1);
XCVS = zeros(L1,1);
XCVI = zeros(L1,1);
XCVIP = zeros(L1,1);

YP = zeros(M1,1);
YV = zeros(M1,1);
YDIF = zeros(M1,1);
YCV = zeros(M1,1);
YCVS = zeros(M1,1);
YCVJ = zeros(M1,1);
YCVJP = zeros(M1,1);

ZP = zeros(N1,1);
ZW = zeros(N1,1);
ZDIF = zeros(N1,1);
ZCV = zeros(N1,1);
ZCVS = zeros(N1,1);
ZCVK = zeros(N1,1);
ZCVKP = zeros(N1,1);

AX = zeros(M1,N1);
AY = zeros(L1,N1);
AZ = zeros(L1,M1);

FX = zeros(L1,1);
FXM = zeros(L1,1);
FY = zeros(M1,1);
FYM = zeros(M1,1);
FZ = zeros(N1,1);
FZM = zeros(N1,1);

XL = 0.0036;
YL = 0.0012;
ZL = 0.0012;

XU(2) = 0;
DX = XL/(L1-2);
for i = 3:L1
    XU(i) = XU(i-1) + DX;
end

YV(2) = 0;
DY = YL/(M1-2);
for j = 3:M1
    YV(j) = YV(j-1) + DY;
end

ZW(2) = 0;
DZ = ZL/(N1-2);
for k = 3:N1
    ZW(k) = ZW(k-1) + DZ;
end




% imin = 2; 
% imax = imin + NI - 1;
% jmin = 2; 
% jmax = jmin + NJ - 1;
% kmin = 2;
% kmax = kmin + NK -1;
% 
% %%%%%%%%%%%%%%%%%%%% Node Positions in X/Y/Z directions %%%%%%%%%%%%%
% 
% XU ( imin : imax+1)=linspace ( 0 , 3.6 , NI +1); %Nodes storing U velocities
% YV ( jmin : jmax+1)=linspace ( 0 , 1.2 , NJ +1); %Nodes storing V velocity
% ZW ( kmin : kmax+1)=linspace ( 0 , 1.2 , NK +1); %Nodes storing W velocity
% XP( imin : imax )= 0.5*( XU ( imin : imax)+ XU ( imin +1: imax + 1 ) ); %Scalar nodes in the X direction
% YP( jmin : jmax )= 0.5*( YV ( jmin : jmax)+ YV ( jmin +1: jmax + 1 ) ); %Scalar nodes in the Y direction
% ZP( kmin : kmax )= 0.5*( ZW ( kmin : kmax)+ ZW ( kmin +1: kmax + 1 ) ); %Scalar nodes in the Z direction
% XP(imax + 1) = XU(imax + 1);
% YP(jmax + 1) = YV(jmax + 1);
% ZP(kmax + 1) = ZW(kmax + 1);
%%%%%% Width of the main C.V and C.V of velocity in X/Y/Z directions %%%%

XP(1) = XU(2);
for i = 2:L2
    XP(i) = 0.5 * (XU(i+1) + XU(i));
end
XP(L1) = XU(L1);

YP(1) = YV(2);
for j = 2:M2
    YP(j) = 0.5 * (YV(j+1) + YV(j));
end
YP(M1) = YV(M1);

ZP(1) = ZW(2);
for k = 2:N2
    ZP(k) = 0.5 * (ZW(k+1) + ZW(k));
end
ZP(N1) = ZW(N1);



for i = 2:L1
    XDIF(i) = XP(i) - XP(i-1); % Distance between grid nodes in the X direction
end
%----------------------------------------------------------------
for i = 2:L2
    XCV(i) = XU(i+1) - XU(i); % The width of the main control volume in the X direction
end

for i = 3:L2
    XCVS(i) = XDIF(i);
end
XCVS(3) = XCVS(3) + XDIF(2);
XCVS(L2) = XCVS(L2) + XDIF(L1);

for i = 3:L3
    XCVI(i) = 0.5*XCV(i);
    XCVIP(i) = XCV(i);
end
XCVIP(2) = XCV(2);
XCVI(L2) = XCV(L2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
for j = 2:M1
    YDIF(j) = YP(j) - YP(j-1); % Distance between grid nodes in the Y direction
end
%-------------------------------------------------------------------
for j = 2:M2
    YCV(j) = YV(j+1) - YV(j); % The width of the main control volume in the Y direction
end

for j = 3:M2
    YCVS(j) = YDIF(j);
end
YCVS(3) = YCVS(3) + YDIF(2);
YCVS(M2) = YCVS(M2) + YDIF(M1);

for j = 3:M3
    YCVJ(j) = 0.5*YCV(j);
    YCVJP(j) = YCVJ(j);
end
YCVJP(2) = YCV(2);
YCVJ(M2) = YCV(M2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 2:N1
    ZDIF(k) = ZP(k) - ZP(k-1); % Distance between grid nodes in the Z direction
end
%-----------------------------------------------------------------
for k = 2:N2
    ZCV(k) = ZW(k+1) - ZW(k); % The width of the main control volume in the Z direction
end

for k = 3:N2
    ZCVS(k) = ZDIF(k);
end
ZCVS(3) = ZCVS(3) + ZDIF(2);
ZCVS(N2) = ZCVS(N2) + ZDIF(N1);

for k = 3:N3
    ZCVK(k) = 0.5*ZCV(k);
    ZCVKP(k) = ZCVK(k);
end
ZCVKP(2) = ZCV(2);
ZCVK(N2) = ZCV(N2);




%%%%%%%%%%%%%%%%%%% Area and mass flow interpolation %%%%%%%%%%%%%%

for k = 2:N2
    for j = 2:M2
        AX(j,k) = YCV(j) * ZCV(k);
    end
    for i = 2:L2
        AY(i,k) = XCV(i) * ZCV(k);
    end
end

for j = 2:M2
    for i = 2:L2
        AZ(i,j) = XCV(i) * YCV(j);
    end
end

%-----------------------------------------------------------------
 
for i = 3:L2
    FX(i) = 0.5 * XCV(i-1)/XDIF(i);
    FXM(i) = 1 - FX(i);
end

FX(2) = 0;
FXM(2) = 1;
FX(L1) = 1;
FXM(L1) = 0;

%------------------------------------------------------------------

for j = 3:M2
    FY(j) = 0.5 * YCV(j-1)/YDIF(j);
    FYM(j) = 1 - FY(j);
end

FY(2) = 0;
FYM(2) = 1;
FY(M1) = 1;
FYM(M1) = 0;

%---------------------------------------------------------------------

for k = 3:N2
    FZ(k) = 0.5 * ZCV(k-1)/ZDIF(k);
    FZM(k) = 1 - FZ(k);
end

FZ(2) = 0;
FZM(2) = 1;
FZ(N1) = 1;
FZM(N1) = 0;

%---------------------------------------------------------------------
 %disp(YP)
 %disp(XU)

 %Initialize fields--------------------------


 end
 %--------------------------------------------------------------------

 



 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,V,W,P,T, Res_max, Res_U, Res_V, Res_W, Res_T, ITER, TIME] = Total_Variation_Diminishing()

NI = 20;
NJ = 10;
NK = 10;

L1 = 22;
L2 = L1 - 1;
L3 = L2 - 1;

M1 = 12;
M2 = M1 - 1;
M3 = M2 - 1;

N1 = 12;
N2 = N1 - 1;
N3 = N2 - 1;

RELAX = 0.6;


[XDIF, XCV, XCVS, XCVI, XCVIP, YDIF, YCV, YCVS, YCVJ, YCVJP, ZDIF, ZCV, ZCVS, ZCVK, ZCVKP, AX, AY, AZ, FX, FXM, FY, FYM, FZ, FZM, XP, YP, ZP] = Initialize_grid();

PC = zeros(L1,M1,N1); %Pressure correction matrix
U = zeros(L1,M1,N1); %U velocity matrix
V = zeros(L1,M1,N1); %V velocity matrix
W = zeros(L1,M1,N1); %W velocity matrix
P = zeros(L1,M1,N1); %Pressure matrix
T = 298.15 .* ones(L1,M1,N1); %Temperature matrix

Fs = zeros(L1,M1,N1);
Fn = zeros(L1,M1,N1);
Fe = zeros(L1,M1,N1);
Fw = zeros(L1,M1,N1);
Ff = zeros(L1,M1,N1);
Fb = zeros(L1,M1,N1);
FLOWSUM = zeros(L1,M1,N1);

REP = zeros(L1,M1,N1);
REM = zeros(L1,M1,N1);
RWP = zeros(L1,M1,N1);
RWM = zeros(L1,M1,N1);
RNP = zeros(L1,M1,N1);
RNM = zeros(L1,M1,N1);
RSP = zeros(L1,M1,N1);
RSM = zeros(L1,M1,N1);
RBP = zeros(L1,M1,N1);
RBM = zeros(L1,M1,N1);
RFP = zeros(L1,M1,N1);
RFM = zeros(L1,M1,N1);

alphae = zeros(L1,M1,N1);
alphaw = zeros(L1,M1,N1);
alphas = zeros(L1,M1,N1);
alphan = zeros(L1,M1,N1);
alphaf = zeros(L1,M1,N1);
alphab = zeros(L1,M1,N1);


Ae = zeros(L1,M1,N1);
Aw = zeros(L1,M1,N1);
As = zeros(L1,M1,N1);
An = zeros(L1,M1,N1);
Af = zeros(L1,M1,N1);
Ab = zeros(L1,M1,N1);

SUDC = zeros(L1,M1,N1);

COFU = zeros(L1,M1,N1,8);
COFV = zeros(L1,M1,N1,8);
COFW = zeros(L1,M1,N1,8);
COFP = zeros(L1,M1,N1,8);
COFT = zeros(L1,M1,N1,8);

RELAX = 0.6;
REL = 1 - RELAX;

DT = 1.0E-3;

TIME = DT;
ITER = 1;

LAST = 1000;

Res_max = 10000;
toler = 1.0e-6;
% cond = 10000;

% toler_time = 1.0e-3;
% 
% TIME_MAX = 2;
% 
% U_old = U;
% V_old = V;
% W_old = W;
% T_old = T;

% RHO = (1990 * (1-0.2)) .* ones(L1,M1,N1);
% CP = (1211 * (1-0.2)) .* ones(L1,M1,N1);

% while (cond > toler_time || TIME <= TIME_MAX)






while ( Res_max > toler || ITER <= LAST)

[~,~,RHO, CP,~,~, U,V,W,T] = Boundary_Conditions (U,V,W,T,ITER,LAST);



if (ITER <= LAST - 650)

NF = 1;
IST = 3;
JST = 2;
KST = 2;

[GAM,CON,AP] = Diffusion_Source(U, V, W, T, NF,ITER,LAST) ;


    for k = 3:N3
        for j = 3:M3
            for i = 3:L2
                if abs(U(i+1,j,k)-U(i,j,k)) > 1.0E-20
                        REP(i,j,k) = (U(i,j,k) - U(i-1,j,k))/(U(i+1,j,k) - U(i,j,k));
                    if i == L2
                        REM(i,j,k) = 1;
                    else
                        REM(i,j,k) = (U(i+2,j,k) - U(i+1,j,k))/(U(i+1,j,k) - U(i,j,k));
                    end
                else
                    REP(i,j,k) = 1;
                    REM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(U(i,j,k)-U(i-1,j,k)) > 1.0E-20
                    if i ==3
                        RWP(i,j,k) = 1;
                    else
                        RWP(i,j,k) = (U(i-1,j,k) - U(i-2,j,k))/(U(i,j,k) - U(i-1,j,k));
                    end
                    RWM(i,j,k) = (U(i+1,j,k) - U(i,j,k))/(U(i,j,k) - U(i-1,j,k));
                else
                    RWP(i,j,k) = 1;
                    RWM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------


               if abs(U(i,j+1,k)-U(i,j,k)) > 1.0E-20
                    RBP(i,j,k) = (U(i,j,k) - U(i,j-1,k))/(U(i,j+1,k) - U(i,j,k));
                    if j == M3
                        RBM(i,j,k) = 1;
                    else
                        RBM(i,j,k) = (U(i,j+2,k) - U(i,j+1,k))/(U(i,j+1,k) - U(i,j,k));
                    end
                else
                    RBP(i,j,k) = 1;
                    RBM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(U(i,j,k)-U(i,j-1,k)) > 1.0E-20
                    if j ==3
                        RFP(i,j,k) = 1;
                    else
                        RFP(i,j,k) = (U(i,j-1,k) - U(i,j-2,k))/(U(i,j,k) - U(i,j-1,k));
                    end
                    RFM(i,j,k) = (U(i,j+1,k) - U(i,j,k))/(U(i,j,k) - U(i,j-1,k));
                else
                    RFP(i,j,k) = 1;
                    RFM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------

                 if abs(U(i,j,k+1)-U(i,j,k)) > 1.0E-20
                    RNP(i,j,k) = (U(i,j,k) - U(i,j,k-1))/(U(i,j,k+1) - U(i,j,k));
                    if k == N3
                        RNM(i,j,k) = 1;
                    else
                        RNM(i,j,k) = (U(i,j,k+2) - U(i,j,k+1))/(U(i,j,k+1) - U(i,j,k));
                    end
                else
                    RNP(i,j,k) = 1;
                    RNM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(U(i,j,k)-U(i,j,k-1)) > 1.0E-20
                    if k ==3
                        RSP(i,j,k) = 1;
                    else
                        RSP(i,j,k) = (U(i,j,k-1) - U(i,j,k-2))/(U(i,j,k) - U(i,j,k-1));
                    end
                    RSM(i,j,k) = (U(i,j,k+1) - U(i,j,k))/(U(i,j,k) - U(i,j,k-1));
                else
                    RSP(i,j,k) = 1;
                    RSM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------
            end
        end
    end
%-----------------------------------------------------------------63
    k = 2;
    for j = 2:M2
        for i = 3:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end


    k = N2;
    for j = 2:M2
        for i = 3:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end
%------------------------------------------------------------
      j = 2;
    for k = 2:N2
        for i = 3:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    j = M2;
    for k = 2:N2
        for i = 3:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    for k = 2:N2
        for i = 3:L2 
            FL = XCVI(i) * V(i,2,k) * RHO(i,1,k);
            FLM = XCVIP(i-1) * V(i-1,2,k) * RHO(i-1,1,k);
            FLOW = ZCV(k) * (FL + FLM);
            DIFF = ZCV(k) * (XCVI(i) * GAM(i,1,k) + XCVIP(i-1) * GAM(i-1,1,k))/YDIF(2);
            Af(i,2,k) = DIFF + max(0,FLOW);
            Ff(i,2,k) = FLOW;
        end
    end 

    for j = 2:M2
        for i = 3:L2
            FL = XCVI(i) * W(i,j,2) * RHO(i,j,1);
            FLM = XCVIP(i-1) * W(i-1,j,2) * RHO(i-1,j,1);
            FLOW = YCV(j) * (FL + FLM);
            DIFF = YCV(j) * (XCVI(i) * GAM(i,j,1) + XCVIP(i-1) * GAM(i-1,j,1))/ZDIF(2);
            As(i,j,2) = DIFF + max(0,FLOW);
            Fs(i,j,2) = FLOW;
        end
    end

    for k = 2:N2
        for j = 2:M2
            FLOW = AX(j,k) * U(2,j,k) * RHO(i,j,k);
            DIFF = AX(j,k) * GAM(1,j,k)/XCV(2);
            Aw(3,j,k) = DIFF + max(0,FLOW);
            Fw(3,j,k) = FLOW;

            for i = 3:L2
                if i == L2
                    FLOW = AX(j,k) * U(L1,j,k) * RHO(L1,j,k);
                    DIFF = AX(j,k) * GAM(L1,j,k)/XCV(L2);

                else
                    FL = U(i,j,k) * (FX(i) * RHO(i,j,k) + FXM(i)*RHO(i-1,j,k));
                    FLP = U(i+1,j,k) * (FX(i+1) * RHO(i+1,j,k) + FXM(i+1)*RHO(i,j,k));
                    FLOW = AX(j,k) * 0.5 * (FL + FLP);
                    DIFF = AX(j,k) * GAM(i,j,k)/XCV(i);
                end
                Aw(i+1,j,k) = DIFF + max(0,FLOW);
                Ae(i,j,k) = DIFF + max(0,-FLOW); %Different in his code
                Fw(i+1,j,k) = FLOW;
                Fe(i,j,k) = FLOW;

                if j == M2
                    FL = XCVI(i) * V(i,M1,k) * RHO(i,M1,k);
                    FLM = XCVIP(i-1) * V(i-1,M1,k) * RHO(i-1,M1,k);
                    DIFF =  ZCV(k) * (XCVI(i) * GAM(i,M1,k) + XCVIP(i-1) * GAM(i-1,M1,k))/YDIF(M1);
                    FLOW = ZCV(k) * (FL + FLM);
                else
                    FL = XCVI(i) * V(i,j+1,k) * (FY(j+1) * RHO(i,j+1,k) + FYM(j+1) * RHO(i,j,k));
                    FLM = XCVIP(i-1) * V(i-1,j+1,k) * (FY(j+1) * RHO(i-1,j+1,k) + FYM(j+1) * RHO(i-1,j,k));
                    FLOW = ZCV(k) * (FL + FLM);
                    GM = GAM(i,j,k) * GAM(i,j+1,k)/(YCV(j) * GAM(i,j+1,k) + YCV(j+1) * GAM(i,j,k)) * XCV(i);
                    GMM = GAM(i-1,j,k) * GAM(i-1,j+1,k)/(YCV(j) * GAM(i-1,j+1,k) + YCV(j+1)*GAM(i-1,j,k)) * XCVIP(i-1);
                    DIFF = ZCV(k) * 2 .* (GM + GMM);
                end
                Af(i,j+1,k) = DIFF + max(0,FLOW);
                Ab(i,j,k) = DIFF + max(0,-FLOW);
                Ff(i,j+1,k) = FLOW;
                Fb(i,j,k) = FLOW;

                if k == N2
                    FL = XCVI(i) * W(i,j,N1) * RHO(i,j,N1);
                    FLM = XCVIP(i-1) * W(i-1,j,N1) * RHO(i-1,j,N1);
                    DIFF =  YCV(j) * (XCVI(i) * GAM(i,j,N1) + XCVIP(i-1) * GAM(i-1,j,N1))/ZDIF(N1);
                    FLOW = YCV(j) * (FL + FLM);
                else
                    FL = XCVI(i) * W(i,j,k+1) * (FZ(k+1) * RHO(i,j,k+1) + FZM(k+1) * RHO(i,j,k));
                    FLM = XCVIP(i-1) * W(i-1,j,k+1) * (FZ(k+1) * RHO(i-1,j,k+1) + FZM(k+1) * RHO(i-1,j,k));
                    FLOW = YCV(j) * (FL + FLM);
                    GM = GAM(i,j,k) * GAM(i,j,k+1)/(ZCV(k) * GAM(i,j,k+1) + ZCV(k+1) * GAM(i,j,k)) * XCV(i);
                    GMM = GAM(i-1,j,k) * GAM(i-1,j,k+1)/(ZCV(k) * GAM(i-1,j,k+1) + ZCV(k+1)*GAM(i-1,j,k)) * XCVIP(i-1);
                    DIFF = YCV(j) * 2 .* (GM + GMM);
                end
                As(i,j,k+1) = DIFF + max(0,FLOW);
                An(i,j,k) = DIFF + max(0,-FLOW);
                Fs(i,j,k+1) = FLOW;
                Fn(i,j,k) = FLOW;

                VOL = AX(j,k) * XCVS(i);
                APT = (RHO(i,j,k)*XCVI(i) + RHO(i-1,j,k)*XCVIP(i-1))/(XCVS(i)*DT);
                AP(i,j,k) = AP(i,j,k) - APT;
                CON(i,j,k) = CON(i,j,k) + APT*U(i,j,k);
                FLOWSUM(i,j,k) = Fe(i,j,k) - Fw(i,j,k) + Fn(i,j,k) - Fs(i,j,k) + Fb(i,j,k) - Ff(i,j,k);
                AP(i,j,k) = (-AP(i,j,k)*VOL + FLOWSUM(i,j,k) + Aw(i,j,k) + Ae(i,j,k) + Af(i,j,k) + Ab(i,j,k) + As(i,j,k) + An(i,j,k))/RELAX;

                alphaw(i,j,k) = 0;
                if (Fw(i,j,k) > 0)
                    alphaw(i,j,k) = 1;
                end

                alphae(i,j,k) = 0;
                if (Fe(i,j,k) > 0)
                    alphae(i,j,k) = 1;
                end

                alphan(i,j,k) = 0;
                if (Fn(i,j,k) > 0)
                    alphan(i,j,k) = 1;
                end

                alphas(i,j,k) = 0;
                if (Fs(i,j,k) > 0)
                    alphas(i,j,k) = 1;
                end

                alphab(i,j,k) = 0;
                if (Fb(i,j,k) > 0)
                    alphab(i,j,k) = 1;
                end

                alphaf(i,j,k) = 0;
                if (Ff(i,j,k) > 0)
                    alphaf(i,j,k) = 1;
                end


                SUDC(i,j,k) = ((1/2)*Fe(i,j,k) * ((1 - alphae(i,j,k)) * PSI(REM(i,j,k)) - alphae(i,j,k) * PSI(REP(i,j,k))) * (U(i+1,j,k)... 
                - U(i,j,k)))+ ((1/2)*Fw(i,j,k) * (alphaw(i,j,k) * PSI(RWP(i,j,k)) - (1 - alphaw(i,j,k)) * PSI(RWM(i,j,k))) * ...
                (U(i,j,k) - U(i-1,j,k)))+ ((1/2)*Fn(i,j,k) * ((1 - alphan(i,j,k)) * PSI(RNM(i,j,k)) - alphan(i,j,k) * PSI(RNP(i,j,k)))...
                * (U(i,j,k+1) - U(i,j,k)))+ ((1/2)*Fs(i,j,k) * (alphas(i,j,k) * PSI(RSP(i,j,k)) - (1 - alphas(i,j,k)) * PSI(RSM(i,j,k)))...
                * (U(i,j,k) - U(i,j,k-1))) + ((1/2)*Ff(i,j,k) * ((1 - alphaf(i,j,k)) * PSI(RFM(i,j,k)) - alphaf(i,j,k) * PSI(RFP(i,j,k)))...
                * (U(i,j+1,k) - U(i,j,k)) + ((1/2)*Fb(i,j,k) * (alphab(i,j,k) * PSI(RBP(i,j,k)) - (1 - alphab(i,j,k)) * PSI(RBM(i,j,k))) * (U(i,j,k) - U(i,j-1,k))));

                CON(i,j,k) = CON(i,j,k) * VOL + REL * AP(i,j,k) * U(i,j,k) + SUDC(i,j,k);
                DU(i,j,k) = VOL/XDIF(i);
%                 DU(i,j,k) = DU(i,j,k)/AP(i,j,k);

                CON(i,j,k) = CON(i,j,k) + DU(i,j,k) * AP(i,j,k) * (P(i-1,j,k) - P(i,j,k));

            end
        end
    end


                
                    for k = 2:N2
                        for j = 2:M2
                            for i = 3:L2
                                COFU(i,j,k,1) = CON(i,j,k);
                                COFU(i,j,k,2) = Ae(i,j,k);
                                COFU(i,j,k,3) = Aw(i,j,k);
                                COFU(i,j,k,4) = Af(i,j,k);
                                COFU(i,j,k,5) = Ab(i,j,k);
                                COFU(i,j,k,6) = As(i,j,k);
                                COFU(i,j,k,7) = An(i,j,k);
                                COFU(i,j,k,8) = AP(i,j,k);
                            end
                        end
                    end


                

            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NF = 2;

IST = 2;
JST = 3;
KST = 2;

 for k = 2:N2
    for j = 2:M2
        for i = 2:L2
            CON(i,j,k) = 0;
            AP(i,j,k) = 0;
        end
    end
 end

 [GAM,CON,AP] = Diffusion_Source(U, V, W, T, NF,ITER,LAST);


for k = 3:N3
        for j = 3:M2
            for i = 3:L3
                if abs(V(i+1,j,k)-V(i,j,k)) > 1.0E-20
                        REP(i,j,k) = (V(i,j,k) - V(i-1,j,k))/(V(i+1,j,k) - V(i,j,k));
                    if i == L3
                        REM(i,j,k) = 1;
                    else
                        REM(i,j,k) = (V(i+2,j,k) - V(i+1,j,k))/(V(i+1,j,k) - V(i,j,k));
                    end
                else
                    REP(i,j,k) = 1;
                    REM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(V(i,j,k)-V(i-1,j,k)) > 1.0E-20
                    if i ==3
                        RWP(i,j,k) = 1;
                    else
                        RWP(i,j,k) = (V(i-1,j,k) - V(i-2,j,k))/(V(i,j,k) - V(i-1,j,k));
                    end
                    RWM(i,j,k) = (V(i+1,j,k) - V(i,j,k))/(V(i,j,k) - V(i-1,j,k));
                else
                    RWP(i,j,k) = 1;
                    RWM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------


               if abs(V(i,j+1,k)-V(i,j,k)) > 1.0E-20
                    RBP(i,j,k) = (V(i,j,k) - V(i,j-1,k))/(V(i,j+1,k) - V(i,j,k));
                    if j == M2
                        RBM(i,j,k) = 1;
                    else
                        RBM(i,j,k) = (V(i,j+2,k) - V(i,j+1,k))/(V(i,j+1,k) - V(i,j,k));
                    end
                else
                    RBP(i,j,k) = 1;
                    RBM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(V(i,j,k)-V(i,j-1,k)) > 1.0E-20
                    if j ==3
                        RFP(i,j,k) = 1;
                    else
                        RFP(i,j,k) = (V(i,j-1,k) - V(i,j-2,k))/(V(i,j,k) - V(i,j-1,k));
                    end
                    RFM(i,j,k) = (V(i,j+1,k) - V(i,j,k))/(V(i,j,k) - V(i,j-1,k));
                else
                    RFP(i,j,k) = 1;
                    RFM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------

                 if abs(V(i,j,k+1)-V(i,j,k)) > 1.0E-20
                    RNP(i,j,k) = (V(i,j,k) - V(i,j,k-1))/(V(i,j,k+1) - V(i,j,k));
                    if k == N3
                        RNM(i,j,k) = 1;
                    else
                        RNM(i,j,k) = (V(i,j,k+2) - V(i,j,k+1))/(V(i,j,k+1) - V(i,j,k));
                    end
                else
                    RNP(i,j,k) = 1;
                    RNM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(V(i,j,k)-V(i,j,k-1)) > 1.0E-20
                    if k ==3
                        RSP(i,j,k) = 1;
                    else
                        RSP(i,j,k) = (V(i,j,k-1) - V(i,j,k-2))/(V(i,j,k) - V(i,j,k-1));
                    end
                    RSM(i,j,k) = (V(i,j,k+1) - V(i,j,k))/(V(i,j,k) - V(i,j,k-1));
                else
                    RSP(i,j,k) = 1;
                    RSM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------
            end
        end
    end
%-----------------------------------------------------------------
    k = 2;
    for j = 3:M2
        for i = 2:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    k = N2;
    for j = 3:M2
        for i = 2:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end
%------------------------------------------------------------
      i = 2;
    for k = 2:N2
        for j = 3:M2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    i = L2;
    for k = 2:N2
        for j = 3:M2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    for k = 2:N2
        for i = 2:L2 
             FLOW = AY(i,k) * V(i,2,k) * RHO(i,1,k);
             DIFF = AY(i,k) * GAM(i,1,k)/YCV(2);
             Af(i,3,k) = DIFF + max(0, FLOW);
             Ff(i,3,k) =  FLOW;
        end
    end 

    for j = 3:M2
        for i = 2:L2
            FL = YCVJ(j) * W(i,j,2) * RHO(i,j,1);
            FLM = YCVJP(j-1) * W(i,j-1,2) * RHO(i,j-1,1);
            FLOW = XCV(i) * (FL + FLM);
            DIFF = XCV(i) * (YCVJ(j) * GAM(i,j,1) + YCVJP(j-1) * GAM(i,j-1,1))/ZDIF(2);
            As(i,j,2) = DIFF + max(0,FLOW);
            Fs(i,j,2) = FLOW;
        end
    end

    for k = 2:N2
        for j = 3:M2
             FL = YCVJ(j) * U(2,j,k) * RHO(1,j,k);
             FLM = YCVJP(j-1) * U(2,j-1,k) * RHO(1,j-1,k);
             FLOW = (FL + FLM) * ZCV(k);
             DIFF = ZCV(k) * (YCVJ(j)*GAM(1,j,k) + YCVJP(j-1)*GAM(1,j-1,k))/XDIF(2);
             Aw(2,j,k) = DIFF + max(0,FLOW);
             Fw(2,j,k) = FLOW;

            for i = 2:L2
                if i == L2
                     FL = YCVJ(j) * U(L1,j,k) * RHO(L1,j,k);
                     FLM = YCVJP(j-1) * U(L1,j-1,k) * RHO(L1,j-1,k);
                     DIFF = ZCV(k) * (YCVJ(j)*GAM(L1,j,k) + YCVJP(j-1)*GAM(L1,j-1,k))/XDIF(L1);
                     FLOW = (FL+FLM)*ZCV(k);

                else
                     FL = YCVJ(j) * U(i+1,j,k)*(FX(i+1)*RHO(i+1,j,k) + FXM(i+1)*RHO(i,j,k));
                     FLM = YCVJP(j-1) * U(i+1,j-1,k)*(FX(i+1)*RHO(i+1,j-1,k) + FXM(i+1)*RHO(i,j-1,k));
                     FLOW = (FL + FLM) * ZCV(k);
                     GM = GAM(i,j,k) * GAM(i+1,j,k)/(XCV(i)*GAM(i+1,j,k) + XCV(i+1)*GAM(i,j,k))*YCVJP(j-1);
                     GMM = GAM(i,j-1,k) * GAM(i+1,j-1,k)/ (XCV(i) * GAM(i+1,j-1,k) + XCV(i+1)*GAM(i,j-1,k))*YCVJP(j-1);
                     DIFF = 2 .* (GM + GMM) * ZCV(k);
                end
                Aw(i+1,j,k) = DIFF + max(0,FLOW);
                Ae(i,j,k) = DIFF + max(0,-FLOW); %Different in his code
                Fw(i+1,j,k) = FLOW;
                Fe(i,j,k) = FLOW;

                if j == M2
                     FLOW = AY(i,k)*V(i,M1,k)*RHO(i,M1,k);
                     DIFF = AY(i,k)*GAM(i,M1,k)/YCV(M2);
                else
                     FL = V(i,j,k)*(FY(j) * RHO(i,j,k) + FYM(j)*RHO(i,j-1,k));
                     FLP = V(i,j+1,k)*(FY(j+1)*RHO(i,j+1,k) + FYM(j+1)*RHO(i,j,k));
                     FLOW = 0.5 * AY(i,k) * (FL + FLP);
                     DIFF = AY(i,k)*GAM(i,j,k)/YCV(j);
                end
                Af(i,j+1,k) = DIFF + max(0,FLOW);
                Ab(i,j,k) = DIFF + max(0,-FLOW);
                Ff(i,j+1,k) = FLOW;
                Fb(i,j,k) = FLOW;

                if k == N2
                    FL = XCVI(i) * W(i,j,N1) * RHO(i,j,N1);
                    FLM = XCVIP(i-1) * W(i-1,j,N1) * RHO(i-1,j,N1);
                    DIFF =  YCV(j) * (XCVI(i) * GAM(i,j,N1) + XCVIP(i-1) * GAM(i-1,j,N1))/ZDIF(N1);
                    FLOW = YCV(j) * (FL + FLM);
                else
                    FL = XCVI(i) * W(i,j,k+1) * (FZ(k+1) * RHO(i,j,k+1) + FZM(k+1) * RHO(i,j,k));
                    FLM = XCVIP(i-1) * W(i-1,j,k+1) * (FZ(k+1) * RHO(i-1,j,k+1) + FZM(k+1) * RHO(i-1,j,k));
                    FLOW = YCV(j) * (FL + FLM);
                    GM = GAM(i,j,k) * GAM(i,j,k+1)/(ZCV(k) * GAM(i,j,k+1) + ZCV(k+1) * GAM(i,j,k)) * XCV(i);
                    GMM = GAM(i-1,j,k) * GAM(i-1,j,k+1)/(ZCV(k) * GAM(i-1,j,k+1) + ZCV(k+1)*GAM(i-1,j,k)) * XCVIP(i-1);
                    DIFF = YCV(j) * 2 .* (GM + GMM);
                end
                As(i,j,k+1) = DIFF + max(0,FLOW);
                An(i,j,k) = DIFF + max(0,-FLOW);
                Fs(i,j,k+1) = FLOW;
                Fn(i,j,k) = FLOW;

                VOL = AY(i,k) * YCVS(j);
                APT = (RHO(i,j,k)*YCVJ(j) + RHO(i,j-1,k)*YCVJP(j-1))/(YCVS(j)*DT);
                AP(i,j,k) = AP(i,j,k) - APT;
                CON(i,j,k) = CON(i,j,k) + APT*U(i,j,k);
                FLOWSUM(i,j,k) = Fe(i,j,k) - Fw(i,j,k) + Fn(i,j,k) - Fs(i,j,k) + Fb(i,j,k) - Ff(i,j,k);
                AP(i,j,k) = (-AP(i,j,k) * VOL + FLOWSUM(i,j,k) + Aw(i,j,k) + Ae(i,j,k) + Af(i,j,k) + Ab(i,j,k) + As(i,j,k) + An(i,j,k))/RELAX;

                alphaw(i,j,k) = 0;
                if (Fw(i,j,k) > 0)
                    alphaw(i,j,k) = 1;
                end

                alphae(i,j,k) = 0;
                if (Fe(i,j,k) > 0)
                    alphae(i,j,k) = 1;
                end

                alphan(i,j,k) = 0;
                if (Fn(i,j,k) > 0)
                    alphan(i,j,k) = 1;
                end

                alphas(i,j,k) = 0;
                if (Fs(i,j,k) > 0)
                    alphas(i,j,k) = 1;
                end

                alphab(i,j,k) = 0;
                if (Fb(i,j,k) > 0)
                    alphab(i,j,k) = 1;
                end

                alphaf(i,j,k) = 0;
                if (Ff(i,j,k) > 0)
                    alphaf(i,j,k) = 1;
                end


                SUDC(i,j,k) = ((1/2)* Fe(i,j,k) * ((1 - alphae(i,j,k)) * PSI(REM(i,j,k)) - alphae(i,j,k) * PSI(REP(i,j,k)))...
                    * (U(i+1,j,k) - U(i,j,k)))+ ((1/2)*Fw(i,j,k) * (alphaw(i,j,k) * PSI(RWP(i,j,k)) - (1 - alphaw(i,j,k))...
                    * PSI(RWM(i,j,k))) * (U(i,j,k) - U(i-1,j,k)))+ ((1/2)*Fn(i,j,k) * ((1 - alphan(i,j,k)) * PSI(RNM(i,j,k))...
                    - alphan(i,j,k) * PSI(RNP(i,j,k))) * (U(i,j,k+1) - U(i,j,k)))+ ((1/2)*Fs(i,j,k) * (alphas(i,j,k) * PSI(RSP(i,j,k))...
                    - (1 - alphas(i,j,k)) * PSI(RSM(i,j,k))) * (U(i,j,k) - U(i,j,k-1))) + ((1/2)*Ff(i,j,k) * ((1 - alphaf(i,j,k)) * PSI(RFM(i,j,k))...
                    - alphaf(i,j,k) * PSI(RFP(i,j,k))) * (U(i,j+1,k) - U(i,j,k)) + ((1/2)*Fb(i,j,k) * (alphab(i,j,k) * PSI(RBP(i,j,k))...
                    - (1 - alphab(i,j,k)) * PSI(RBM(i,j,k))) * (U(i,j,k) - U(i,j-1,k))));

                CON(i,j,k) = CON(i,j,k)*VOL + REL * AP(i,j,k) * V(i,j,k) + SUDC(i,j,k);
                DV(i,j,k) = VOL/YDIF(j);
                DV(i,j,k) = DV(i,j,k)/AP(i,j,k);

                CON(i,j,k) = CON(i,j,k) + DV(i,j,k) * AP(i,j,k) * (P(i,j-1,k) - P(i,j,k));

            end
        end
    end

                
                    for k = 2:N2
                        for j = 3:M2
                            for i = 2:L2
                                COFV(i,j,k,1) = CON(i,j,k);
                                COFV(i,j,k,2) = Ae(i,j,k);
                                COFV(i,j,k,3) = Aw(i,j,k);
                                COFV(i,j,k,4) = Af(i,j,k);
                                COFV(i,j,k,5) = Ab(i,j,k);
                                COFV(i,j,k,6) = As(i,j,k);
                                COFV(i,j,k,7) = An(i,j,k);
                                COFV(i,j,k,8) = AP(i,j,k);
                            end
                        end
                    end


           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NF = 3;
IST = 2;
JST = 2;
KST = 3;


for k = 2:N2
    for j = 2:M2
        for i = 2:L2
            CON(i,j,k) = 0;
            AP(i,j,k) = 0;
        end
    end
end

[GAM,CON,AP] = Diffusion_Source(U, V, W, T, NF,ITER,LAST);

 for k = 3:N2
        for j = 3:M3
            for i = 3:L3
                if abs(W(i+1,j,k)-W(i,j,k)) > 1.0E-20
                        REP(i,j,k) = (W(i,j,k) - W(i-1,j,k))/(W(i+1,j,k) - W(i,j,k));
                    if i == L3
                        REM(i,j,k) = 1;
                    else
                        REM(i,j,k) = (W(i+2,j,k) - W(i+1,j,k))/(W(i+1,j,k) - W(i,j,k));
                    end
                else
                    REP(i,j,k) = 1;
                    REM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(W(i,j,k)-W(i-1,j,k)) > 1.0E-20
                    if i ==3
                        RWP(i,j,k) = 1;
                    else
                        RWP(i,j,k) = (W(i-1,j,k) - W(i-2,j,k))/(W(i,j,k) - W(i-1,j,k));
                    end
                    RWM(i,j,k) = (W(i+1,j,k) - W(i,j,k))/(W(i,j,k) - W(i-1,j,k));
                else
                    RWP(i,j,k) = 1;
                    RWM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------


               if abs(W(i,j+1,k)-W(i,j,k)) > 1.0E-20
                    RBP(i,j,k) = (W(i,j,k) - W(i,j-1,k))/(W(i,j+1,k) - W(i,j,k));
                    if j == M3
                        RBM(i,j,k) = 1;
                    else
                        RBM(i,j,k) = (W(i,j+2,k) - W(i,j+1,k))/(W(i,j+1,k) - W(i,j,k));
                    end
                else
                    RBP(i,j,k) = 1;
                    RBM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(W(i,j,k)-W(i,j-1,k)) > 1.0E-20
                    if j ==3
                        RFP(i,j,k) = 1;
                    else
                        RFP(i,j,k) = (W(i,j-1,k) - W(i,j-2,k))/(W(i,j,k) - W(i,j-1,k));
                    end
                    RFM(i,j,k) = (W(i,j+1,k) - W(i,j,k))/(W(i,j,k) - W(i,j-1,k));
                else
                    RFP(i,j,k) = 1;
                    RFM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------

                 if abs(W(i,j,k+1)-W(i,j,k)) > 1.0E-20
                    RNP(i,j,k) = (W(i,j,k) - W(i,j,k-1))/(W(i,j,k+1) - W(i,j,k));
                    if k == N2
                        RNM(i,j,k) = 1;
                    else
                        RNM(i,j,k) = (W(i,j,k+2) - W(i,j,k+1))/(W(i,j,k+1) - W(i,j,k));
                    end
                else
                    RNP(i,j,k) = 1;
                    RNM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(W(i,j,k)-W(i,j,k-1)) > 1.0E-20
                    if k ==3
                        RSP(i,j,k) = 1;
                    else
                        RSP(i,j,k) = (W(i,j,k-1) - W(i,j,k-2))/(W(i,j,k) - W(i,j,k-1));
                    end
                    RSM(i,j,k) = (W(i,j,k+1) - W(i,j,k))/(W(i,j,k) - W(i,j,k-1));
                else
                    RSP(i,j,k) = 1;
                    RSM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------
            end
        end
    end
%-----------------------------------------------------------------63
    i = 2;
    for j = 2:M2
        for k = 3:N2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    i = L2;
    for j = 2:M2
        for k = 3:N2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end
%------------------------------------------------------------
      j = 2;
    for k = 3:N2
        for i = 2:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    j = M2;
    for k = 3:N2
        for i = 2:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end


    for k = 3:N2
        for i = 2:L2
            FL = ZCVK(k) * V(i,2,k) * RHO(i,1,k);
            FLM = ZCVKP(k-1) * V(i,2,k-1) * RHO(i,1,k-1);
            FLOW = XCV(i) * (FL + FLM);
            DIFF = XCV(i) * (ZCVK(k) * GAM(i,1,k) + ZCVKP(k-1) * GAM(i,1,k-1))/YDIF(2);
            Af(i,2,k) = DIFF + max(0,FLOW);
            Ff(i,2,k) = FLOW;
        end
    end

    for j = 2:M2
        for i = 2:L2
            FLOW = AZ(i,j) * W(i,j,2) * RHO(i,j,1);
            DIFF = AZ(i,j) * GAM(i,j,1) /ZCV(2);
            As(i,j,3) = DIFF + max(0, FLOW);
            Fs(i,j,3) = FLOW;
        end 
    end

    for k = 3:N2
        for j = 2:M2
            FL = ZCVK(k) * U(2,j,k) * RHO(1,j,k);
            FLM = ZCVKP(k-1) * U(2,j,k-1) * RHO(1,j,k-1);
            FLOW = (FL + FLM) * YCV(j);
            DIFF = YCV(j) * (ZCVK(k) * GAM(1,j,k) + ZCVKP(k-1) * GAM(1,j,k-1))/XDIF(2);
            Aw(2,j,k) = DIFF + max(0,FLOW);
            Fw(2,j,k) = FLOW;

            for i = 2:L2
                if (i==L2)
                    FL = ZCVK(k) * U(L1,j,k) * RHO(L1,j,k);
                    FLM = ZCVKP(k-1) * U(L1,j,k-1) * RHO(L1,j,k-1);
                    DIFF = YCV(j) * (ZCVK(k) * GAM(L1,j,k) + ZCVKP(k-1) * GAM(L1,j,k-1))/XDIF(L1);
                    FLOW = (FL + FLM) * YCV(j);

                else
                    FL = ZCVK(k) * U(i+1,j,k) * (FX(i+1) * RHO(i+1,j,k) + FXM(i+1) * RHO(i,j,k));
                    FLM = ZCVKP(k-1) * U(i+1,j,k-1) * (FX(i+1) * RHO(i+1,j,k-1) + FXM(i+1) * RHO(i,j,k-1));
                    FLOW = (FL + FLM) * YCV(j);
                    GM = GAM(i,j,k) * GAM(i+1,j,k) / (XCV(i) * GAM(i+1,j,k) + XCV(i+1) * GAM(i,j,k)) * ZCVK(k);
                    GMM = GAM(i,j,k-1) * GAM(i+1,j,k-1) / (XCV(i) * GAM(i+1,j,k-1) + XCV(i+1) * GAM(i,j,k-1)) * ZCVKP(k-1);
                    DIFF = 2 .* (GM + GMM) * YCV(j);
                end

                Aw(i+1,j,k) = DIFF + max(0,FLOW);
                Ae(i,j,k) = DIFF + max(0,-FLOW);
                Fw(i+1,j,k) = FLOW;
                Fe(i,j,k) = FLOW;

                if (k == N2)
                    FLOW = AZ(i,j) * W(i,j,N1) * RHO(i,j,N1);
                    DIFF = AZ(i,j) * GAM(i,j,N1) / ZCV(N2);

                else
                    FL = W(i,j,k) * (FZ(k) * RHO(i,j,k) + FZM(k) * RHO(i,j,k-1));
                    FLP = W(i,j,k+1) * (FZ(k+1) * RHO(i,j,k+1) + FZM(k+1) * RHO(i,j,k));
                    FLOW = 0.5 * AZ(i,j) * (FL + FLP);
                    DIFF = AZ(i,j) * GAM(i,j,k) / ZCV(k);
                end

                As(i,j,k+1) = DIFF + max(0,FLOW);
                An(i,j,k) = DIFF + max(0,-FLOW);
                Fs(i,j,k+1) = FLOW;
                Fn(i,j,k) = FLOW;

                if (j == M2)
                    FL = ZCVK(k) * V(i,M1,k) * RHO(i,M1,k);
                    FLM = ZCVKP(k-1) * V(i,M1,k-1) * RHO(i,M1,k-1);
                    DIFF = XCV(i) * (ZCVK(k) * GAM(i,M1,k) + ZCVKP(k-1) * GAM(i,M1,k-1)) / YDIF(M1);
                    FLOW = XCV(i) * (FL + FLM);

                else
                    FL = ZCVK(k) * V(i,j+1,k) * (FY(j+1) * RHO(i,j+1,k) + FYM(j+1) * RHO(i,j,k));
                    FLM = ZCVKP(k-1) * V(i,j+1,k-1) * (FY(j+1) * RHO(i,j+1,k-1) + FYM(j+1) * RHO(i,j,k-1));
                    FLOW = XCV(i) * (FL + FLM);
                    GM = GAM(i,j,k) * GAM(i,j+1,k) / (YCV(j) * GAM(i,j+1,k) + YCV(j+1) * GAM(i,j,k)) * ZCVK(k);
                    GMM = GAM(i,j,k-1) * GAM(i,j+1,k-1) /(YCV(j) * GAM(i,j+1,k-1) + YCV(j+1) * GAM(i,j,k-1)) * ZCVKP(k-1);
                    DIFF = XCV(i) * 2 .* (GM + GMM);

                end

                Af(i,j+1,k) = DIFF + max(0,FLOW);
                Ab(i,j,k) = DIFF + max(0,-FLOW);
                Ff(i,j+1,k) = FLOW;
                Fb(i,j,k) = FLOW;

                VOL = AZ(i,j) * ZCVS(k);
                APT = (RHO(i,j,k)*ZCVK(k) + RHO(i,j,k-1)*ZCVKP(k-1))/(ZCVS(k)*DT);
                AP(i,j,k) = AP(i,j,k) - APT;
                CON(i,j,k) = CON(i,j,k) + APT*W(i,j,k);
                FLOWSUM(i,j,k) = Fe(i,j,k) - Fw(i,j,k) + Fn(i,j,k) - Fs(i,j,k) + Fb(i,j,k) - Ff(i,j,k);
                AP(i,j,k) = (-AP(i,j,k) * VOL + FLOWSUM(i,j,k) + Aw(i,j,k) + Ae(i,j,k) + Af(i,j,k) + Ab(i,j,k) + As(i,j,k) + An(i,j,k))/RELAX;

                alphaw(i,j,k) = 0;
                if (Fw(i,j,k) > 0)
                    alphaw(i,j,k) = 1;
                end

                alphae(i,j,k) = 0;
                if (Fe(i,j,k) > 0)
                    alphae(i,j,k) = 1;
                end

                alphan(i,j,k) = 0;
                if (Fn(i,j,k) > 0)
                    alphan(i,j,k) = 1;
                end

                alphas(i,j,k) = 0;
                if (Fs(i,j,k) > 0)
                    alphas(i,j,k) = 1;
                end

                alphab(i,j,k) = 0;
                if (Fb(i,j,k) > 0)
                    alphab(i,j,k) = 1;
                end

                alphaf(i,j,k) = 0;
                if (Ff(i,j,k) > 0)
                    alphaf(i,j,k) = 1;
                end


                SUDC(i,j,k) = ((1/2)*Fe(i,j,k) * ((1 - alphae(i,j,k)) * PSI(REM(i,j,k)) - alphae(i,j,k) * PSI(REP(i,j,k))) * (U(i+1,j,k) - U(i,j,k)))+ ((1/2)*Fw(i,j,k) * (alphaw(i,j,k) * PSI(RWP(i,j,k))- (1 - alphaw(i,j,k)) * PSI(RWM(i,j,k))) * (U(i,j,k) - U(i-1,j,k)))+ ((1/2)*Fn(i,j,k) * ((1 - alphan(i,j,k)) * PSI(RNM(i,j,k)) - alphan(i,j,k) * PSI(RNP(i,j,k))) * (U(i,j,k+1) - U(i,j,k)))+ ((1/2)*Fs(i,j,k) * (alphas(i,j,k) * PSI(RSP(i,j,k)) - (1 - alphas(i,j,k)) * PSI(RSM(i,j,k))) * (U(i,j,k) - U(i,j,k-1))) + ((1/2)*Ff(i,j,k) * ((1 - alphaf(i,j,k)) * PSI(RFM(i,j,k)) - alphaf(i,j,k) * PSI(RFP(i,j,k))) * (U(i,j+1,k) - U(i,j,k)) + ((1/2)*Fb(i,j,k) * (alphab(i,j,k) * PSI(RBP(i,j,k)) - (1 - alphab(i,j,k)) * PSI(RBM(i,j,k))) * (U(i,j,k) - U(i,j-1,k))));
                CON(i,j,k) = CON(i,j,k)*VOL + REL * AP(i,j,k) * W(i,j,k) + SUDC(i,j,k);
                DW(i,j,k) = VOL/ZDIF(k);
                DW(i,j,k) = DW(i,j,k)/AP(i,j,k);

                CON(i,j,k) = CON(i,j,k) + DW(i,j,k) * AP(i,j,k) * (P(i,j,k-1) - P(i,j,k));

            end
        end
    end

                
                    for k = 3:N2
                        for j = 2:M2
                            for i = 2:L2
                                COFW(i,j,k,1) = CON(i,j,k);
                                COFW(i,j,k,2) = Ae(i,j,k);
                                COFW(i,j,k,3) = Aw(i,j,k);
                                COFW(i,j,k,4) = Af(i,j,k);
                                COFW(i,j,k,5) = Ab(i,j,k);
                                COFW(i,j,k,6) = As(i,j,k);
                                COFW(i,j,k,7) = An(i,j,k);
                                COFW(i,j,k,8) = AP(i,j,k);
                            end
                        end
                    end

            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for k = 2:N2
        for j = 2:M2
            for i = 3:L2
                UHAT(i,j,k) = (COFU(i,j,k,2) * U(i+1,j,k) + COFU(i,j,k,3) * U(i-1,j,k) + COFU(i,j,k,4) * U(i,j+1,k) + COFU(i,j,k,5) * U(i,j-1,k) + COFU(i,j,k,6) * U(i,j,k+1) + COFU(i,j,k,7) * U(i,j,k-1) + COFU(i,j,k,1))/ (COFU(i,j,k,8));
            end
        end
    end

    for k = 2:N2
        for j = 3:M2
            for i = 2:L2
                VHAT(i,j,k) = (COFV(i,j,k,2) * V(i+1,j,k) + COFV(i,j,k,3) * V(i-1,j,k) + COFV(i,j,k,4) * V(i,j+1,k) + COFV(i,j,k,5) * V(i,j-1,k) + COFV(i,j,k,6) * V(i,j,k+1) + COFV(i,j,k,7) * V(i,j,k-1) + COFV(i,j,k,1))/ (COFV(i,j,k,8));
            end
        end
    end

    for k = 3:N2
        for j = 2:M2
            for i = 2:L2
                WHAT(i,j,k) = (COFW(i,j,k,2) * W(i+1,j,k) + COFW(i,j,k,3) * W(i-1,j,k) + COFW(i,j,k,4) * W(i,j+1,k) + COFW(i,j,k,5) * W(i,j-1,k) + COFW(i,j,k,6) * W(i,j,k+1) + COFW(i,j,k,7) * W(i,j,k-1) + COFW(i,j,k,1))/ (COFW(i,j,k,8));
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%% coefficients for Pressure equation %%%%%%%%%%%%

    NF = 4;

    IST = 2;
    JST = 2;
    KST = 2;


for k = 2:N2
    for j = 2:M2
        for i = 2:L2
            CON(i,j,k) = 0;
            AP(i,j,k) = 0;
        end
    end
end

[GAM,CON,AP] = Diffusion_Source(U, V, W, T, NF,ITER,LAST) ;

for k = 2:N2
    for j = 2:M2
        for i = 2:L2
            VOL = YCV(j) * XCV(i) * ZCV(k);
            CON(i,j,k) = CON(i,j,k) * VOL;
        end
    end
end

for i = 2:L2
    for k = 2:N2
        ARHO = AY(i,k) * RHO(i,1,k);
        CON(i,2,k) = CON(i,2,k) + ARHO * V(i,2,k);
        Ab(i,2,k) = 0;
    end

    for j = 2:M2
        ARHO = AZ(i,j) * RHO(i,j,1);
        CON(i,j,2) = CON(i,j,2) + ARHO * W(i,j,2);
        As(i,j,2) = 0;
    end
end

for k = 2:N2
    for j = 2:M2
        ARHO = AX(j,k) * RHO(1,j,k);
        CON(2,j,k) = CON(2,j,k) + ARHO * U(2,j,k);
        Aw(2,j,k) = 0;
        for i = 2:L2

            if i == L2
                ARHO = AX(j,k) * RHO(L1,j,k);
                CON(L2,j,k) = CON(L2,j,k) - ARHO * U(L1,j,k);
                Ae(i,j,k) = 0;
            else
                ARHO = AX(j,k) * (FX(i+1) * RHO(i+1,j,k) + FXM(i+1) * RHO(i,j,k));
                FLOW = ARHO * UHAT(i+1,j,k);
                CON(i,j,k) = CON(i,j,k) - FLOW;
                CON(i+1,j,k) = CON(i+1,j,k) + FLOW;
                Ae(i,j,k) = ARHO * DU(i+1,j,k);
                Aw(i+1,j,k) = Ae(i,j,k);
            end

            if j == M2
                ARHO = AY(i,k) * RHO(i,M1,k);
                CON(i,M2,k) = CON(i,M2,k) - ARHO * V(i,M1,k);
                Af(i,j,k) = 0;
            else
                ARHO = AY(i,k) * (FY(j+1) * RHO(i,j+1,k) + FYM(j+1) * RHO(i,j,k));
                FLOW = ARHO * VHAT(i,j+1,k);
                CON(i,j,k) = CON(i,j,k) - FLOW;
                CON(i,j+1,k) = CON(i,j+1,k) + FLOW;
                Af(i,j,k) = ARHO * DV(i,j+1,k);
                Ab(i,j+1,k) = Af(i,j,k);
            end

             if k == N2
                ARHO = AZ(i,j) * RHO(i,j,N1);
                CON(i,j,N2) = CON(i,j,N2) - ARHO * W(i,j,N1);
                An(i,j,k) = 0;
            else
                ARHO = AZ(i,j) * (FZ(k+1) * RHO(i,j,k+1) + FZM(k+1) * RHO(i,j,k));
                FLOW = ARHO * WHAT(i,j,k+1);
                CON(i,j,k) = CON(i,j,k) - FLOW;
                CON(i,j,k+1) = CON(i,j,k+1) + FLOW;
                An(i,j,k) = ARHO * DW(i,j,k+1);
                As(i,j,k+1) = An(i,j,k);
             end

             AP(i,j,k) = Ae(i,j,k) + Aw(i,j,k) + Ab(i,j,k) + Af(i,j,k) + An(i,j,k) + As(i,j,k);
             
        end
    end
end

for k = 2:N2
    for j = 2:M2
        for i = 2:L2
            AP(i,j,k) = AP(i,j,k) / RELAX;
            CON(i,j,k) = CON(i,j,k) + (1 - RELAX)*AP(i,j,k) * P(i,j,k);
        end
    end
end

for i = 2:L2
    for j = 2:M2
        for k = 2:N2
            COFP(i,j,k,2) = Ae(i,j,k);
            COFP(i,j,k,3) = Aw(i,j,k);
            COFP(i,j,k,4) = Af(i,j,k);
            COFP(i,j,k,5) = Ab(i,j,k);
            COFP(i,j,k,6) = As(i,j,k);
            COFP(i,j,k,7) = An(i,j,k);
        end
    end
end

[P] = SOLVE(P, IST, JST, KST, AP, Af, Ab, Ae, Aw, An, As, CON, ITER);

%%%%%%%% Calculate Velocity U %%%%%%%

IST = 3;
JST = 2;
KST = 2;

for i = 3:L2
    for j = 2:M2
        for k = 2:N2
            CON(i,j,k) = COFU(i,j,k,1);
            Ae(i,j,k) = COFU(i,j,k,2);
            Aw(i,j,k) = COFU(i,j,k,3);
            Af(i,j,k) = COFU(i,j,k,4);
            Ab(i,j,k) = COFU(i,j,k,5);
            As(i,j,k) = COFU(i,j,k,6);
            An(i,j,k) = COFU(i,j,k,7);
            AP(i,j,k) = COFU(i,j,k,8);
        end
    end
end


for i = IST:L2
    for j = JST:M2
        for k = KST:N2
            CON(i,j,k) = CON(i,j,k) + DU(i,j,k) * AP(i,j,k) * (P(i-1,j,k) - P(i,j,k));
            COFU(i,j,k,1) = CON(i,j,k);
        end
    end
end




[U] = SOLVE(U, IST, JST, KST, AP, Af, Ab, Ae, Aw, An, As, CON, ITER); %Call Solve

%%%%%%%% Calculate Velocity V %%%%%%%

IST = 2;
JST = 3;
KST = 2;


for i = 2:L2
    for j = 3:M2
        for k = 2:N2
            CON(i,j,k) = COFV(i,j,k,1);
            Ae(i,j,k) = COFV(i,j,k,2);
            Aw(i,j,k) = COFV(i,j,k,3);
            Af(i,j,k) = COFV(i,j,k,4);
            Ab(i,j,k) = COFV(i,j,k,5);
            As(i,j,k) = COFV(i,j,k,6);
            An(i,j,k) = COFV(i,j,k,7);
            AP(i,j,k) = COFV(i,j,k,8);
        end
    end
end


for k = KST:N2
    for j = JST:M2
        for i = IST:L2
            CON(i,j,k) = CON(i,j,k) + DV(i,j,k) * AP(i,j,k) * (P(i,j-1,k) - P(i,j,k));
            COFV(i,j,k,1) = CON(i,j,k);
        end
    end
end



[V] = SOLVE(V, IST, JST, KST, AP, Af, Ab, Ae, Aw, An, As, CON, ITER);%%%%%%%%%CALL SOLVE

%%%%%%%% Calculate Velocity W %%%%%%%

IST = 2;
JST = 2;
KST = 3;

for i = 2:L2
    for j = 2:M2
        for k = 3:N2
            CON(i,j,k) = COFW(i,j,k,1);
            Ae(i,j,k) = COFW(i,j,k,2);
            Aw(i,j,k) = COFW(i,j,k,3);
            Af(i,j,k) = COFW(i,j,k,4);
            Ab(i,j,k) = COFW(i,j,k,5);
            As(i,j,k) = COFW(i,j,k,6);
            An(i,j,k) = COFW(i,j,k,7);
            AP(i,j,k) = COFW(i,j,k,8);
        end
    end
end



for k = KST:N2
    for j = JST:M2
        for i = IST:L2
            CON(i,j,k) = CON(i,j,k) + DW(i,j,k) * AP(i,j,k) * (P(i,j,k-1) - P(i,j,k));
            COFW(i,j,k,1) = CON(i,j,k);
        end
    end
end

 

[W] = SOLVE(W, IST, JST, KST, AP, Af, Ab, Ae, Aw, An, As, CON, ITER);%%%%%%%%%CALL SOLVE

% Coefficients for the pressure correction equation------------


NF = 4;
IST = 2;
JST = 2;
KST = 2;



for k = 2:N2
    for j = 2:M2
        for i = 2:L2
            CON(i,j,k) = 0; %%%%%%%%%%%Reset
            AP(i,j,k) = 0;
        end
    end
end

for i = 2:L2
    for j = 2:M2
        for k = 2:N2
             
            Ae(i,j,k) = COFP(i,j,k,2);
            Aw(i,j,k) = COFP(i,j,k,3);
            Af(i,j,k) = COFP(i,j,k,4);
            Ab(i,j,k) = COFP(i,j,k,5);
            As(i,j,k) = COFP(i,j,k,6);
            An(i,j,k) = COFP(i,j,k,7);
            
        end
    end
end

[GAM,CON,AP] = Diffusion_Source(U, V, W, T, NF,ITER,LAST); %%%%%CALL GAMSOR

for k = 2:N2
    for j = 2:M2
        for i = 2:L2
            VOL = YCV(j) * XCV(i) * ZCV(k);
            CON(i,j,k) = CON(i,j,k) * VOL;
        end
    end
end

for i = 2:L2
    for k = 2:N2
        ARHO = AY(i,k) * RHO(i,1,k);
        CON(i,2,k) = CON(i,2,k) + ARHO * V(i,2,k);
    end

    for j = 2:M2
        ARHO = AZ(i,j) * RHO(i,j,1);
        CON(i,j,2) = CON(i,j,2) + ARHO * W(i,j,2);
    end
end

for k = 2:N2
    for j = 2:M2
        ARHO = AX(j,k) * RHO(1,j,k);
        CON(2,j,k) = CON(2,j,k) + ARHO * U(2,j,k);
        for i = 2:L2

            if i == L2
                ARHO = AX(j,k) * RHO(L1,j,k);
                CON(L2,j,k) = CON(L2,j,k) - ARHO * U(L1,j,k);
            else
                ARHO = AX(j,k) * (FX(i+1) * RHO(i+1,j,k) + FXM(i+1) * RHO(i,j,k));
                FLOW = ARHO * U(i+1,j,k);
                CON(i,j,k) = CON(i,j,k) - FLOW;
                CON(i+1,j,k) = CON(i+1,j,k) + FLOW;
            end

            if j == M2
                ARHO = AY(i,k) * RHO(i,M1,k);
                CON(i,M2,k) = CON(i,M2,k) - ARHO * V(i,M1,k);
            else
                ARHO = AY(i,k) * (FY(j+1) * RHO(i,j+1,k) + FYM(j+1) * RHO(i,j,k));
                FLOW = ARHO * V(i,j+1,k);
                CON(i,j,k) = CON(i,j,k) - FLOW;
                CON(i,j+1,k) = CON(i,j+1,k) + FLOW;
            end

            if k == N2
                ARHO = AZ(i,j) * RHO(i,j,N1);
                CON(i,j,N2) = CON(i,j,N2) - ARHO * W(i,j,N1);
            else
                ARHO = AZ(i,j) * (FZ(k+1) * RHO(i,j,k+1) + FZM(k+1) * RHO(i,j,k));
                FLOW = ARHO * W(i,j,k+1);
                CON(i,j,k) = CON(i,j,k) - FLOW;
                CON(i,j,k+1) = CON(i,j,k+1) + FLOW;
            end

            AP(i,j,k) = Ae(i,j,k) + Aw(i,j,k) + Af(i,j,k) + Ab(i,j,k) + An(i,j,k) + As(i,j,k);
            PC(i,j,k) = 0;

        end
    end
end

[PC] = SOLVE(PC, IST, JST, KST, AP, Af, Ab, Ae, Aw, An, As, CON, ITER); %SOLVE-----

%-------------Correct the Velocities------------------------

for k = 2:N2
    for j = 2:M2
        for i = 2:L2
            if i ~= 2
                U(i,j,k) = U(i,j,k) + DU(i,j,k) * (PC(i-1,j,k) - PC(i,j,k));
            end

            if j ~= 2
                V(i,j,k) = V(i,j,k) + DV(i,j,k) * (PC(i,j-1,k) - PC(i,j,k));
            end

            if k ~= 2
                W(i,j,k) = W(i,j,k) + DW(i,j,k) * (PC(i,j,k-1) - PC(i,j,k));
            end

        end
    end
end

% Residuals for U,V,W-----------------


%-----U---------------------------------

NF = 1;
IST = 3;
JST = 2;
KST = 2;

Res = 0;
R_sum = 0;

for i = 3:L2
    for j = 2:M2
        for k = 2:N2
            CON(i,j,k) = COFU(i,j,k,1);
            Ae(i,j,k) = COFU(i,j,k,2);
            Aw(i,j,k) = COFU(i,j,k,3);
            Af(i,j,k) = COFU(i,j,k,4);
            Ab(i,j,k) = COFU(i,j,k,5);
            As(i,j,k) = COFU(i,j,k,6);
            An(i,j,k) = COFU(i,j,k,7);
            AP(i,j,k) = COFU(i,j,k,8);
        end
    end
end

for i = IST:L2
    for j = JST:M2
        for k = KST:N2

            Res = (AP(i,j,k)*U(i,j,k) - (Ae(i,j,k) * U(i+1,j,k) + Aw(i,j,k) * U(i-1,j,k)...
                + Ab(i,j,k) * U(i,j+1,k) + Af(i,j,k) * U(i,j-1,k)...
                + An(i,j,k) * U(i,j,k+1) + As(i,j,k) * U(i,j,k-1) + CON(i,j,k)) ) ^ 2; 

            R_sum = R_sum + Res;
        end
    end
end

Res_U = sqrt((R_sum/((L2-IST+1) * (M2 - JST + 1) * (N2 - KST + 1))));

%----------------------V-------------------------------

NF = 2;
IST = 2;
JST = 3;
KST = 2;

Res = 0;
R_sum = 0;

for i = 3:L2
    for j = 2:M2
        for k = 2:N2
            CON(i,j,k) = COFV(i,j,k,1);
            Ae(i,j,k) = COFV(i,j,k,2);
            Aw(i,j,k) = COFV(i,j,k,3);
            Af(i,j,k) = COFV(i,j,k,4);
            Ab(i,j,k) = COFV(i,j,k,5);
            As(i,j,k) = COFV(i,j,k,6);
            An(i,j,k) = COFV(i,j,k,7);
            AP(i,j,k) = COFV(i,j,k,8);
        end
    end
end

for k = KST:N2
    for j = JST:M2
        for i = IST:L2

            Res = (AP(i,j,k)*V(i,j,k) - (Ae(i,j,k) * V(i+1,j,k) + Aw(i,j,k) * V(i-1,j,k)...
                + Ab(i,j,k) * V(i,j+1,k) + Af(i,j,k) * V(i,j-1,k)...
                + An(i,j,k) * V(i,j,k+1) + As(i,j,k) * V(i,j,k-1) + CON(i,j,k)) ) ^ 2; 

            R_sum = R_sum + Res;
        end
    end
end

Res_V = sqrt((R_sum/((L2-IST+1) * (M2 - JST + 1) * (N2 - KST + 1))));

%----------------------W-------------------------------

NF = 3;
IST = 2;
JST = 2;
KST = 3;

Res = 0;
R_sum = 0;

for i = 3:L2
    for j = 2:M2
        for k = 2:N2
            CON(i,j,k) = COFW(i,j,k,1);
            Ae(i,j,k) = COFW(i,j,k,2);
            Aw(i,j,k) = COFW(i,j,k,3);
            Af(i,j,k) = COFW(i,j,k,4);
            Ab(i,j,k) = COFW(i,j,k,5);
            As(i,j,k) = COFW(i,j,k,6);
            An(i,j,k) = COFW(i,j,k,7);
            AP(i,j,k) = COFW(i,j,k,8);
        end
    end
end

for k = KST:N2
    for j = JST:M2
        for i = IST:L2

            Res = (AP(i,j,k)*W(i,j,k) - (Ae(i,j,k) * W(i+1,j,k) + Aw(i,j,k) * W(i-1,j,k)...
                + Ab(i,j,k) * W(i,j+1,k) + Af(i,j,k) * W(i,j-1,k)...
                + An(i,j,k) * W(i,j,k+1) + As(i,j,k) * W(i,j,k-1) + CON(i,j,k)) ) ^ 2; 

            R_sum = R_sum + Res;
        end
    end
end

Res_W = sqrt((R_sum/((L2-IST+1) * (M2 - JST + 1) * (N2 - KST + 1))));

end



%------------Coefficients of Temperature----------------

NF = 5;
IST = 2;
JST = 2;
KST = 2;




for k = 2:N2
    for j = 2:M2
        for i = 2:L2
            CON(i,j,k) = 0; %RESET
            AP(i,j,k) = 0;
        end
    end
end

[GAM,CON,AP] = Diffusion_Source(U, V, W, T, NF,ITER,LAST); %CALL GAMSOR

for k = 1:N1
    for i = 1:L1
        for j = 1:M1
            RHO(i,j,k) = RHO(i,j,k) * CP(i,j,k);
        end
    end
end
 

for k = 3:N3
        for j = 3:M3
            for i = 3:L3
                if abs(T(i+1,j,k)-T(i,j,k)) > 1.0E-20
                        REP(i,j,k) = (T(i,j,k) - T(i-1,j,k))/(T(i+1,j,k) - T(i,j,k));
                    if i == L3
                        REM(i,j,k) = 1;
                    else
                        REM(i,j,k) = (T(i+2,j,k) - T(i+1,j,k))/(T(i+1,j,k) - T(i,j,k));
                    end
                else
                    REP(i,j,k) = 1;
                    REM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(T(i,j,k)-T(i-1,j,k)) > 1.0E-20
                    if i ==3
                        RWP(i,j,k) = 1;
                    else
                        RWP(i,j,k) = (T(i-1,j,k) - T(i-2,j,k))/(T(i,j,k) - T(i-1,j,k));
                    end
                    RWM(i,j,k) = (T(i+1,j,k) - T(i,j,k))/(T(i,j,k) - T(i-1,j,k));
                else
                    RWP(i,j,k) = 1;
                    RWM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------


               if abs(T(i,j+1,k)-T(i,j,k)) > 1.0E-20
                    RBP(i,j,k) = (T(i,j,k) - T(i,j-1,k))/(T(i,j+1,k) - T(i,j,k));
                    if j == M2
                        RBM(i,j,k) = 1;
                    else
                        RBM(i,j,k) = (T(i,j+2,k) - T(i,j+1,k))/(T(i,j+1,k) - T(i,j,k));
                    end
                else
                    RBP(i,j,k) = 1;
                    RBM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(T(i,j,k)-T(i,j-1,k)) > 1.0E-20
                    if j ==3
                        RFP(i,j,k) = 1;
                    else
                        RFP(i,j,k) = (T(i,j-1,k) - T(i,j-2,k))/(T(i,j,k) - T(i,j-1,k));
                    end
                    RFM(i,j,k) = (T(i,j+1,k) - T(i,j,k))/(T(i,j,k) - T(i,j-1,k));
                else
                    RFP(i,j,k) = 1;
                    RFM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------

                 if abs(T(i,j,k+1)-T(i,j,k)) > 1.0E-20
                    RNP(i,j,k) = (T(i,j,k) - T(i,j,k-1))/(T(i,j,k+1) - T(i,j,k));
                    if k == N3
                        RNM(i,j,k) = 1;
                    else
                        RNM(i,j,k) = (T(i,j,k+2) - T(i,j,k+1))/(T(i,j,k+1) - T(i,j,k));
                    end
                else
                    RNP(i,j,k) = 1;
                    RNM(i,j,k) = 1;
                end
%-------------------------------------------------------
                if abs(T(i,j,k)-T(i,j,k-1)) > 1.0E-20
                    if k ==3
                        RSP(i,j,k) = 1;
                    else
                        RSP(i,j,k) = (T(i,j,k-1) - T(i,j,k-2))/(T(i,j,k) - T(i,j,k-1));
                    end
                    RSM(i,j,k) = (T(i,j,k+1) - T(i,j,k))/(T(i,j,k) - T(i,j,k-1));
                else
                    RSP(i,j,k) = 1;
                    RSM(i,j,k) = 1;
                end
%-----------------------------------------------------------------
%----------------------------------------------------------------
            end
        end
    end
%-----------------------------------------------------------------
    i = 2;
    for j = 2:M2
        for k = 2:N2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    i = L2;
    for j = 2:M2
        for k = 2:N2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end
%------------------------------------------------------------
      j = 2;
    for k = 2:N2
        for i = 2:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    j = M2;
    for k = 2:N2
        for i = 2:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    k = 2;
    for j = 2:M2
        for i = 2:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end

    k = N2;
    for j = 2:M2
        for i = 2:L2
            REP(i,j,k) = 1;
            REM(i,j,k) = 1;
            RWP(i,j,k) = 1;
            RWM(i,j,k) = 1;

            RBP(i,j,k) = 1;
            RBM(i,j,k) = 1;
            RFP(i,j,k) = 1;
            RFM(i,j,k) = 1;

            RNP(i,j,k) = 1;
            RNM(i,j,k) = 1;
            RSP(i,j,k) = 1;
            RSM(i,j,k) = 1;
        end
    end


    for i = 2:L2
        for k = 2:N2
            FLOW = AY(i,k) * V(i,2,k) * RHO(i,1,k);
            DIFF = AY(i,k) * GAM(i,1,k)/YDIF(2);
            Af(i,2,k) = DIFF + max(0,FLOW);
            Ff(i,2,k) = FLOW;
        end

        for j = 2:M2
            FLOW = AZ(i,j) * W(i,j,2) * RHO(i,j,1);
            DIFF = AZ(i,j) * GAM(i,j,1)/ZDIF(2);
            As(i,j,2) = DIFF + max(0,FLOW);
            Fs(i,j,2) = FLOW;
        end
    end

    for k = 2:N2
        for j = 2:M2
            FLOW = AX(j,k) * RHO(1,j,k) * U(2,j,k);
            DIFF = AX(j,k) * GAM(1,j,k)/XDIF(2);
            Aw(2,j,k) = DIFF + max(0,FLOW);
            Fw(2,j,k) = FLOW;

            for i = 2:L2
                if i == L2
                    FLOW = AX(j,k) * RHO(L1,j,k) * U(L1,j,k);
                    DIFF = AX(j,k) * GAM(L1,j,k)/XDIF(L1);
                else
                    FLOW = AX(j,k) * U(i+1,j,k) * (FX(i+1) * RHO(i+1,j,k) + FXM(i+1) * RHO(i,j,k));
                    DIFF = AX(j,k) * 2 .* GAM(i,j,k) * GAM(i+1,j,k) / (XCV(i) * GAM(i+1,j,k) + XCV(i+1) * GAM(i,j,k));
                end

                Aw(i+1,j,k) = DIFF + max(0,FLOW);
                Ae(i,j,k) = DIFF + max(0,-FLOW);
                Fw(i+1,j,k) = FLOW;
                Fe(i,j,k) = FLOW;

                if j == M2

                    FLOW = AY(i,k) * RHO(i,M1,k) * V(i,M1,k);
                    DIFF = AY(i,k) * GAM(i,M1,k)/YDIF(M1);
                else
                    FLOW = AY(i,k) * V(i,j+1,k) * (FY(j+1) * RHO(i,j+1,k) + FYM(j+1) * RHO(i,j,k));
                    DIFF = AY(i,k) * 2 .* GAM(i,j,k) * GAM(i,j+1,k)/(YCV(j) * GAM(i,j+1,k) + YCV(j+1) * GAM(i,j,k));
                end

                Af(i,j+1,k) = DIFF + max(0,FLOW);
                Ab(i,j,k) = DIFF + max(0,-FLOW);
                Ff(i,j+1,k) = FLOW;
                Fb(i,j,k) = FLOW;

                if k == N2
                    FLOW = AZ(i,j) * RHO(i,j,N1) * W(i,j,N1);
                    DIFF = AZ(i,j) * GAM(i,j,N1)/ZDIF(N1);

                else
                    FLOW = AZ(i,j) * W(i,j,k+1) * (FZ(k+1) * RHO(i,j,k+1) + FZM(k+1) * RHO(i,j,k));
                    DIFF = AZ(i,j) * 2 .* GAM(i,j,k) * GAM(i,j,k+1)/(ZCV(k) * GAM(i,j,k+1) + ZCV(k+1) * GAM(i,j,k));
                end

                As(i,j,k+1) = DIFF + max(0,FLOW);
                An(i,j,k) = DIFF + max(0,-FLOW);
                Fs(i,j,k+1) = FLOW;
                Fn(i,j,k) = FLOW;

                VOL = XCV(i) * YCV(j) * ZCV(k);
                APT = RHO(i,j,k)/DT;
                AP(i,j,k) = AP(i,j,k) - APT;
                CON(i,j,k) = CON(i,j,k) + APT * T(i,j,k);

                FLOWSUM(i,j,k) = Fe(i,j,k) - Fw(i,j,k) + Fn(i,j,k) - Fs(i,j,k) + Fb(i,j,k) - Ff(i,j,k);

                AP(i,j,k) = (-AP(i,j,k)*VOL + FLOWSUM(i,j,k) + Ae(i,j,k) + Aw(i,j,k) + An(i,j,k) + As(i,j,k) + Af(i,j,k) + Ab(i,j,k))/RELAX;

                alphaw(i,j,k) = 0;
                if (Fw(i,j,k) > 0)
                    alphaw(i,j,k) = 1;
                end

                alphae(i,j,k) = 0;
                if (Fe(i,j,k) > 0)
                    alphae(i,j,k) = 1;
                end

                alphan(i,j,k) = 0;
                if (Fn(i,j,k) > 0)
                    alphan(i,j,k) = 1;
                end

                alphas(i,j,k) = 0;
                if (Fs(i,j,k) > 0)
                    alphas(i,j,k) = 1;
                end

                alphab(i,j,k) = 0;
                if (Fb(i,j,k) > 0)
                    alphab(i,j,k) = 1;
                end

                alphaf(i,j,k) = 0;
                if (Ff(i,j,k) > 0)
                    alphaf(i,j,k) = 1;
                end


                SUDC(i,j,k) = ((1/2)*Fe(i,j,k) * ((1 - alphae(i,j,k)) * PSI(REM(i,j,k)) - alphae(i,j,k) * PSI(REP(i,j,k))) * (U(i+1,j,k) - U(i,j,k)))+ ((1/2)*Fw(i,j,k) * (alphaw(i,j,k) * PSI(RWP(i,j,k)) - (1 - alphaw(i,j,k)) * PSI(RWM(i,j,k))) * (U(i,j,k) - U(i-1,j,k)))+ ((1/2)*Fn(i,j,k) * ((1 - alphan(i,j,k)) * PSI(RNM(i,j,k)) - alphan(i,j,k) * PSI(RNP(i,j,k))) * (U(i,j,k+1) - U(i,j,k)))+ ((1/2)*Fs(i,j,k) * (alphas(i,j,k) * PSI(RSP(i,j,k)) - (1 - alphas(i,j,k)) * PSI(RSM(i,j,k))) * (U(i,j,k) - U(i,j,k-1))) + ((1/2)*Ff(i,j,k) * ((1 - alphaf(i,j,k)) * PSI(RFM(i,j,k)) - alphaf(i,j,k) * PSI(RFP(i,j,k))) * (U(i,j+1,k) - U(i,j,k)) + ((1/2)*Fb(i,j,k) * (alphab(i,j,k) * PSI(RBP(i,j,k)) - (1 - alphab(i,j,k)) * PSI(RBM(i,j,k))) * (U(i,j,k) - U(i,j-1,k))));
                CON(i,j,k) = CON(i,j,k)*VOL + REL * AP(i,j,k) * T(i,j,k) + SUDC(i,j,k);
                
            end
        end
    end

                   for k = 2:N2
                        for j = 2:M2
                            for i = 2:L2
                                COFT(i,j,k,1) = CON(i,j,k);
                                COFT(i,j,k,2) = Ae(i,j,k);
                                COFT(i,j,k,3) = Aw(i,j,k);
                                COFT(i,j,k,4) = Af(i,j,k);
                                COFT(i,j,k,5) = Ab(i,j,k);
                                COFT(i,j,k,6) = As(i,j,k);
                                COFT(i,j,k,7) = An(i,j,k);
                                COFT(i,j,k,8) = AP(i,j,k);
                            end
                        end
                    end


    [T] = SOLVE(T, IST, JST, KST, AP, Af, Ab, Ae, Aw, An, As, CON, ITER); %%%%%%%%%%CALL SOLVE

    for k = 1:N1
        for i = 1:L1
            for j = 1:M1
                RHO(i,j,k) = RHO(i,j,k)/CP(i,j,k);
            end
        end
    end

    %----------------------T-------------------------------

NF = 5;
IST = 2;
JST = 2;
KST = 2;

Res = 0;
R_sum = 0;

for i = 2:L2
    for j = 2:M2
        for k = 2:N2
            CON(i,j,k) = COFT(i,j,k,1);
            Ae(i,j,k) = COFT(i,j,k,2);
            Aw(i,j,k) = COFT(i,j,k,3);
            Af(i,j,k) = COFT(i,j,k,4);
            Ab(i,j,k) = COFT(i,j,k,5);
            As(i,j,k) = COFT(i,j,k,6);
            An(i,j,k) = COFT(i,j,k,7);
            AP(i,j,k) = COFT(i,j,k,8);
        end
    end
end

for k = KST:N2
    for j = JST:M2
        for i = IST:L2

            Res = (AP(i,j,k)*T(i,j,k) - (Ae(i,j,k) * T(i+1,j,k) + Aw(i,j,k) * T(i-1,j,k)...
                + Ab(i,j,k) * T(i,j+1,k) + Af(i,j,k) * T(i,j-1,k)...
                + An(i,j,k) * T(i,j,k+1) + As(i,j,k) * T(i,j,k-1) + CON(i,j,k)) ) ^ 2; 

            R_sum = R_sum + Res;
        end
    end
end

Res_T = sqrt((R_sum/((L2-IST+1) * (M2 - JST + 1) * (N2 - KST + 1))));

Res_max = max([Res_U Res_V Res_W Res_T]);

if ( (mod(ITER,10)==0)||(ITER==LAST) )

     fprintf('%d   %e   %e   %e   %e   %e  %e  %e  %e  %e  %e \n', ITER, TIME, DT, Res_U, Res_V, Res_W, Res_T, U(6,6,6), V(6,6,6), W(6,6,6), T(6,6,6));

fp3 = fopen('residual_history.dat' , 'w' ) ;
fprintf (fp3, 'TITLE = "Residual history"\n' ) ;
fprintf (fp3, ' variables="Iterations", "Time taken", "Time step", "Residual_U", "Residual_V", "Residual_W", "Residual_T"\n');
% fprintf (fp1, ' zone T="n=%d"\n' , Nc) ;
fprintf (fp3, ' I= %d J= %d K=%d\n' ,202, 92, 92) ;
 

fprintf (fp3, '%e, %e, %e, %e, %e, %e, %e\n' , ITER , TIME, DT, Res_U, Res_V, Res_W, Res_T);


 
end
% 
    TIME = TIME + DT;
    ITER = ITER + 1;
%%%%%%%%%%%%%%% Check Iterative Convergence %%%%%%%%%%
% disp(ITER)
% disp(TIME)

fp1 = fopen('velocity1.dat' , 'w' ) ;
fprintf (fp1, 'TITLE = "Field Data"\n' ) ;
fprintf (fp1, ' variables= "x" (m), "y (m)", "z(m)", "u(m/s)", "v(m/s)", "w(m/s)"\n');
% fprintf (fp1, ' zone T="n=%d"\n' , Nc) ;
fprintf (fp1, ' zone T="n=%d", T = "n=%d", T = "n=%d"\n' ,22, 12, 12) ;
fprintf (fp1, ' I= %d J=%d K=%d\n' , 22, 12, 12) ;
for i =1: 22
    for j = 1:12
        for k = 1:12

            if (i == 1)
                uc(i,j,k) = U(2,j,k);
            elseif (i == 22)
                uc(i,j,k) = U(22,j,k);
            else
                uc(i,j,k) = 0.5*(U(i,j,k) + U(i+1,j,k));
            end

             if (j == 1)
                vc(i,j,k) = V(i,2,k);
            elseif (j == 12)
                vc(i,j,k) = V(i,12,k);
            else
                vc(i,j,k) = 0.5*(V(i,j,k) + V(i,j+1,k));
             end

             if (k == 1)
                wc(i,j,k) = W(i,j,2);
            elseif (k == 12)
                wc(i,j,k) = W(i,j,12);
            else
                wc(i,j,k) = 0.5*(W(i,j,k) + W(i,j,k+1));
             end


    fprintf (fp1, '%e, %e, %e, %e, %e, %e\n' , i , j, k,  uc(i,j,k) , vc(i,j,k), wc(i,j,k) ) ;

        end
    end
end
fclose(fp1);

fp2 = fopen('temperature.dat' , 'w' ) ;
fprintf (fp2, 'TITLE = "Field temperature"\n' ) ;
fprintf (fp2, ' variables="x (m)", "y(m)", "z(m)", "T(k)"\n');
fprintf (fp2, ' zone T="n=%d", T = "n=%d" , T = "n=%d"\n' , 22, 12, 12) ;
fprintf (fp2, ' I= %d J = %d K=%d\n' ,22, 12, 12) ;
for i =1: 22
    for j = 1:12
        for k = 1:12

            fprintf (fp2, '%e, %e, %e, %e\n' , i , j, k, T(i,j,k)) ;

        end
    end
end
fclose(fp2);



fp4 = fopen('pressure.dat' , 'w' ) ;
fprintf (fp4, 'TITLE = "Field pressure"\n' ) ;
fprintf (fp4, ' variables="x (m)", "y(m)", "z(m)", "P(N/m^2)"\n');
% fprintf (fp1, ' zone T="n=%d"\n' , Nc) ;
fprintf (fp4, ' zone T="n=%d", T = "n=%d", T = "n=%d"\n' , 22, 12, 12) ;
fprintf (fp4, ' I= %d J=%d K=%d\n' ,22, 12, 12) ;
for i =1: 22
    for j = 1:12
        for k = 1:12

            fprintf (fp4, '%e, %e, %e, %e\n' , i , j, k, P(i,j,k)) ;

        end
    end
end
fclose(fp4);

T_y = zeros(L1,N1);

for i = 1:L1
    for k = 1:N1
        T_y (i,k) = T(i,1,k);
    end
end
% 
% [X,Z] = meshgrid(XP,ZP);

u_y = zeros(L1,N1);

for i = 1:L1
    for k = 1:N1
        u_y (i,k) = uc(i,1,k);
    end
end

v_y = zeros(L1,N1);

for i = 1:L1
    for k = 1:N1
        v_y (i,k) = vc(i,1,k);
    end
end

w_y = zeros(L1,N1);

for i = 1:L1
    for k = 1:N1
        w_y (i,k) = wc(i,1,k);
    end
end

T_rate(ITER) = ((T_y(16,12) - T_y(16,11)) / ( (1.2/11) * (10^-3) ) ) * 0.1;

% % fp2 = fopen('temperature.dat' , 'w' ) ;
% % fprintf (fp2, 'TITLE = "Field temperature"\n' ) ;
% % fprintf (fp2, ' variables="x (m)", "z(m)", "T(k)"\n');
% % fprintf (fp2, ' zone T="n=%d", T = "n=%d"\n' , 22, 12) ;
% % fprintf (fp2, ' I= %d K=%d\n' ,22, 12) ;
% % for i =1: 22
% %     for j = 1:12
% %         for k = 1:12
% % 
% %             fprintf (fp2, '%e, %e, %e\n' , i , k, T_y(i,k)) ;
% % 
% %         end
% %     end
% % end
% % fclose(fp2);
% % 
% % fp1 = fopen('velocity.dat' , 'w' ) ;
% % fprintf (fp1, 'TITLE = "Field Data"\n' ) ;
% % fprintf (fp1, ' variables= "x" (m), "z(m)", "u(m/s)", "w(m/s)"\n');
% % fprintf (fp1, ' zone T="n=%d"\n' , Nc) ;
% % fprintf (fp1, ' zone T="n=%d",T = "n=%d"\n' ,22, 12) ;
% % fprintf (fp1, ' I= %d K=%d\n' , 22, 12) ;
% % for i = 1:22
% %     for k = 1:12
% %         fprintf(fp1, '%e, %e, %e, %e\n', i, k, u_y(i,k), w_y(i,k));
% %     end
% % end
% % fclose(fp1);
% % 
% % fp4 = fopen('pressure.dat' , 'w' ) ;
% % fprintf (fp4, 'TITLE = "Field pressure"\n' ) ;
% % fprintf (fp4, ' variables="x (m)", "y(m)", "z(m)", "P(N/m^2)"\n');
% % fprintf (fp1, ' zone T="n=%d"\n' , Nc) ;
% % fprintf (fp4, ' zone T="n=%d", T = "n=%d", T = "n=%d"\n' , 22, 12, 12) ;
% % fprintf (fp4, ' I= %d J=%d K=%d\n' ,22, 12, 12) ;
% % for i =1: 22
% %     for j = 1:12
% %         for k = 1:12
% % 
% %             fprintf (fp4, '%e, %e, %e, %e\n' , i , j, k, P(i,j,k)) ;
% % 
% %         end
% %     end
% % end
% % fclose(fp4);

% if (ITER == 350)

% figure(1), clf(1)
% % surf(X,Z,T_y)
% % contourf(T_y)
% contourf(T_y, 'LineStyle', 'none')
% view(0,90), xlabel('Z'), ylabel('X')
% c = colorbar;
% grid off
% title ('TEMPERATURE Y = 1 plane')
% colormap jet
% shading interp
% c.Label.String = 'Temperature(K)';
% 
% 
% 
% % figure(2)
% % plot(T(:,:,12))
% % xlabel('X coordinate')
% % ylabel('Temperature at different Y')
% % title('Temperature at Z = 12')
% % 
% figure(3), clf(3)
% % surf(X,Z,T_y) 
% % contourf(T_y)
% contourf(u_y, 'LineStyle', 'none')
% view(0,90), xlabel('Z'), ylabel('X')
% c = colorbar;
% grid off
% title ('U VELOCITY Y = 1 plane')
% colormap jet
% shading interp
% c.Label.String = ' X-Velocity(m/s)';
% % % 
% % figure(4), clf(4)
% % % surf(X,Z,T_y)
% % % contourf(T_y)
% % contourf(v_y, 'LineStyle', 'none')
% % view(0,90), xlabel('X'), ylabel('Z')
% % c = colorbar;
% % grid off
% % title ('V VELOCITY Y = 1 plane')
% % colormap jet
% % shading interp
% % c.Label.String = ' Y-Velocity(m/s)';
% figure(8)
% plot(T_rate)
% xlabel('Number of iterations')
% ylabel('Cooling rate')
% title('Cooling rate at Y = 1 plane')
% 
% % 
% % 
% figure(4), clf(4)
% % surf(X,Z,T_y)
% % contourf(T_y)
% contourf(w_y, 'LineStyle', 'none')
% view(0,90), xlabel('Z'), ylabel('X')
% c = colorbar;
% grid off
% title ('W VELOCITY Y = 1 plane')
% colormap jet
% shading interp
% c.Label.String = ' Z-Velocity(m/s)';
% % 
% figure(5)
% plot(u_y(:,:))
% xlabel('X coordinate')
% ylabel('U velocity at different Y')
% title('X- velocity at Z = 12')
% % 
% figure(6)
% plot(w_y(:,:))
% xlabel('X coordinate')
% ylabel('W velocity at different Y')
% title('Z- velocity at Z = 12')

% 
% % figure(7), clf(5)
% % % surf(X,Z,T_y)
% % % contourf(T_y)
% % contourf(p_y, 'LineStyle', 'none')
% % view(0,90), xlabel('X'), ylabel('Z')
% % c = colorbar;
% % grid off
% % title ('W VELOCITY Y = 1 plane')
% % colormap jet
% % shading interp
% % c.Label.String = ' Z-Velocity(m/s)';
% 
% % end
% 
% 


 

end

% err_u=(U-U_old).^2;
% error_u=(sum(err_u( :, :, : ))/(L1 *M1* N1)) .^ 0.5; %L2 norm of U velocity
% 
% err_v=(V-V_old).^2;
% error_v=(sum(err_v( :, :, : ))/(L1 *M1* N1)) .^ 0.5; %L2 norm of V velocity
% 
% err_w=(W-W_old).^2;
% error_w=(sum(err_w( :, :, : ))/(L1 *M1* N1)) .^ 0.5; %L2 norm of W velocity
% 
% err_T=(T-T_old).^2;
% error_T=(sum(err_T( :, :, : ))/(L1 *M1* N1)) .^ 0.5; %L2 norm of T velocity
% 
% cond = max(error_u, error_v, error_w, error_T);
% 
% TIME = TIME + DT;
% 
%  fprintf('%d   %e   %e   %e   %e   %e  %e  %e  %e  %e  %e \n', ITER, TIME, DT, error_U, error_V, error_W, error_T, U(6,6,6), V(6,6,6), W(6,6,6), T(6,6,6));
% 
%  U_old = U;
%  V_old = V;
%  W_old = W;
%  T_old = T;
% 
% 
% 
 %end
 fclose(fp3);


end



% ------------Extrapolate the pressure-------------------

function [P] = Extrapolate(P)

NI = 200;
NJ = 90;
NK = 90;

L1 = 22;
L2 = L1 - 1;
L3 = L2 - 1;

M1 = 12;
M2 = M1 - 1;
M3 = M2 - 1;

N1 = 12;
N2 = N1 - 1;
N3 = N2 - 1;


[XDIF, XCV, XCVS, XCVI, XCVIP, YDIF, YCV, YCVS, YCVJ, YCVJP, ZDIF, ZCV, ZCVS, ZCVK, ZCVKP, AX, AY, AZ, FX, FXM, FY, FYM, FZ, FZM, XP, YP, ZP] = Initialize_grid();

for k = 2:N2
    for j = 2:M2
        P(1,j,k) = (P(2,j,k) * XCVS(3) - P(3,j,k) * XDIF(2))/XDIF(3);
        P(L1,j,k) = (P(L2,j,k) * XCVS(L2) - P(L3,j,k)*XDIF(L1))/XDIF(L2);
    end
end

for k = 2:N2
    for i = 2:L2
        P(i,1,k) = (P(i,2,k) * YCVS(3) - P(i,3,k) * YDIF(2))/YDIF(3);
        P(i,M1,k) = (P(i,M2,k) * YCVS(M2) - P(i,M3,k) * YDIF(M1))/YDIF(M2);
    end
end

for j = 2:M2
    for i = 2:L2
        P(i,j,1) = (P(i,j,2) * ZCVS(3) - P(i,k,3) * ZDIF(2))/ZDIF(3);
        P(i,k,M1) = (P(i,j,N2) * ZCVS(N2) - P(i,j,N3) * ZDIF(N1))/ZDIF(N2);
    end
end

P(1,1,1) = (P(2,1,1) + P(1,2,1) + P(1,1,2) - P(2,2,2))/2;
P(L1,1,1) = (P(L2,1,1) + P(L1,2,1) + P(L1,1,2) - P(L2,2,2))/2;          
P(1,M1,1) = (P(1,M2,1) + P(2,M1,1) + P(1,M1,2) - P(2,M2,2))/2;
P(1,1,N1) = (P(1,1,N2) + P(2,1,N1) + P(1,2,N1) - P(2,2,N2))/2;
P(L1,M1,1) = (P(L1,M1,2) + P(L2,M1,1) + P(L1,M2,1) - P(L2,M2,2))/2;
P(L1,1,N1) = (P(L1,1,N2) + P(L2,1,N1) + P(L1,2,N1) - P(L2,2,N2))/2;
P(1,M1,N1) = (P(1,M1,N2) + P(2,M1,N1) + P(1,M2,N1) -   P(2,M2,N2))/2;
P(L1,M1,N1) = (P(L1,M1,N2) + P(L2,M1,N1) + P(L1,M2,N1) - P(L2,M2,N2))/2;

PREF = P(1,1,1);

for k = 1:N1
    for j = 1:M1
        for i = 1:L1
            P(i,j,k) = P(i,j,k) - PREF;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F] = SOLVE(F, IST, JST, KST, AP, Af, Ab, Ae, Aw, An, As, CON, ITER)

iter1 = ITER;

LAST = 1000;

NI = 20;
NJ = 10;
NK = 10;

L1 = 22;
L2 = L1 - 1;
L3 = L2 - 1;

M1 = 12;
M2 = M1 - 1;
M3 = M2 - 1;

N1 = 12;
N2 = N1 - 1;
N3 = N2 - 1;


ISTF = IST - 1;
JSTF = JST - 1;
KSTF = KST - 1;
IT1 = L2 + IST;
IT2 = L3 + IST;
JT1 = M2 + JST;
JT2 = M3 + JST;
KT1 = N2 + KST;
KT2 = N3 + KST;



PT(ISTF) = 0;
QT(ISTF) = 0;

if ((iter1 < LAST - 650) || (iter1 > LAST - 650))

%------------West to East-----------------%

for i = IST : L2
    BL = 0;
    BLP = 0;
    BLM = 0;
    BLC = 0;

    for j = JST:M2
        for k = KST:N2
            BL = BL + AP(i,j,k);
            if j ~= M2
                BL = BL - Ab(i,j,k);
            end
             if j ~= JST
                 BL = BL - Af(i,j,k);
             end
             if k ~= N2
                 BL = BL - An(i,j,k);
             end
             if k ~= KST
                 BL = BL - As(i,j,k);
             end

             BLP = BLP + Ae(i,j,k);
             BLM = BLM + Aw(i,j,k);
             BLC = BLC + CON(i,j,k) + Ae(i,j,k) * F(i+1,j,k) + Aw(i,j,k) * F(i-1,j,k) + Ab(i,j,k) * F(i,j+1,k) + Af(i,j,k) * F(i,j-1,k) + An(i,j,k) * F(i,j,k+1) + As(i,j,k) * F(i,j,k-1) - AP(i,j,k) * F(i,j,k);

        end
    end

    DENOM = BL - PT(i-1) * BLM;
    if abs(DENOM/BL) < 1.0e-10
         
        DENOM = 1.0e25;
    end
    PT(i) = BLP/DENOM;
    QT(i) = (BLC + BLM*QT(i-1))/DENOM;
end

BL = 0;

for i1 = IST:L2
    i = IT1 - i1;
    BL = BL * PT(i) + QT(i);
    for j = JST:M2
        for k = KST:N2
            F(i,j,k) = F(i,j,k) + BL;
        end
    end
end

%-----------------Front to back------------------

PT(JSTF) = 0;
QT(JSTF) = 0;

for j = JST : M2
    BL = 0;
    BLP = 0;
    BLM = 0;
    BLC = 0;

    for i = IST:L2
        for k = KST:N2
            BL = BL + AP(i,j,k);
            if i ~= L2
                BL = BL - Ae(i,j,k);
            end
             if i ~= IST
                 BL = BL - Aw(i,j,k);
             end
             if k ~= N2
                 BL = BL - An(i,j,k);
             end
             if k ~= KST
                 BL = BL - As(i,j,k);
             end

             BLP = BLP + Ab(i,j,k);
             BLM = BLM + Af(i,j,k);
             BLC = BLC + CON(i,j,k) + Ae(i,j,k) * F(i+1,j,k) + Aw(i,j,k) * F(i-1,j,k) + Ab(i,j,k) * F(i,j+1,k) + Af(i,j,k) * F(i,j-1,k) + An(i,j,k) * F(i,j,k+1) + As(i,j,k) * F(i,j,k-1) - AP(i,j,k) * F(i,j,k);

        end
    end

    DENOM = BL - PT(j-1) * BLM;
    if abs(DENOM/BL) < 1.0e-10
         
        DENOM = 1.0e25;
    end
    PT(j) = BLP/DENOM;
    QT(j) = (BLC + BLM*QT(j-1))/DENOM;
end

BL = 0;

for j1 = JST:M2
    j = JT1 - j1;
    BL = BL * PT(j) + QT(j);
    for i = IST:L2
        for k = KST:N2
            F(i,j,k) = F(i,j,k) + BL;
        end
    end
end

%----------South to north-------------------

PT(KSTF) = 0;
QT(KSTF) = 0;

for k = KST : N2
    BL = 0;
    BLP = 0;
    BLM = 0;
    BLC = 0;

    for j = JST:M2
        for i = IST:L2
            BL = BL + AP(i,j,k);
            if i ~= L2
                BL = BL - Ae(i,j,k);
            end
             if i ~= IST
                 BL = BL - Aw(i,j,k);
             end
             if j ~= M2
                 BL = BL - Ab(i,j,k);
             end
             if j ~= JST
                 BL = BL - Af(i,j,k);
             end

             BLP = BLP + An(i,j,k);
             BLM = BLM + As(i,j,k);
             BLC = BLC + CON(i,j,k) + Ae(i,j,k) * F(i+1,j,k) + Aw(i,j,k) * F(i-1,j,k) + Ab(i,j,k) * F(i,j+1,k) + Af(i,j,k) * F(i,j-1,k) + An(i,j,k) * F(i,j,k+1) + As(i,j,k) * F(i,j,k-1) - AP(i,j,k) * F(i,j,k);

        end
    end

    DENOM = BL - PT(k-1) * BLM;
    if abs(DENOM/BL) < 1.0e-10
         
        DENOM = 1.0e25;
    end
    PT(k) = BLP/DENOM;
    QT(k) = (BLC + BLM*QT(k-1))/DENOM;
end

BL = 0;

for k1 = KST:N2
    k = KT1 - k1;
    BL = BL * PT(k) + QT(k);
    for i = IST:L2
        for j = JST:M2
            F(i,j,k) = F(i,j,k) + BL;
        end
    end
end

end
% --------------------------------------------------------------
% --------------ADI---------------------------------------------
% 
% ---------i direction from west to east-------------------

for k = KST:N2
    for j = JST:M2
        PT(ISTF) = 0;
        QT(ISTF) = F(ISTF,j,k);
        for i = IST:L2
            DENOM = AP(i,j,k) - PT(i-1) * Aw(i,j,k);
            PT(i) = Ae(i,j,k) / (DENOM + 1.0e-30);
            TEMP = CON(i,j,k) + Ab(i,j,k) * F(i,j+1,k) + Af(i,j,k) * F(i,j-1,k) + An(i,j,k) * F(i,j,k+1) + As(i,j,k) * F(i,j,k-1);
            QT(i) = (TEMP + Aw(i,j,k) * QT(i-1))/(DENOM + 1.0e-30);
        end
        
        for i1 = IST:L2
            i = IT1 - i1;
            F(i,j,k) = F(i+1,j,k) * PT(i) + QT(i);
        end

    end
end

%---------i direction from east to west--------------

for k1 = KST:N3
    k = KT2 - k1;
    for j1 = JST:M3
        j = JT2 - j1;
        PT(ISTF) = 0;
        QT(ISTF) = F(ISTF,j,k);
        
        for i = IST:L2
            DENOM = AP(i,j,k) - PT(i-1) * Aw(i,j,k);
            PT(i) = Ae(i,j,k)/(DENOM + 1.0e-30);
            TEMP = CON(i,j,k) + Ab(i,j,k) * F(i,j+1,k) + Af(i,j,k) * F(i,j-1,k) + An(i,j,k) * F(i,j,k+1) + As(i,j,k) * F(i,j,k-1);
            QT(i) = (TEMP + Aw(i,j,k) * QT(i-1))/(DENOM + 1.0e-30);
           
        end
            
        for i1 = IST:L2
            i = IT1 - i1;
            F(i,j,k) = F(i+1,j,k) * PT(i) + QT(i);
        end
    end
end

%----------J direction from front to back------------

for k = KST:N2
    for i = IST:L2
        PT(JSTF) = 0;
        QT(JSTF) = F(i,JSTF,k);
        for j = JST:M2
            DENOM = AP(i,j,k) - PT(j-1) * Af(i,j,k);
            PT(j) = Ab(i,j,k) / (DENOM + 1.0e-30);
            TEMP = CON(i,j,k) + Ae(i,j,k) * F(i+1,j,k) + Aw(i,j,k) * F(i-1,j,k) + An(i,j,k) * F(i,j,k+1) + As(i,j,k) * F(i,j,k-1);
            QT(j) = (TEMP + Af(i,j,k) * QT(j-1))/(DENOM + 1.0e-30);
        end

        for j1 = JST:M2
            j = JT1 - j1;
            F(i,j,k) = F(i,j+1,k) * PT(j) + QT(j);
        end

    end
end

%------------J direction from back to front------------

for k1 = KST:N3
    k = KT2 - k1;
    for i1 = IST:L3
        i = IT2 - i1;
        PT(JSTF) = 0;
        QT(JSTF) = F(i,JSTF,k);

        for j = JST:M2
            DENOM = AP(i,j,k) - PT(j-1) * Af(i,j,k);
            PT(j) = Ab(i,j,k)/(DENOM + 1.0e-30);
            TEMP = CON(i,j,k) + Ae(i,j,k) * F(i+1,j,k) + Aw(i,j,k) * F(i-1,j,k) + An(i,j,k) * F(i,j,k+1) + As(i,j,k) * F(i,j,k-1);
            QT(i) = (TEMP + Af(i,j,k) * QT(j-1))/(DENOM + 1.0e-30);
        end

        for j1 = JST:M2
            j = JT1 - j1;
            F(i,j,k) = F(i,j+1,k) * PT(j) + QT(j);
        end
    end
end

%----------K direction from south to north--------------

for j = JST:M2
    for i = IST:L2
        PT(KSTF) = 0;
        QT(KSTF) = F(i,j,KSTF);
        for k = KST:N2
            DENOM = AP(i,j,k) - PT(k-1) * As(i,j,k);
            PT(k) = An(i,j,k) / (DENOM + 1.0e-30);
            TEMP = CON(i,j,k) + Ae(i,j,k) * F(i+1,j,k) + Aw(i,j,k) * F(i-1,j,k) + Ab(i,j,k) * F(i,j+1,k) + Af(i,j,k) * F(i,j-1,k);
            QT(k) = (TEMP + As(i,j,k) * QT(k-1))/(DENOM + 1.0e-30);
        end

        for k1 = KST:N2
            k = KT1 - k1;
            F(i,j,k) = F(i,j,k+1) * PT(k) + QT(k);
        end

    end
end

%-----------K direction from north to south-------------

 for j1 = JST:M3
    j = JT2 - j1;
    for i1 = IST:L3
        i = IT2 - i1;
        PT(KSTF) = 0;
        QT(KSTF) = F(i,j,KSTF);

        for k = KST:N2
            DENOM = AP(i,j,k) - PT(k-1) * As(i,j,k);
            PT(k) = An(i,j,k)/(DENOM + 1.0e-30);
            TEMP = CON(i,j,k) + Ae(i,j,k) * F(i+1,j,k) + Aw(i,j,k) * F(i-1,j,k) + Ab(i,j,k) * F(i,j+1,k) + Af(i,j,k) * F(i,j-1,k);
            QT(k) = (TEMP + As(i,j,k) * QT(k-1))/(DENOM + 1.0e-30);
        end

        for k1 = KST:N2
            k = KT1 - k1;
            F(i,j,k) = F(i,j,k+1) * PT(k) + QT(k);
        end
    end
 end

end

            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [U, V, W, GAM, CON, AP, CP, RHO] = Setup(ITER, LAST, U, V, W, T, NF)
% 
% % function START
% 
% NI = 20;
% NJ = 10;
% NK = 10;
% 
% L1 = 22;
% L2 = L1 - 1;
% L3 = L2 - 1;
% 
% M1 = 12;
% M2 = M1 - 1;
% M3 = M2 - 1;
% 
% N1 = 12;
% N2 = N1 - 1;
% N3 = N2 - 1;
% 
% 
% [XDIF, XCV, XCVS, XCVI, XCVIP, YDIF, YCV, YCVS, YCVJ, YCVJP, ZDIF, ZCV, ZCVS, ZCVK, ZCVKP, AX, AY, AZ, FX, FXM, FY, FYM, FZ, FZM, XP, YP, ZP] = Initialize_grid();
% 
% RHO = zeros(L1,M1,N1);%Density matrix
% CON = zeros(L1,M1,N1); %Source 
% AP = zeros(L1,M1,N1);
% CP = zeros(L1,M1,N1);
% 
% 
% Gravity = 9.8;
% Ub = 0.1;
% AMUS = 1.0e4;
% AMUL = 5.0e-3;
% KS = 3.5;
% KL = 5.0;
% 
% DS = 4.8e-10;
% DL = 4.8e-9;
% 
% Latent = 4.5e5;
% BetaT = 7.5e-6;
% BetaS = 2.6e-6;
% 
% Tmelt = 1375;
% DeltaT = 200;
% KKK = 0.8;
% Tamb = 298.15;
% Tf = 298.15;
% Hconv = 10;
% Rlaser = 0.0005;
% DeltT = -1.0e-5;
% Qlaser = 35.186;
% AbCoff = 0.5;
% Boltz = 5.67e-8;
% EPS = 0.2;
% 
% XL = 0.0036;
% YL = 0.0012;
% ZL = 0.0012;
% 
% AMUPLUS = zeros(L1,M1,N1);
% KPLUS = zeros(L1,M1,N1);
% 
% MassS = zeros(L1,M1,N1);
% MassL = zeros(L1,M1,N1);
% 
% 
% 
% for k = 1:N1
%     for j = 1:M1
%         for i = 1: L1
%             MassS(i,j,k) = 1.0;
%             MassL(i,j,k) = 1 - MassS(i,j,k);
%         end
%     end
% end
% 
% for j = 1:M1
%     for i = 1:L1
% 
%         SN1(i,j) = 0;
%         SN2(i,j) = 0;
%         SN3(i,j) = 0;
%         WS(i,j) = 0;
% 
%     end
% end
% 
% for j = 1:M1
%     SSN1(j) = 0;
%     LN1(j) = 0;
% end
% 
% 
% % for k = 1:N1
% %     for j = 1:M1
% %         for i = 1:L1
% %             U(i,j,k) = 0;
% %             V(i,j,k) = 0;
% %             W(i,j,k) = 0;
% %             T(i,j,k) = Tamb;
% %         end
% %     end
% % end
% 
% %-------------------------END---------------------%
% 
% %----------START BOUNDARY CONDITIONS-----------------------------%
% 
% end

function [AMUPLUS, KPLUS,RHO, CP, MassL, MassS, U,V,W,T] = Boundary_Conditions (U,V,W,T,ITER,LAST)

NI = 20;
NJ = 10;
NK = 10;

L1 = 22;
L2 = L1 - 1;
L3 = L2 - 1;

M1 = 12;
M2 = M1 - 1;
M3 = M2 - 1;

N1 = 12;
N2 = N1 - 1;
N3 = N2 - 1;

RHO = zeros(L1,M1,N1);%Density matrix

CP = zeros(L1,M1,N1);

EPS = 0.2;





[XDIF, XCV, XCVS, XCVI, XCVIP, YDIF, YCV, YCVS, YCVJ, YCVJP, ZDIF, ZCV, ZCVS, ZCVK, ZCVKP, AX, AY, AZ, FX, FXM, FY, FYM, FZ, FZM, XP, YP, ZP] = Initialize_grid();

Tmelt = 1375;
DeltaT = 200;
AMUS = 1.0e4;
AMUL = 5.0e-3;
KS = 3.5;
KL = 5.0;

Hconv = 10;
AbCoff = 0.5;
Boltz = 5.67e-8;
Tf = 298.15;
Tamb = 298.15;

Rlaser = 0.0005;
Qlaser = 35;

XL = 0.0036;
YL = 0.0012;
ZL = 0.0012;

Ub = 0.1;

DeltT = -1.0e-5;


MassS = zeros(L1,M1,N1);
MassL = zeros(L1,M1,N1);

for k = 1:N1
    for j = 1:M1
        for i = 1: L1
            MassS(i,j,k) = 1.0;
            MassL(i,j,k) = 1 - MassS(i,j,k);
        end
    end
end

for j = 1:M1
    for i = 1:L1

        SN1(i,j) = 0;
        SN2(i,j) = 0;
        SN3(i,j) = 0;
        WS(i,j) = 0;

    end
end

for j = 1:M1
    SSN1(j) = 0;
    LN1(j) = 0;
end
% 
  RHO = (1990 * (1-0.2)) .* ones(L1,M1,N1);
 CP = (1211 * (1-0.2)) .* ones(L1,M1,N1);
% %  
if (ITER < LAST - 650) %%%%%%%%%%%%%%%%%%%%%%%%%%   1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   RHO = (1990 * (1-0.2)) .* ones(L1,M1,N1);
%  CP = (1211 * (1-0.2)) .* ones(L1,M1,N1);

for k = 1:N1
    for j = 1:M1
        for i = 1:L1

            if (T(i,j,k) > Tmelt)

                MassL(i,j,k) = 1;
                MassS(i,j,k) = 1 - MassL(i,j,k);

            elseif (T(i,j,k) < Tmelt - DeltaT)

                MassL(i,j,k) = 0;
                MassS(i,j,k) = 1 - MassL(i,j,k);

            else

                MassL(i,j,k) = ((Tmelt - T(i,j,k))/DeltaT);
                MassS(i,j,k) = 1 - MassL(i,j,k);

            end

            AMUPLUS(i,j,k) = AMUS * MassS(i,j,k) + AMUL*MassL(i,j,k);
            KPLUS(i,j,k) = KS * MassS(i,j,k) + KL * MassL(i,j,k);

        end
    end
end

for j = 1:M1
    for i = 1:L1
        QQ(i,j) = 0;
        HEFF(i,j) = Hconv + AbCoff * Boltz * (T(i,j,N1) + Tf) * (T(i,j,N1) ^ 2 + Tf^2);
        GZ(i,j) = KPLUS(i,j,N1)/ZDIF(N1);
    end
end

for j = 1:M1
    for i = 1:L1
        D2 = ((XP(i) - XL*(3/4))^2 + (YP(j))^2);
        if (D2 < (Rlaser ^ 2))
            QQ(i,j) = (Qlaser/(3.1416 * (Rlaser)^2)) * exp(-((XP(i) - XL*(3/4))^2 + (YP(j))^2)/(Rlaser)^2);
        else
            QQ(i,j) = 0;
        end
    end
end

% North Boundary  -------------------------------
for j = 1:M1
    for i = 1:L1
        DELT = DeltT * (T(i,j,N1) - T(i,j,N2))/ZDIF(N1);
        W(i,j,1) = 0;

        if U(i,j,N1) > 0

            U(i,j,N1) = DELT/AMUPLUS(i,j,N2) * ZDIF(N1) + U(i,j,N2);

        end

        if U(i,j,N1) < 0

            U(i,j,N1) = -DELT/AMUPLUS(i,j,N2) * ZDIF(N1) + U(i,j,N2);

        end


        if V(i,j,N1) > 0

            V(i,j,N1) = DELT/AMUPLUS(i,j,N2) * ZDIF(N1) + V(i,j,N2);

        end

        if V(i,j,N1) < 0

            V(i,j,N1) = -DELT/AMUPLUS(i,j,N2) * ZDIF(N1) + V(i,j,N2);

        end

        T(i,j,N1) = (HEFF(i,j) * Tf + GZ(i,j) * T(i,j,N2) + QQ(i,j))/(HEFF(i,j) + GZ(i,j));

    end
end

%  West and east boundaries  -----------------------------------

for k = 1:N1
    for j = 1:M1

        V(1,j,k) = 0;
        U(1,j,k) = 0;
        W(i,j,k) = 0;

        V(L1,j,k) = 0;
        U(L1,j,k) = 0;
        W(L1,j,k) = 0;

        T(1,j,k) = Tamb;
        T(L1,j,k) = Tamb;
    end
end

% South Boundary ----------------------------

for j = 1:M1
    for i = 1:L1

        V(i,j,1) = 0;
        U(i,j,1) = 0;
        W(i,j,1) = 0;
        T(i,j,1) = Tamb;
    end
end

%  Back and front boundary --------------------------------

for k = 1:N1
    for i = 1:L1
        V(i,1,k) = 0;
        U(i,1,k) = U(i,2,k);
        W(i,1,k) = W(i,2,k);

        V(i,M1,k) = 0;
        U(i,M1,k) = 0;
        W(i,M1,k) = 0;

       T(i,1,k) = T(i,2,k);
       T(i,M1,k) =  Tamb;
    end
end

% RHO = (1990 * (1-0.2)) .* ones(L1,M1,N1);
% CP = (1211 * (1-0.2)) .* ones(L1,M1,N1);


% CALL PRINT     
%ITER == LAST - 100

elseif   (ITER == LAST - 650)

for j = 1:M1
    for i = 1:L1
        DELTAS1(i,j) = 0;
        DELTAS2(i,j) = 0;
        DELTAS3(i,j) = 0;
    end
end 

for j = 1:M1
    for i = 1:L1
        for k = 1:N1

            if (T(i,j,k) <= (Tmelt - DeltaT * EPS/(1-EPS)) && T(i,j,k) > (Tmelt - DeltaT))
                DELTAS1(i,j) = (1-EPS) * MassL(i,j,k) + DELTAS1(i,j);
                DELTAS3(i,j) = DELTAS3(i,j) + 1 - (1-EPS) * MassL(i,j,k);
            end

            if (T(i,j,k) > (Tmelt - DeltaT * EPS/(1 - EPS)) )
                DELTAS1(i,j) = EPS + DELTAS1(i,j);
                DELTAS2(i,j) = 1 - EPS + DELTAS2(i,j);
            end

        end
    end
end

for j = 1:M1
    for i = 1:L1
        SN1(i,j) = ceil(DELTAS1(i,j) - 1.0e-20);
        SN2(i,j) = ceil(DELTAS2(i,j) + DELTAS1(i,j) - 1.0e-20);
        SN3(i,j) = ceil(DELTAS3(i,j) + DELTAS2(i,j) + DELTAS1(i,j) - 1.0e-20);
    end
end

for j = 1:M1
    for i = 1:L1
        SSN1(j) = max(SSN1(j), SN1(i,j));
    end
end

for j = 1:M1
    for i = 1:L1
        if (SN1(i,j) == SSN1(j))

            LN1(j) = i;
        end
    end
end

for j = 1:M1
    for i = 1:L1
        if (i <= LN1(j))
            SN1(i,j) = SSN1(j);
            SN2(i,j) = SN2(LN1(j), j);
            SN3(i,j) = SN3(LN1(j), j);
        end
    end
end


for k = 1:N1
    for j = 1:M1
        for i = 1:L1
            if (T(i,j,k) > Tmelt)
                MassL(i,j,k) = 1;
                MassS(i,j,k) = 1 - MassL(i,j,k);
            elseif (T(i,j,k) < (Tmelt - DeltaT))
                MassL(i,j,k) = 0;
                MassS(i,j,k) = 1 - MassL(i,j,k);
            else
                MassL(i,j,k) = (Tmelt - T(i,j,k))/DeltaT;
                MassS(i,j,k) = 1 - MassL(i,j,k);
            end

            AMUPLUS(i,j,k) = AMUS * MassS(i,j,k) + AMUL * MassL(i,j,k);
            KPLUS(i,j,k) = KS * MassS(i,j,k) + KL * MassL(i,j,k);
        end
    end
end

for j = 1:M1
    for i = 1:L1
        QQ(i,j) = 0;
        HEFF(i,j) = Hconv + AbCoff * Boltz * (T(i,j,N1 - SN1(i,j)) + Tf) * (T(i,j,N1-SN1(i,j))^2 + Tf^2);
        GZ(i,j) = KPLUS(i,j,N1 - SN1(i,j))/ZDIF(N1);
    end
end

for j = 1:M1
    for i = 1:L1
        D2 = ((XP(i) - XL*(3/4))^2 + (YP(j))^2);
        if (D2 < (Rlaser ^ 2))
            QQ(i,j) = (Qlaser/(3.1416 * (Rlaser)^2)) * exp(-((XP(i) - XL*(3/4))^2 + (YP(j))^2)/(Rlaser)^2);
        else
            QQ(i,j) = 0;
        end

    end
end

% North Boundary  -------------------------------
for j = 1:M1
    for i = 1:L1
        if (i>1 && i<L1)
            WS(i,j) = 0 - EPS * Ub * (SN1(i+1,j) - SN1(i-1,j))*ZDIF(N1)/XCV(i);
        DELT = DeltT * (T(i,j,N1 - SN1(i,j)) - T(i,j,N2 - SN1(i,j)))/ZDIF(N1);
        end
        W(i,j,N1) = WS(i,j);

        if U(i,j,N1) > 0

            U(i,j,N1) = DELT/AMUPLUS(i,j,N2) * ZDIF(N1) + U(i,j,N2);

        end

        if U(i,j,N1) < 0

            U(i,j,N1) = -DELT/AMUPLUS(i,j,N2) * ZDIF(N1) + U(i,j,N2);

        end


        if V(i,j,N1) > 0

            V(i,j,N1) = DELT/AMUPLUS(i,j,N2) * ZDIF(N1) + V(i,j,N2);

        end

        if V(i,j,N1) < 0

            V(i,j,N1) = -DELT/AMUPLUS(i,j,N2) * ZDIF(N1) + V(i,j,N2);

        end

        T(i,j,N1) = (HEFF(i,j) * Tf + GZ(i,j) * T(i,j,N2) + QQ(i,j))/(HEFF(i,j) + GZ(i,j));

    end
end


for k=1:N1
    for j=1:M1
        for i=1:L1
                                                          
                  if(k >= (N1 - SN1(i,j))) % Set a few parameters 
                     AMUPLUS(i,j,k)=1.8e-5;
                     RHO(i,j,k)=1.225;
                     KPLUS(i,j,k)=0.0001;

                     

                     CP(i,j,k) = 718;
                     
                     MassL(i,j,k)=1;
                     MassS(i,j,k)=0;

                  end

                  if (k < (N1 - SN1(i,j)) && k >= (N1 - SN2(i,j)))
                     AMUPLUS(i,j,k) = AMUL;
                     RHO(i,j,k) = 1990;
                     KPLUS(i,j,k) = KL;
                     
                     CP(i,j,k) = 1211;
                     
                     MassL(i,j,k) = 1.0;
                     MassS(i,j,k) = 0;
                  end
                  
                  if (k < (N1-SN2(i,j)) && k >= (N1-SN3(i,j)))
                     RHO(i,j,k) = 1990 * (1.0-EPS) * MassS(i,j,k) + 1990 * MassL(i,j,k) * EPS;
                  end

        end
    end
end 

for k = 1:N1
    for j = 1:M1
        for i = 1:L1

            if ( k > SN1(i,j))

                T(i,j,k-SN1(i,j)) = T(i,j,k);
                U(i,j,k-SN1(i,j)) = U(i,j,k);
                V(i,j,k-SN1(i,j)) = V(i,j,k);
                W(i,j,k-SN1(i,j)) = W(i,j,k);
            end

        end
    end
end


for k = 1:N1
    for j = 1:M1
        for i = 1:L1

            if (k > (N1 - SN1(i,j)) )
                T(i,j,k) = Tamb;
                U(i,j,k) = 0;
                V(i,j,k) = 0;
                W(i,j,k) = 0;
            end

        end
    end
end

%-------North Boundary------------   

for j = 1:M1
    for i = 1:L1


        T(i,j,N1) = T(i,j,N1-1);
        U(i,j,N1) = 0;
        V(i,j,N1) = 0;
        W(i,j,N1) = 0;
    end
end

%-------- West and east boundaries -------------------

for k = 1:N1
    for j = 1:M1

        V(1,j,k) = 0;
        U(1,j,k) = 0;
        W(1,j,k) = 0;

        V(L1,j,k) = 0;
        U(L1,j,k) = 0;
        W(L1,j,k) = 0;

        T(1,j,k) = Tamb;
        T(L1,j,k) = Tamb;

    end
end

% ----------- South boundary-------------

for j = 1:M1
    for i = 1:L1

        V(i,j,1) = 0;
        U(i,j,1) = 0;
        W(i,j,1) = 0;

        T(i,j,1) = Tamb;

    end
end

%----------Back and front boundary-------------

for k = 1:N1
    for i = 1:L1
        V(i,1,k) = 0;
        U(i,1,k) = U(i,2,k);
        W(i,1,k) = W(i,2,k);

        V(i,M1,k) = 0;
        U(i,M1,k) = 0;
        W(i,M1,k) = 0;

        T(i,1,k) = T(i,2,k);
        T(i,M1,k) = Tamb;
    end
end

%------Call Print--------------

else

 for k = 1:N1
    for j = 1:M1
        for i = 1:L1

            if (T(i,j,k) > Tmelt)

                MassL(i,j,k) = 1;
                MassS(i,j,k) = 1 - MassL(i,j,k);

            elseif (T(i,j,k) < Tmelt - DeltaT)

                MassL(i,j,k) = 0;
                MassS(i,j,k) = 1 - MassL(i,j,k);

            else

                MassL(i,j,k) = ((Tmelt - T(i,j,k))/DeltaT);
                MassS(i,j,k) = 1 - MassL(i,j,k);

            end

            AMUPLUS(i,j,k) = AMUS * MassS(i,j,k) + AMUL * MassL(i,j,k);
            KPLUS(i,j,k) = KS * MassS(i,j,k) + KL * MassL(i,j,k);

            if (k >= (N1 - SN1(i,j)))

                AMUPLUS(i,j,k) = 1.81e-5;
                RHO(i,j,k) = 1.225;
                KPLUS(i,j,k) = 0.0001;
                CP(i,j,k) = 718;

                MassL(i,j,k) = 1;
                MassS(i,j,k) = 0;

            end

            if (k < (N1 - SN1(i,j)) && k >= (N1 - SN2(i,j)) )

                AMUPLUS(i,j,k) = AMUL;
                RHO(i,j,k) = 1990;
                KPLUS(i,j,k) = KL;
                CP(i,j,k) = 1211;

                MassL(i,j,k) = 1;
                MassS(i,j,k) = 0;
            end

            if (k < (N1 - SN2(i,j)) && k >= (N1 - SN3(i,j)))

                RHO(i,j,k) = 1990 * (1-EPS) * MassS(i,j,k) + 1990 * MassL(i,j,k) * EPS;
            end

        end
    end
 end

 for j = 1:M1
     for i = 1:L1
         QQ(i,j) = 0;
         HEFF(i,j) = Hconv + AbCoff * Boltz * (T(i,j,N1 - SN1(i,j)) + Tf) * (T(i,j,N1-SN1(i,j))^2 + Tf^2);
         GZ(i,j) = KPLUS(i,j,N1 - SN1(i,j))/ZDIF(N1);
    end
 end

 for j = 1:M1
    for i = 1:L1
        D2 = ((XP(i) - XL*(3/4))^2 + (YP(j))^2);
        if (D2 < (Rlaser ^ 2))
            QQ(i,j) = (Qlaser/(3.1416 * (Rlaser)^2)) * exp(-((XP(i) - XL*(3/4))^2 + (YP(j))^2)/(Rlaser)^2);
        else
            QQ(i,j) = 0;
        end

    end
end


% North Boundary  -------------------------------
for j = 1:M1
    for i = 1:L1
        if (i>1 && i<L1)
            WS(i,j) = 0 - EPS * Ub * (SN1(i+1,j) - SN1(i-1,j))*ZDIF(N1)/XCV(i);
        DELT = DeltT * (T(i,j,N1 - SN1(i,j)) - T(i,j,N2 - SN1(i,j)))/ZDIF(N1);
        end
        W(i,j,N1 - SN1(i,j)) = WS(i,j);

        if U(i,j,N1 - SN1(i,j)) > 0

            U(i,j,N1 - SN1(i,j)) = DELT/AMUPLUS(i,j,N2 - SN1(i,j)) * ZDIF(N1) + U(i,j,N2 - SN1(i,j));

        end

        if U(i,j,N1 - SN1(i,j)) < 0

            U(i,j,N1 - SN1(i,j)) = -DELT/AMUPLUS(i,j,N2 - SN1(i,j)) * ZDIF(N1) + U(i,j,N2 - SN1(i,j));

        end


        if V(i,j,N1 - SN1(i,j)) > 0

            V(i,j,N1 - SN1(i,j)) = DELT/AMUPLUS(i,j,N2 - SN1(i,j)) * ZDIF(N1) + V(i,j,N2 - SN1(i,j));

        end

        if V(i,j,N1 - SN1(i,j)) < 0

            V(i,j,N1 - SN1(i,j)) = -DELT/AMUPLUS(i,j,N2 - SN1(i,j)) * ZDIF(N1) + V(i,j,N2 - SN1(i,j));

        end
%         W(i,j,N1 - SN1(i,j)) = WS(i,j);

        T(i,j,N1 - SN1(i,j)) = (HEFF(i,j) * Tf + GZ(i,j) * T(i,j,N2 - SN1(i,j)) + QQ(i,j))/(HEFF(i,j) + GZ(i,j));

    end
end

for k = 1:N1
    for j = 1:M1
        for i = 1:L1

            if (k > (N1-SN1(i,j)))

                T(i,j,k) = Tamb;
                U(i,j,k) = 0;
                V(i,j,k) = 0;
                W(i,j,k) = 0;

            end

        end
    end
end

for j = 1:M1
    for i = 1:L1

        T(i,j,N1) = T(i,j,N1-1);
        U(i,j,N1) = 0;
        V(i,j,N1) = 0;
        W(i,j,N1) = 0;

    end
end

for k = 1:N1
    for j = 1:M1

        V(1,j,k) = 0;
        U(1,j,k) = 0;
        W(1,j,k) = 0;

        V(L1,j,k) = 0;
        U(L1,j,k) = 0;
        W(L1,j,k) = 0;

        T(1,j,k) = Tamb;
        T(L1,j,k) =  Tamb;

    end
end

for j = 1:M1
    for i = 1:L1
        V(i,j,1) = 0;
        U(i,j,1) = 0;
        W(i,j,1) = 0;

        T(i,j,1) = Tamb;
    end
end

for k = 1:N1
    for i = 1:L1

        V(i,1,k) = 0;
        U(i,1,k) = U(i,2,k);
        W(i,1,k) = W(i,2,k);

        V(i,M1,k) = 0;
        U(i,M1,k) = 0;
        W(i,M1,k) = 0;

        T(i,1,k) = T(i,2,k);
        T(i,M1,k) = Tamb;
    end
end


end

end

function [GAM,CON,AP] = Diffusion_Source(U, V, W, T, NF,ITER,LAST)

NI = 20;
NJ = 10;
NK = 10;

L1 = 22;
L2 = L1 - 1;
L3 = L2 - 1;

M1 = 12;
M2 = M1 - 1;
M3 = M2 - 1;

N1 = 12;
N2 = N1 - 1;
N3 = N2 - 1;


CON = zeros(L1,M1,N1); %Source 
AP = zeros(L1,M1,N1);

GAM = zeros(L1,M1,N1);
 

[AMUPLUS, KPLUS,RHO, CP, MassL, MassS, ~,~,~,~] = Boundary_Conditions (U,V,W,T,ITER,LAST);

Gravity = 9.8;
Ub = 0.1;
BetaT = 7.5e-6;
BetaS = 2.6e-6;
Tamb = 298.15;
Latent = 4.5e5;


[XDIF, XCV, XCVS, XCVI, XCVIP, YDIF, YCV, YCVS, YCVJ, YCVJP, ZDIF, ZCV, ZCVS, ZCVK, ZCVKP, AX, AY, AZ, FX, FXM, FY, FYM, FZ, FZM, XP, YP, ZP] = Initialize_grid();


if (NF == 1 || NF == 2 || NF == 3)

    for k = 1:N1
        for j = 1:M1
            for i = 1:L1

                GAM(i,j,k) = AMUPLUS(i,j,k);

            end
        end
    end
end


if (NF == 4)
    GAM(1:L1,1:M1,1:N1) = 0;
end

if (NF == 5)

        for k = 1:N1
            for j = 1:M1
                for i = 1:L1
                    GAM(i,j,k) = KPLUS(i,j,k);
                end
            end
        end
end

if (NF == 1)
    for k = 2:N2
        for j = 2:M2
            for i = 3:L2

                CON(i,j,k) = RHO(i,j,k) * Ub * (U(i+1,j,k) - U(i,j,k))/XCV(i)...
                    + (RHO(i+1,j,k) - RHO(i,j,k)) * Ub * U(i,j,k) / XDIF(i);

                AP(i,j,k) = 0;
            end
        end
    end
end

if (NF == 2)
    for k = 2:N2
        for j = 3:M2
            for i = 2:L2

                CON(i,j,k) = RHO(i,j,k) * Ub * (V(i+1,j,k) - V(i,j,k))/XCV(i)...
                    + (RHO(i+1,j,k) - RHO(i,j,k)) * Ub * V(i,j,k) / XDIF(i);

                AP(i,j,k) = 0;    
            end
        end
    end
end

if (NF == 3)
    for k = 3:N2
        for j = 2:M2
            for i = 2:L2

                CON(i,j,k) = RHO(i,j,k) * Gravity * BetaT * (T(i,j,k) - Tamb) + ...
                    RHO(i,j,k) * Ub * (W(i+1,j,k) - W(i,j,k))/XCV(i) + ...
                    (RHO(i+1,j,k) - RHO(i,j,k)) * Ub * W(i,j,k)/XDIF(i);

                AP(i,j,k) = 0;

            end
        end
    end
end

if (NF == 5)
    for k = 2:N2
        for j = 2:M2
            for i = 2:L2

                CON(i,j,k) = 0 - (RHO(i+1,j,k) * Latent * MassL(i+1,j,k) * (U(i+1,j,k) - Ub)...
                    - RHO(i,j,k) * Latent * MassL(i,j,k) * (U(i,j,k) - Ub))/ XCV(i)...
                    - (RHO(i,j+1,k) * Latent * MassL(i,j+1,k) * V(i,j+1,k)...
                    - RHO(i,j,k) * Latent * MassL(i,j,k) * V(i,j,k))/YCV(j)...
                    - (RHO(i,j,k+1) * Latent * MassL(i,j,k+1) * W(i,j,k+1)...
                    - RHO(i,j,k) * Latent * MassL(i,j,k) * W(i,j,k))/ZCV(k)...
                    + RHO(i,j,k) * Ub * CP(i,j,k) * (T(i+1,j,k) - T(i,j,k))/XDIF(i)...
                    + (RHO(i+1,j,k) - RHO(i,j,k)) * Ub * CP(i,j,k) * T(i,j,k)/XDIF(i);

                if (ITER > LAST - 650)
                    CON(i,j,k) = 0;
                end

                AP(i,j,k) = 0;
            end
        end
    end
end

if (NF == 4)
    CON(1:L1,1:M1,1:N1) = 0;
end






end
    

% Flux limiter function
function [Psi] =  PSI(A)

Psi = (A + (A^2))/(1 + (A^2));

end
 
   














